#include <stdio.h>

#include "NetworkCheckpoint.hpp"
#include "io/binary_io.hpp"
#include "io/file_io.hpp"

using namespace std;

void NetworkCheckpoint::save_ml_network()
{
  if (ml_networks.empty() || loglh() > ml_networks.best_score())
    best_models = models;
  ml_networks.push_back(loglh(), network);
}

void NetworkCheckpoint::save_bs_network()
{
  bs_networks.push_back(loglh(), network);
}

void NetworkCheckpointManager::write(const std::string& ckp_fname) const
{
  backup();

  BinaryFileStream fs(ckp_fname, std::ios::out);

  fs << _checkp;

  remove_backup();
}

bool NetworkCheckpointManager::read(const std::string& ckp_fname)
{
  if (sysutil_file_exists(ckp_fname))
  {
    try
    {
      BinaryFileStream fs(ckp_fname, std::ios::in);

      fs >> _checkp;

      return true;
    }
    catch (runtime_error& e)
    {
      _checkp.search_state = SearchState();
      LOG_DEBUG << "Error reading checkpoint: " << e.what() << endl;
      return false;
    }
  }
  else
    return false;
}

void NetworkCheckpointManager::remove()
{
  if (sysutil_file_exists(_ckp_fname))
    std::remove(_ckp_fname.c_str());
}

void NetworkCheckpointManager::backup() const
{
  if (sysutil_file_exists(_ckp_fname))
    std::rename(_ckp_fname.c_str(), backup_fname().c_str());
}

void NetworkCheckpointManager::remove_backup() const
{
  if (sysutil_file_exists(backup_fname()))
    std::remove(backup_fname().c_str());
}

SearchState& NetworkCheckpointManager::search_state()
{
  if (_active)
    return _checkp.search_state;
  else
  {
    _empty_search_state = SearchState();
    return _empty_search_state;
  }
};

void NetworkCheckpointManager::reset_search_state()
{
  ParallelContext::thread_barrier();

  if (ParallelContext::master_thread())
    _checkp.search_state = SearchState();

  ParallelContext::thread_barrier();
};

void NetworkCheckpointManager::save_ml_network()
{
  if (ParallelContext::master_thread())
  {
    _checkp.save_ml_network();
    if (_active)
      write();
  }
}

void NetworkCheckpointManager::save_bs_network()
{
  if (ParallelContext::master_thread())
  {
    _checkp.save_bs_network();
    if (_active)
      write();
  }
}

void NetworkCheckpointManager::update_and_write(const NetworkInfo& networkinfo)
{
  if (!_active)
    return;

  if (ParallelContext::master_thread())
    _updated_models.clear();

  ParallelContext::barrier();

  for (auto p: networkinfo.parts_master())
  {
    /* we will modify a global map -> define critical section */
    ParallelContext::UniqueLock lock;

    assign(_checkp.models.at(p), networkinfo, p);

    /* remember which models were updated but this rank ->
     * will be used later to collect them at the master */
    _updated_models.insert(p);
  }

  ParallelContext::barrier();

  if (ParallelContext::num_ranks() > 1)
    gather_model_params();

  if (ParallelContext::master())
  {
    assign_network(_checkp, networkinfo);
    write();
  }
}

void NetworkCheckpointManager::gather_model_params()
{
  /* send callback -> worker ranks */
  auto worker_cb = [this](void * buf, size_t buf_size) -> int
      {
        BinaryStream bs((char*) buf, buf_size);
        bs << _updated_models.size();
        for (auto p: _updated_models)
        {
          bs << p << _checkp.models.at(p);
        }
        return (int) bs.pos();
      };

  /* receive callback -> master rank */
  auto master_cb = [this](void * buf, size_t buf_size)
     {
       BinaryStream bs((char*) buf, buf_size);
       auto model_count = bs.get<size_t>();
       for (size_t m = 0; m < model_count; ++m)
       {
         size_t part_id;
         bs >> part_id;

         // read parameter estimates from binary stream
         bs >> _checkp.models[part_id];
       }
     };

  ParallelContext::mpi_gather_custom(worker_cb, master_cb);
}

BasicBinaryStream& operator<<(BasicBinaryStream& stream, const NetworkCheckpoint& ckp)
{
  stream << ckp.version;

  // NB: accumulated runtime from past runs + current elapsed time
  stream << ckp.elapsed_seconds + global_timer().elapsed_seconds();

  stream << ckp.search_state;

  stream << ckp.network.topology();

  stream << ckp.models.size();
  for (const auto& m: ckp.models)
    stream << m.first << m.second;

  stream << ckp.ml_networks;

  stream << ckp.bs_networks;

  return stream;
}

BasicBinaryStream& operator>>(BasicBinaryStream& stream, NetworkCheckpoint& ckp)
{
  stream >> ckp.version;

  if (ckp.version < CKP_MIN_SUPPORTED_VERSION)
  {
    throw runtime_error("Unsupported checkpoint file version!");
  }

  stream >> ckp.elapsed_seconds;

  stream >> ckp.search_state;

//  auto topol = fs.get<NetworkTopology>();
//
//  printf("READ topology size: %u\n", topol.size());

  ckp.network.topology(stream.get<NetworkTopology>());

  size_t num_models, part_id;
  stream >> num_models;
  assert(num_models == ckp.models.size());
  for (size_t m = 0; m < num_models; ++m)
  {
    stream >> part_id;
    stream >> ckp.models[part_id];
  }

  stream >> ckp.ml_networks;

  stream >> ckp.bs_networks;

  return stream;
}

void assign_network(NetworkCheckpoint& ckp, const NetworkInfo& networkinfo)
{
  ckp.network = networkinfo.network();
}

void assign_model(NetworkCheckpoint& ckp, const NetworkInfo& networkinfo, size_t index)
{
  assign(ckp.models.at(index), networkinfo, index);
}

void assign_models(NetworkCheckpoint& ckp, const NetworkInfo& networkinfo)
{
  for (auto p: networkinfo.parts_master())
    assign_model(ckp, networkinfo, p);
}

void assign_models(NetworkInfo& networkinfo, const NetworkCheckpoint& ckp)
{
  const pllmod_networkinfo_t& pll_networkinfo = networkinfo.pll_networkinfo();
  for (auto& m: ckp.models)
  {
    if (!pll_networkinfo.partitions[m.first])
      continue;

    networkinfo.model(m.first, m.second);
  }
}

void assign(NetworkCheckpoint& ckp, const NetworkInfo& networkinfo)
{
  assign_network(ckp, networkinfo);
  assign_models(ckp, networkinfo);
}

void assign(NetworkInfo& networkinfo, const NetworkCheckpoint& ckp)
{
  // TODO: it is currently not possible to change tree after pll_treeinfo has been created
  // this should be fixed while doing refactoring to change pll_unode_t -> pll_tree_t
  assert(0);
  networkinfo.network(ckp.network);

  assign_models(networkinfo, ckp);
}

