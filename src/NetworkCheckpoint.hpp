#ifndef RAXML_NETWORK_CHECKPOINT_HPP_
#define RAXML_NETWORK_CHECKPOINT_HPP_

#include "common.h"
#include "NetworkInfo.hpp"
#include "io/binary_io.hpp"
#include "CheckpointCommon.hpp"

struct NetworkCheckpoint
{
	NetworkCheckpoint() : version(CKP_VERSION), elapsed_seconds(0.), search_state(), network(), models() {}

	NetworkCheckpoint(const NetworkCheckpoint&) = delete;
	NetworkCheckpoint& operator=(const NetworkCheckpoint&) = delete;
	NetworkCheckpoint(NetworkCheckpoint&&) = default;
	NetworkCheckpoint& operator=(NetworkCheckpoint&&) = default;

  int version;

  double elapsed_seconds;

  SearchState search_state;

  Network network;
  ModelMap models;
  ModelMap best_models;  /* model parameters for the best-scoring ML tree */

  NetworkCollection ml_networks;
  NetworkCollection bs_networks;

  double loglh() const { return search_state.loglh; }

  void save_ml_network();
  void save_bs_network();
};

class NetworkCheckpointManager
{
public:
	NetworkCheckpointManager(const std::string& ckp_fname) : _active(true), _ckp_fname(ckp_fname) {}

  const NetworkCheckpoint& checkpoint() { return _checkp; }
  void checkpoint(NetworkCheckpoint&& ckp) { _checkp = std::move(ckp); }

  //TODO: this is not very elegant, but should do the job for now
  SearchState& search_state();
  void reset_search_state();

  void enable() { _active = true; }
  void disable() { _active = false; }

  void update_and_write(const NetworkInfo& networkinfo);

  void save_ml_tree();
  void save_bs_tree();

  bool read() { return read(_ckp_fname); }
  bool read(const std::string& ckp_fname);
  void write() const { write(_ckp_fname); }
  void write(const std::string& ckp_fname) const;

  void remove();
  void backup() const;
  void remove_backup() const;

private:
  bool _active;
  std::string _ckp_fname;
  NetworkCheckpoint _checkp;
  IDSet _updated_models;
  SearchState _empty_search_state;

  void gather_model_params();
  std::string backup_fname() const { return _ckp_fname + ".bk"; }
};

BasicBinaryStream& operator<<(BasicBinaryStream& stream, const NetworkCheckpoint& ckp);
BasicBinaryStream& operator>>(BasicBinaryStream& stream, NetworkCheckpoint& ckp);

void assign_network(NetworkCheckpoint& ckp, const NetworkInfo& networkinfo);
void assign_models(NetworkCheckpoint& ckp, const NetworkInfo& networkinfo);
void assign_models(NetworkInfo& networkinfo, const NetworkCheckpoint& ckp);

void assign(NetworkCheckpoint& ckp, const NetworkInfo& networkinfo);
void assign(NetworkInfo& networkinfo, const NetworkCheckpoint& ckp);

#endif /* RAXML_NETWORK_CHECKPOINT_HPP_ */
