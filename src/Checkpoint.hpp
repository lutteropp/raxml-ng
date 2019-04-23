#ifndef RAXML_CHECKPOINT_HPP_
#define RAXML_CHECKPOINT_HPP_

#include "common.h"
#include "TreeInfo.hpp"
#include "io/binary_io.hpp"
#include "CheckpointCommon.hpp"

struct Checkpoint
{
  Checkpoint() : version(CKP_VERSION), elapsed_seconds(0.), search_state(), tree(), models() {}

  Checkpoint(const Checkpoint&) = delete;
  Checkpoint& operator=(const Checkpoint&) = delete;
  Checkpoint(Checkpoint&&) = default;
  Checkpoint& operator=(Checkpoint&&) = default;

  int version;

  double elapsed_seconds;

  SearchState search_state;

  Tree tree;
  ModelMap models;
  ModelMap best_models;  /* model parameters for the best-scoring ML tree */

  TreeCollection ml_trees;
  TreeCollection bs_trees;

  double loglh() const { return search_state.loglh; }

  void save_ml_tree();
  void save_bs_tree();
};

class CheckpointManager
{
public:
  CheckpointManager(const std::string& ckp_fname) : _active(true), _ckp_fname(ckp_fname) {}

  const Checkpoint& checkpoint() { return _checkp; }
  void checkpoint(Checkpoint&& ckp) { _checkp = std::move(ckp); }

  //TODO: this is not very elegant, but should do the job for now
  SearchState& search_state();
  void reset_search_state();

  void enable() { _active = true; }
  void disable() { _active = false; }

  void update_and_write(const TreeInfo& treeinfo);

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
  Checkpoint _checkp;
  IDSet _updated_models;
  SearchState _empty_search_state;

  void gather_model_params();
  std::string backup_fname() const { return _ckp_fname + ".bk"; }
};

BasicBinaryStream& operator<<(BasicBinaryStream& stream, const Checkpoint& ckp);
BasicBinaryStream& operator>>(BasicBinaryStream& stream, Checkpoint& ckp);

void assign_tree(Checkpoint& ckp, const TreeInfo& treeinfo);
void assign_models(Checkpoint& ckp, const TreeInfo& treeinfo);
void assign_models(TreeInfo& treeinfo, const Checkpoint& ckp);

void assign(Checkpoint& ckp, const TreeInfo& treeinfo);
void assign(TreeInfo& treeinfo, const Checkpoint& ckp);

#endif /* RAXML_CHECKPOINT_HPP_ */
