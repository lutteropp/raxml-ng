/*
 * CheckpointCommon.hpp
 *
 *  Created on: Apr 23, 2019
 *      Author: sarah
 */

#ifndef SRC_CHECKPOINTCOMMON_HPP_
#define SRC_CHECKPOINTCOMMON_HPP_


constexpr int CKP_VERSION = 1;
constexpr int CKP_MIN_SUPPORTED_VERSION = 1;

typedef std::unordered_map<size_t, Model> ModelMap;

enum class CheckpointStep
{
  start,
  brlenOpt,
  modOpt1,
  radiusDetect,
  modOpt2,
  fastSPR,
  modOpt3,
  slowSPR,
  modOpt4,
  finish
};

struct SearchState
{
  SearchState() : step(CheckpointStep::start), loglh(0.), iteration(0), fast_spr_radius(0) {}

  CheckpointStep step;
  double loglh;

  int iteration;
  spr_round_params spr_params;
  int fast_spr_radius;
};


#endif /* SRC_CHECKPOINTCOMMON_HPP_ */
