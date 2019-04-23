/*
 * SPRRoundParams.hpp
 *
 *  Created on: Apr 23, 2019
 *      Author: sarah
 */

#ifndef SRC_SPRROUNDPARAMS_HPP_
#define SRC_SPRROUNDPARAMS_HPP_

struct spr_round_params
{
  bool thorough;
  int radius_min;
  int radius_max;
  int ntopol_keep;
  double subtree_cutoff;
  cutoff_info_t cutoff_info;

  void reset_cutoff_info(double loglh)
  {
    cutoff_info.lh_dec_count = 0;
    cutoff_info.lh_dec_sum = 0.;
    cutoff_info.lh_cutoff = loglh / -1000.0;
  }
};


#endif /* SRC_SPRROUNDPARAMS_HPP_ */
