//
//  ref_crit.h
//  swept_volume
//
//  Created by Yiwen Ju on 12/10/24.
//

#ifndef ref_crit_h
#define ref_crit_h

#include <cmath>

#include "adaptive_column_grid.h"
#include "timer.h"

bool refineFt(
              const std::array<vertex4d, 5>& verts,
              const double threshold,
              bool& inside,
              bool& choice,
              bool& zeroX,
              std::array<double, timer_amount>& profileTimer);

bool refineCap(
               const std::array<vertex4d, 4> verts,
               const double threshold,
               bool& zeroX);

bool refine3D(std::array<Eigen::RowVector4d, 4> pts,
              Eigen::RowVector<double, 4> vals,
              std::array<Eigen::RowVector4d, 4> grads,
              bool caps,
              const double threshold);
#endif /* ref_crit_h */
