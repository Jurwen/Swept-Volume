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
#include <convex_hull_membership/contains.h>

bool refine4D(
//              const std::array<Eigen::RowVector4d, 5>& pts,
//              const std::array<double, 5>& vals,
//              const std::array<Eigen::RowVector4d, 5>& grads,
              const std::array<vertex4d, 5> verts,
              const double threshold,
              bool& inside);

bool refine4D_test(
              const std::array<Eigen::RowVector4d, 5>& pts,
              const std::array<double, 5>& vals,
              const std::array<Eigen::RowVector4d, 5>& grads,
              const double threshold,
                   bool& inside);


#endif /* ref_crit_h */
