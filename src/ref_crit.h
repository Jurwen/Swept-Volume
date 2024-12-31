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

bool refineContour(
                   const std::array<vertex4d, 5> verts,
                   const double threshold,
                   bool& inside,
                   bool& choice,
                   bool& zeroX);

bool refineCap(
               const std::array<vertex4d, 4> verts,
               const double threshold,
               bool& zeroX);
#endif /* ref_crit_h */
