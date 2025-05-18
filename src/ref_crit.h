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

/// 4D (un-projected) distance check for the time derivative function for a given 4D 5-cell. The default threshold here is 0.01, and it will not change with user's input
/// @param[in] verts         An array of 5 4D vertices. This vertex4d data structure also stores this vertex value and gradient.
/// @param[in] threshold            Threshold for the maxmimum allowed distance error for this 5-cell not to be subdivided.
/// @param[out] inside          An tag that checks if this 5-cell is completely covered by the sweep. The check is conservative that it is only true if the bezier ordiantes of either top face or the bottom face of the 5-cell is completely inside(negative).
/// @param[out] choice         A tag that outputs the choice of subdivision for this 5-cell. It compares the distance error between bezier control points on the temporal edge with points on the top and bottom faces.
/// @param[out] zeroX           A tag that shows if this 5-cell contains the intersection of time derivative function and the sweep function.
/// @param[out] profileTimer            The time profiler. Details can be found in `timer.h`
/// @return         A boolean that decides if this 5-cell needs to be subdivided.
bool refineFt(
              const std::array<vertex4d, 5>& verts,
              const double threshold,
              bool& inside,
              bool& choice,
              bool& zeroX,
              std::array<double, timer_amount>& profileTimer);

/// 3D (projected) distance check for the sweep function for 3D tetrahedra lives in the 4D domain, i.e., tets with 4D vertices. These tets are the results after the first stage subdivision.
/// @param[in] pts         An array of the coordinates of 4 4D vertices.
/// @param[in] vals            A list of 4 function values of these 4D vertices.
/// @param[in] grads          An array of 4 directional derivatives of the sweep functions for these 4D vertices.
/// @param[in] threshold            Threshold for the maxmimum allowed projected distance error for this tetrahedra.
/// @param[out] zeroX           A tag that shows if this 5-cell contains the intersection of time derivative function and the sweep function.
/// @return         A boolean that decides if this tetrahedra needs to be subdivided.
bool refine3D(std::array<Eigen::RowVector4d, 4> pts,
              Eigen::RowVector<double, 4> vals,
              std::array<Eigen::RowVector4d, 4> grads,
              const double threshold);
#endif /* ref_crit_h */
