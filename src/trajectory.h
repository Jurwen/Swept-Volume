//
//  trajectory.h
//  adaptive_column_grid
//
//  Created by Yiwen Ju on 12/5/24.
//

#ifndef trajectory_h
#define trajectory_h

#include "adaptive_column_grid.h"

double sweptFunc();

Eigen::RowVector3d trajLine(double t, std::span<const Scalar, 3> coords){
    Eigen::RowVector3d translation = (1 - t) * Eigen::RowVector3d{1/3, 0., 0.} + t * Eigen::RowVector3d{0, 1/3, 0.} + Eigen::RowVector3d{1/3, 1/3 , 1/2};
    return Eigen::Map<const Eigen::RowVector3d>(coords.data()) - translation;
}



#endif /* trajectory_h */
