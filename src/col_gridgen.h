//
//  col_gridgen.h
//  swept_volume
//
//  Created by Yiwen Ju on 12/4/24.
//

#ifndef col_gridgen_h
#define col_gridgen_h

#include "adaptive_column_grid.h"
#include "ref_crit.h"
#include "timer.h"

bool gridRefine(mtet::MTetMesh &grid,
                vertExtrude &vertexMap,
                tetExtrude &cell5Map,
                const std::function<std::pair<Scalar, Eigen::RowVector4d>(Eigen::RowVector4d)> func,
                const double threshold,
                const int max_splits,
                std::array<double, timer_amount>& profileTimer);

#endif /* col_gridgen_h */
