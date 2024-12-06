//
//  col_gridgen.h
//  swept_volume
//
//  Created by Yiwen Ju on 12/4/24.
//

#ifndef col_gridgen_h
#define col_gridgen_h

#include "adaptive_column_grid.h"

bool gridRefine(mtet::MTetMesh grid,
                vertexCol &timeMap,
                tetCol &cell5Map,
                const std::function<Eigen::RowVector4d(std::span<const Scalar, 4>)> func,
                const double threshold);

#endif /* col_gridgen_h */
