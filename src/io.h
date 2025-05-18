//
//  io.h
//  adaptive_colume_grid
//
//  Created by Yiwen Ju on 12/3/24.
//

#ifndef io_h
#define io_h


#include "adaptive_column_grid.h"

void convert_4d_grid_col(mtet::MTetMesh grid,
                         vertExtrude vertexMap,
                         std::vector<std::array<double, 3>> &verts,
                         std::vector<std::array<size_t, 4>> &simps,
                         std::vector<std::vector<double>> &time,
                         std::vector<std::vector<double>> &values);

void convert_4d_grid_mtetcol(mtet::MTetMesh grid,
                             vertExtrude vertexMap,
                             std::vector<double> &verts,
                             std::vector<uint32_t> &simps,
                             std::vector<std::vector<double>> &time,
                             std::vector<std::vector<double>> &values,
                             bool cyclic);

#endif /* io_h */
