//
//  col_gridgen.h
//  swept_volume
//
//  Created by Yiwen Ju on 12/4/24.
//

#ifndef col_gridgen_h
#define col_gridgen_h
#include <nanothread/nanothread.h>
#include <iostream>
#include <mtetcol/simplicial_column.h>
#include "adaptive_column_grid.h"
#include "ref_crit.h"
#include "timer.h"

void initNewVert(vertexCol& newVert, std::span<const mtet::Scalar, 3> data, const std::function<std::pair<Scalar, Eigen::RowVector4d>(Eigen::RowVector4d)> func, const double threshold);

/// Create a column-wise 4D grid data structure. It contains a base 3D tetrahedra grid. On top of each vertex and tetrahedra, there sits a column of 4D data. For a 3D vertex, there is a list of 4D vertices that only differ in the fourth coordinate. For a 3D tet, there is a list of 4D 5-cell/simplex that two of its tetrahedra faces can be projected down to the same 3D tet, and three of the faces are projected down to 2D triangular faces.
/// @param[out] grid         The base 3D tetrahedra grid.
/// @param[out] vertexMap            This maps a 3D vertex to a 4D vertex column. Details of the data structure can be found in `adaptive_column_grid.h`.
/// @param[out] insideMap         This maps a 3D tet to an inside-ness tag. Details can be found `adaptive_column_grid.h`.
/// @param[in] func         The implicit function that represents the sweep. It takes in a 4D coordinate and outputs a Scalar of value and a size 4 vector of gradient.
/// @param[in] threshold            The threshold for the refinement critieria of each 4D simplex.
/// @param[in] max_splits           Max number of splits of the grid
/// @param[out] profileTimer            The time profiler. Details can be found in `timer.h`
bool gridRefine(mtet::MTetMesh &grid,
                vertExtrude &vertexMap,
                insidenessMap &insideMap,
                const std::function<std::pair<Scalar, Eigen::RowVector4d>(Eigen::RowVector4d)> func,
                const double threshold,
                const int max_splits,
                std::array<double, timer_amount>& profileTimer);

#endif /* col_gridgen_h */
