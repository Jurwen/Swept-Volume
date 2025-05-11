//
//  io.h
//  adaptive_colume_grid
//
//  Created by Yiwen Ju on 12/3/24.
//

#ifndef io_h
#define io_h


#include "adaptive_column_grid.h"


bool save_mesh_json(const std::string& filename,
                    const mtet::MTetMesh mesh);

void convert_4d_grid(mtet::MTetMesh grid,
                     vertExtrude vertexMap,
                     tetExtrude cell5Map,
                     std::vector<std::array<double, 4>> &verts,
                     std::vector<std::array<size_t, 5>> &simps,
                     std::vector<std::array<size_t, 4>> &ulsimp,
                     std::vector<std::array<size_t, 4>> &llsimp,
                     std::vector<double> &values);

bool save_4d_grid(const std::string filename,
                  const std::vector<std::array<double, 4>> verts,
                  const std::vector<std::array<size_t, 5>> simps,
                  const std::vector<std::array<size_t, 4>> ulsimp,
                  const std::vector<std::array<size_t, 4>> llsimp);

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

bool save_surface_mesh(const std::string filename,
                       const std::vector<std::array<double, 3>> output_vertices,
                       const std::vector<std::array<size_t, 3>> output_triangles);

#endif /* io_h */
