//
//  io.h
//  adaptive_colume_grid
//
//  Created by Yiwen Ju on 12/3/24.
//

#ifndef io_h
#define io_h

#include <ankerl/unordered_dense.h>
#include "adaptive_column_grid.h"

using namespace mtet;

bool save_mesh_json(const std::string& filename,
                    const mtet::MTetMesh mesh);

#endif /* io_h */
