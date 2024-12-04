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

bool save_4d_grid(const std::string& filename,
                  const mtet::MTetMesh grid,
                  const vertexCol timeMap,
                  const tetCol cell5Map);

#endif /* io_h */
