//
//  io.cpp
//  adaptive_colume_grid
//
//  Created by Yiwen Ju on 12/3/24.
//
#include "io.h"


void convert_4d_grid_col(mtet::MTetMesh grid,
                         vertExtrude vertexMap,
                         std::vector<std::array<double, 3>> &verts,
                         std::vector<std::array<size_t, 4>> &simps,
                         std::vector<std::vector<double>> &time,
                         std::vector<std::vector<double>> &values){
    size_t vert_num = grid.get_num_vertices();
    size_t tet_num = grid.get_num_tets();
    verts.reserve(vert_num);
    simps.reserve(tet_num);
    time.reserve(vert_num);
    values.reserve(vert_num);
    size_t vertIt = 0;
    using IndexMap = ankerl::unordered_dense::map<uint64_t, size_t>;
    IndexMap ind4DMap;
    grid.seq_foreach_vertex([&](VertexId vid, std::span<const Scalar, 3> data){
        verts.emplace_back(std::array<double, 3>{data[0], data[1], data[2]});
        vertexCol::vert4d_list vert4dList = vertexMap[value_of(vid)].vert4dList;
        ind4DMap[value_of(vid)] = vertIt;
        values.emplace_back(std::vector<double>{});
        time.emplace_back(std::vector<double>{});
        values[vertIt].reserve(vert4dList.size());
        time[vertIt].reserve(vert4dList.size());
        for (size_t i = 0; i < vert4dList.size(); i ++){
            Eigen::RowVector4d coord = vert4dList[i].coord;
            values[vertIt].emplace_back(vert4dList[i].valGradList.second[3]);
            time[vertIt].emplace_back(vert4dList[i].coord(3));
        }
        vertIt++;
    });
    grid.seq_foreach_tet([&](TetId tid, [[maybe_unused]] std::span<const VertexId, 4> data) {
        std::span<VertexId, 4> vs = grid.get_tet(tid);
        simps.emplace_back(std::array<size_t, 4>{ind4DMap[value_of(vs[0])],
            ind4DMap[value_of(vs[1])],
            ind4DMap[value_of(vs[2])],
            ind4DMap[value_of(vs[3])]});
    });
}

void convert_4d_grid_mtetcol(mtet::MTetMesh grid,
                             vertExtrude vertexMap,
                             std::vector<double> &verts,
                             std::vector<uint32_t> &simps,
                             std::vector<std::vector<double>> &time,
                             std::vector<std::vector<double>> &values,
                             bool cyclic){
    size_t vert_num = grid.get_num_vertices();
    size_t tet_num = grid.get_num_tets();
    size_t tet4d_num = 0, vert4d_num = 0;
    verts.reserve(vert_num * 3);
    simps.reserve(tet_num * 4);
    time.reserve(vert_num);
    values.reserve(vert_num);
    size_t vertIt = 0;
    using IndexMap = ankerl::unordered_dense::map<uint64_t, size_t>;
    IndexMap ind4DMap;
    grid.seq_foreach_vertex([&](VertexId vid, std::span<const Scalar, 3> data){
        verts.emplace_back(static_cast<double>(data[0]));
        verts.emplace_back(static_cast<double>(data[1]));
        verts.emplace_back(static_cast<double>(data[2]));
        vertexCol::vert4d_list vert4dList = vertexMap[value_of(vid)].vert4dList;
        ind4DMap[value_of(vid)] = vertIt;
        values.emplace_back(std::vector<double>{});
        time.emplace_back(std::vector<double>{});
        values[vertIt].reserve(vert4dList.size());
        time[vertIt].reserve(vert4dList.size());
        vert4d_num += vert4dList.size();
        for (size_t i = 0; i < vert4dList.size(); i ++){
            Eigen::RowVector4d coord = vert4dList[i].coord;
            values[vertIt].emplace_back(vert4dList[i].valGradList.second[3]);
            time[vertIt].emplace_back(vert4dList[i].coord(3));
        }
        if (cyclic){
            values[vertIt].back() = values[vertIt].front();
        }
        vertIt++;
    });
    grid.seq_foreach_tet([&](TetId tid, [[maybe_unused]] std::span<const VertexId, 4> data) {
        std::span<VertexId, 4> vs = grid.get_tet(tid);
        simps.emplace_back(static_cast<uint32_t>(ind4DMap[value_of(vs[0])]));
        simps.emplace_back(static_cast<uint32_t>(ind4DMap[value_of(vs[1])]));
        simps.emplace_back(static_cast<uint32_t>(ind4DMap[value_of(vs[2])]));
        simps.emplace_back(static_cast<uint32_t>(ind4DMap[value_of(vs[3])]));
        tet4d_num += vertexMap[value_of(vs[0])].vert4dList.size();
        tet4d_num += vertexMap[value_of(vs[1])].vert4dList.size();
        tet4d_num += vertexMap[value_of(vs[2])].vert4dList.size();
        tet4d_num += vertexMap[value_of(vs[3])].vert4dList.size();
        tet4d_num -= 4;
    });
    std::cout << "4D Vertex Number: " << vert4d_num << " 4D Tetrahedra Number: " << tet4d_num << std::endl;
}
