//
//  io.cpp
//  adaptive_colume_grid
//
//  Created by Yiwen Ju on 12/3/24.
//
#include "io.h"

bool save_mesh_json(const std::string& filename,
                    const mtet::MTetMesh mesh)
{
    std::vector<std::array<double, 3>> vertices((int)mesh.get_num_vertices());
    std::vector<std::array<size_t, 4>> tets((int)mesh.get_num_tets());
    ankerl::unordered_dense::map<uint64_t, size_t> vertex_tag_map;
    vertex_tag_map.reserve(mesh.get_num_vertices());
    int counter = 0;
    mesh.seq_foreach_vertex([&](VertexId vid, std::span<const Scalar, 3> data){
        size_t vertex_tag = vertex_tag_map.size() + 1;
        vertex_tag_map[value_of(vid)] = vertex_tag;
        vertices[counter] = {data[0], data[1], data[2]};
        counter ++;
    });
    counter = 0;
    mesh.seq_foreach_tet([&](TetId, std::span<const VertexId, 4> data) {
        tets[counter] = {vertex_tag_map[value_of(data[0])] - 1, vertex_tag_map[value_of(data[1])] - 1, vertex_tag_map[value_of(data[2])] - 1, vertex_tag_map[value_of(data[3])] - 1};
        counter ++;
    });
    if (std::filesystem::exists(filename.c_str())){
        std::filesystem::remove(filename.c_str());
    }
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str(),std::ios::app);
    json jOut;
    jOut.push_back(json(vertices));
    jOut.push_back(json(tets));
    fout << jOut.dump(4, ' ', true, json::error_handler_t::replace) << std::endl;
    fout.close();
    return true;
}

bool save_4d_grid(const std::string& filename,
                  mtet::MTetMesh grid,
                  vertexCol timeMap,
                  tetCol cell5Map){
    std::vector<std::array<double, 4>> verts;
    verts.reserve(grid.get_num_vertices() * 256);
    std::vector<std::array<size_t, 5>> simps;
    int simps_reserved = grid.get_num_tets() * 256;
    simps.reserve(simps_reserved);
    using IndexMap = ankerl::unordered_dense::map<uint64_t, size_t>;
    IndexMap ind4DMap;
    std::vector<size_t> vertHashHead;
    size_t vertHashIt = 0;
    size_t curSum = 0;
    vertHashHead.reserve(grid.get_num_vertices());
    grid.seq_foreach_vertex([&](VertexId vid, std::span<const Scalar, 3> data){
        std::vector<vertex4d> timeStamp = timeMap[value_of(vid)];
        ind4DMap[value_of(vid)] = vertHashIt;
        for (size_t i = 0; i < timeStamp.size(); i ++){
            double time = (double)timeStamp[i].time / 1024;
            verts.push_back(std::array<double, 4>{data[0], data[1], data[2], time});
        }
        vertHashHead[vertHashIt] = curSum;
        vertHashIt++;
        curSum += timeStamp.size();
    });
    int simpNum = 0;
    std::cout << "start simp output" << std::endl;
    grid.seq_foreach_tet([&](TetId tid, [[maybe_unused]] std::span<const VertexId, 4> data) {
        std::span<VertexId, 4> vs = grid.get_tet(tid);
        std::vector<cell5> cell5Col = cell5Map[vs];
        std::array<size_t, 4> headList;
        headList[0] = vertHashHead[ind4DMap[value_of(vs[0])]];
        headList[1] = vertHashHead[ind4DMap[value_of(vs[1])]];
        headList[2] = vertHashHead[ind4DMap[value_of(vs[2])]];
        headList[3] = vertHashHead[ind4DMap[value_of(vs[3])]];
        for (size_t i = 0; i < cell5Col.size(); i++){
            simpNum ++;
            if (simpNum > simps_reserved){
                simps_reserved *= 2;
                simps.reserve(simps_reserved);
            }
            std::array<size_t, 5> curSimp;
            std::array<int, 5> simpHash = cell5Col[i].hash;
            curSimp[0] = headList[0] + (int)simpHash[0];
            curSimp[1] = headList[1] + (int)simpHash[1];
            curSimp[2] = headList[2] + (int)simpHash[2];
            curSimp[3] = headList[3] + (int)simpHash[3];
            std::cout << " " << simpHash[0]<< " " << simpHash[1]<< " " << simpHash[2]<< " " << simpHash[3]<< " " << simpHash[4] << std::endl;
            curSimp[4] = curSimp[simpHash[4]] - 1;
            simps.push_back(curSimp);
        }
    });
    if (std::filesystem::exists(filename.c_str())){
        std::filesystem::remove(filename.c_str());
    }
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str(),std::ios::app);
    json jOut;
    jOut.push_back(json(verts));
    jOut.push_back(json(simps));
    fout << jOut.dump(4, ' ', true, json::error_handler_t::replace) << std::endl;
    fout.close();
    return true;
}
