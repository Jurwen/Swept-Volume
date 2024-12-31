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
                  vertExtrude vertexMap,
                  tetExtrude cell5Map){
    std::vector<Eigen::RowVector4d> verts;
    verts.reserve(grid.get_num_vertices() * 256);
    std::vector<std::array<size_t, 5>> simps;
    int simps_reserved = grid.get_num_tets() * 256;
    std::vector<std::array<size_t, 4>> ulsimp, llsimp;
    ulsimp.reserve(grid.get_num_tets());
    llsimp.reserve(grid.get_num_tets());
    simps.reserve(simps_reserved);
    using IndexMap = ankerl::unordered_dense::map<uint64_t, size_t>;
    IndexMap ind4DMap;
    std::vector<size_t> vertHashHead;
    size_t vertHashIt = 0;
    size_t curSum = 0;
    vertHashHead.reserve(grid.get_num_vertices());
    grid.seq_foreach_vertex([&](VertexId vid, std::span<const Scalar, 3> data){
        vertexCol::vert4d_list timeStamp = vertexMap[value_of(vid)].timeStamp;
        ind4DMap[value_of(vid)] = vertHashIt;
        for (size_t i = 0; i < timeStamp.size(); i ++){
            verts.emplace_back(timeStamp[i].coord);
        }
        vertHashHead[vertHashIt] = curSum;
        vertHashIt++;
        curSum += timeStamp.size();
    });
    std::cout << grid.get_num_vertices() << " " << verts.size() << " ";
    int simpNum = 0;
    grid.seq_foreach_tet([&](TetId tid, [[maybe_unused]] std::span<const VertexId, 4> data) {
        std::span<VertexId, 4> vs = grid.get_tet(tid);
        simpCol::cell5_list cell5Col = cell5Map[vs].cell5Col;
        if (cell5Col.size() == 0){
            return;
        }
        std::array<size_t, 4> headList;
        headList[0] = vertHashHead[ind4DMap[value_of(vs[0])]];
        headList[1] = vertHashHead[ind4DMap[value_of(vs[1])]];
        headList[2] = vertHashHead[ind4DMap[value_of(vs[2])]];
        headList[3] = vertHashHead[ind4DMap[value_of(vs[3])]];
        size_t i = 0;
        std::array<int, 5> simpHash;
        std::array<size_t, 5> curSimp;
        while (i < cell5Col.size()){
            if (!cell5Map[vs].cell5_eval[i]){
                i++;
                continue;
            }
            simpNum ++;
            if (simpNum > simps_reserved){
                simps_reserved *= 2;
                simps.reserve(simps_reserved);
            }
            simpHash = cell5Col[i]->hash;
            curSimp[0] = headList[0] + (size_t)simpHash[0];
            curSimp[1] = headList[1] + (size_t)simpHash[1];
            curSimp[2] = headList[2] + (size_t)simpHash[2];
            curSimp[3] = headList[3] + (size_t)simpHash[3];
            curSimp[4] = curSimp[simpHash[4]] - 1;
            simps.emplace_back(curSimp);
            i++;
        }
        if (cell5Map[vs].cell5_eval.front()){
            ulsimp.emplace_back(headList);
        }
        if (cell5Map[vs].cell5_eval.back()){
            std::array<int, 5> llSimpHash = cell5Col[cell5Col.size() - 1]->hash;
            std::array<size_t, 4> llface;
            llface[0] = curSimp[0];
            llface[1] = curSimp[1];
            llface[2] = curSimp[2];
            llface[3] = curSimp[3];
            llsimp.emplace_back(llface);
        }
    });
    std::cout << grid.get_num_tets() << " " << simps.size() << std::endl;
    if (std::filesystem::exists(filename.c_str())){
        std::filesystem::remove(filename.c_str());
    }
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str(),std::ios::app);
    json jOut;
    jOut.push_back(json(verts));
    jOut.push_back(json(simps));
    jOut.push_back(json(ulsimp));
    jOut.push_back(json(llsimp));
    fout << jOut.dump(4, ' ', true, json::error_handler_t::replace) << std::endl;
    fout.close();
    return true;
}
