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

void convert_4d_grid(mtet::MTetMesh grid,
                     vertExtrude vertexMap,
                     tetExtrude cell5Map,
                     std::vector<std::array<double, 4>> &verts,
                     std::vector<std::array<size_t, 5>> &simps,
                     std::vector<std::array<size_t, 4>> &ulsimp,
                     std::vector<std::array<size_t, 4>> &llsimp,
                     std::vector<double> &values){
    int vert_num = grid.get_num_vertices();
    int tet_num = grid.get_num_tets();
    verts.reserve(vert_num * 256);
    values.reserve(vert_num * 256);
    int simps_reserved = tet_num * 256;
    ulsimp.reserve(tet_num);
    llsimp.reserve(tet_num);
    simps.reserve(simps_reserved);
    using IndexMap = ankerl::unordered_dense::map<uint64_t, size_t>;
    IndexMap ind4DMap;
    std::vector<size_t> vertHashHead;
    size_t vertHashIt = 0;
    size_t curSum = 0;
    vertHashHead.reserve(vert_num);
    grid.seq_foreach_vertex([&](VertexId vid, std::span<const Scalar, 3> data){
        vertexCol::vert4d_list vert4dList = vertexMap[value_of(vid)].vert4dList;
        ind4DMap[value_of(vid)] = vertHashIt;
        for (size_t i = 0; i < vert4dList.size(); i ++){
            Eigen::RowVector4d coord = vert4dList[i].coord;
            verts.emplace_back(std::array<double, 4>{coord[0], coord[1], coord[2], coord[3]});
            values.emplace_back(vert4dList[i].valGradList.second[3]);
        }
        vertHashHead[vertHashIt] = curSum;
        vertHashIt++;
        curSum += vert4dList.size();
    });
    std::cout << vert_num << " " << verts.size() << " ";
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
//            if (!cell5Map[vs].cell5_eval[i]){
//                i++;
//                continue;
//            }
            simpNum ++;
            if (simpNum > simps_reserved){
                simps_reserved *= 2;
                simps.reserve(simps_reserved);
            }
            simpHash = cell5Col[i].hash;
            curSimp[0] = headList[0] + (size_t)simpHash[0];
            curSimp[1] = headList[1] + (size_t)simpHash[1];
            curSimp[2] = headList[2] + (size_t)simpHash[2];
            curSimp[3] = headList[3] + (size_t)simpHash[3];
            curSimp[4] = curSimp[simpHash[4]] - 1;
            simps.emplace_back(curSimp);
            i++;
        }
        ulsimp.emplace_back(headList);
        std::array<int, 5> llSimpHash = cell5Col[cell5Col.size() - 1].hash;
        std::array<size_t, 4> llface;
        llface[0] = curSimp[0];
        llface[1] = curSimp[1];
        llface[2] = curSimp[2];
        llface[3] = curSimp[3];
        llsimp.emplace_back(llface);
    });
    std::cout << tet_num << " " << simps.size() << " ";
}

bool save_4d_grid(const std::string filename,
                  const std::vector<std::array<double, 4>> verts,
                  const std::vector<std::array<size_t, 5>> simps,
                  const std::vector<std::array<size_t, 4>> ulsimp,
                  const std::vector<std::array<size_t, 4>> llsimp)
{
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

bool save_surface_mesh(const std::string filename,
                       const std::vector<std::array<double, 3>> output_vertices,
                       const std::vector<std::array<size_t, 3>> output_triangles)
{
    if (std::filesystem::exists(filename.c_str())){
        std::filesystem::remove(filename.c_str());
    }
    using json = nlohmann::json;
    std::ofstream fout(filename.c_str(),std::ios::app);
    json jOut;
    jOut.push_back(json(output_vertices));
    jOut.push_back(json(output_triangles));
    fout << jOut.dump(4, ' ', true, json::error_handler_t::replace) << std::endl;
    fout.close();
    return true;
}
