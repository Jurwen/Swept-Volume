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
    const int zeros = 4;
    std::vector<size_t> ind_list;
    std::vector<std::vector<std::array<double, 5>>> vert_lists;
//    std::vector<std::array<std::vector<std::array<double, 4>>, 4>> tet4D_list;
    std::vector<std::array<std::vector<std::array<double, 5>>, 4>> tet4D_list_complex;
    auto get_sign = [&](double x) {
        return x > 0;
    };
    auto get_vert_lists = [&](vertexCol::vert4d_list vert4dList){
        std::vector<std::array<double, 5>> ret;
        ret.reserve(vert4dList.size());
        for (size_t i = 0; i < vert4dList.size(); i++){
            Eigen::RowVector4d coord = vert4dList[i].coord;
            ret.emplace_back(std::array<double, 5>{coord[0], coord[1], coord[2], coord[3], (double)vert4dList[i].inherit});
        }
        return ret;
    };
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
        bool sign = get_sign(vert4dList[0].valGradList.first);
        bool d_sign = get_sign(vert4dList[0].valGradList.second[3]);
        //        int counter = 0;
        for (size_t i = 0; i < vert4dList.size(); i ++){
            Eigen::RowVector4d coord = vert4dList[i].coord;
            verts.emplace_back(std::array<double, 4>{coord[0], coord[1], coord[2], coord[3]});
            values.emplace_back(vert4dList[i].valGradList.second[3]);
//            if (sign != get_sign(vert4dList[i].valGradList.first) && d_sign != get_sign(vert4dList[i].valGradList.second[3])){
//                std::vector<std::array<double, 5>> vert_col =get_vert_lists(vert4dList);
//                vert_lists.push_back(vert_col);
//                //                sign = get_sign(vert4dList[i].valGradList.first);
//                //                counter++;
//            }
//            sign = get_sign(vert4dList[i].valGradList.first);
//            d_sign = get_sign(vert4dList[i].valGradList.second[3]);
        }
        if (vert4dList[0].coord[0] > 0.425 && vert4dList[0].coord[1] > 0.35 && vert4dList[0].coord[2] > 0.2 && vert4dList[0].coord[0] < 0.525 && vert4dList[0].coord[1] < 0.45 && vert4dList[0].coord[2] < 0.3){
            std::vector<std::array<double, 5>> vert_col =get_vert_lists(vert4dList);
            vert_lists.push_back(vert_col);
            ind_list.push_back(curSum);
        }
        //        if (counter <= zeros && ind_list[counter - 1] == 0){
        //            ind_list[counter - 1] = curSum;
        //        }
        vertHashHead[vertHashIt] = curSum;
        vertHashIt++;
        curSum += vert4dList.size();
    });
//    for (size_t i = 0; i < zeros; i++){
//        std::cout << ind_list[i] << std::endl;
//    }
    std::cout << " vertex list size: "<< vert_lists.size() << std::endl;
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
//        if (headList[0] == ind_list[0] || headList[1] == ind_list[0] || headList[2] == ind_list[0] || headList[3] == ind_list[0]){
//            std::cout << "here" << std::endl;
//            vertexCol::vert4d_list ti = vertexMap[value_of(vs[0])].vert4dList;
//            vertexCol::vert4d_list tj = vertexMap[value_of(vs[1])].vert4dList;
//            vertexCol::vert4d_list tk = vertexMap[value_of(vs[2])].vert4dList;
//            vertexCol::vert4d_list tl = vertexMap[value_of(vs[3])].vert4dList;
//            std::array<std::vector<std::array<double, 4>>, 4> tet_4d_verts =
//            {get_vert_lists(ti),get_vert_lists(tj),get_vert_lists(tk),get_vert_lists(tl)};
//            tet4D_list.push_back(tet_4d_verts);
//        }
//        if (headList[0] == ind_list[3] || headList[1] == ind_list[3] || headList[2] == ind_list[3] || headList[3] == ind_list[3]){
        if (std::count(ind_list.begin(), ind_list.end(), headList[0]) > 0 &&
            std::count(ind_list.begin(), ind_list.end(), headList[1]) > 0 &&
            std::count(ind_list.begin(), ind_list.end(), headList[2]) > 0 &&
            std::count(ind_list.begin(), ind_list.end(), headList[3]) > 0){
//            std::cout << "here" << std::endl;
            vertexCol::vert4d_list ti = vertexMap[value_of(vs[0])].vert4dList;
            vertexCol::vert4d_list tj = vertexMap[value_of(vs[1])].vert4dList;
            vertexCol::vert4d_list tk = vertexMap[value_of(vs[2])].vert4dList;
            vertexCol::vert4d_list tl = vertexMap[value_of(vs[3])].vert4dList;
            std::array<std::vector<std::array<double, 5>>, 4> tet_4d_verts =
            {get_vert_lists(ti),get_vert_lists(tj),get_vert_lists(tk),get_vert_lists(tl)};
            tet4D_list_complex.push_back(tet_4d_verts);
        }
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
    std::string time_sample_file = "vert_time_sample.json";
    if (std::filesystem::exists(time_sample_file.c_str())){
        std::filesystem::remove(time_sample_file.c_str());
    }
    using json = nlohmann::json;
    std::ofstream fout(time_sample_file.c_str(),std::ios::app);
    json jOut;
//    jOut.push_back(json(tet4D_list));
    jOut.push_back(json(tet4D_list_complex));
//    jOut.push_back(json(vert_lists));
    fout << jOut.dump(4, ' ', true, json::error_handler_t::replace) << std::endl;
    fout.close();
    std::cout << tet_num << " " << simps.size() << " ";
}

bool save_1_tier_grid(const std::string filename,
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
