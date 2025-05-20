//
//  post_processing.h
//  general_sweep
//
//  Created by Yiwen Ju on 5/13/25.
//

#ifndef post_processing_h
#define post_processing_h
#include <numeric>   // for std::iota
#include <utility>   // for std::pair
#include <limits>    // for std::numeric_limits
#include <arrangement/Arrangement.h>
const double threshold = 1e-5;

struct PairHasher {
    std::size_t operator()(const std::pair<int,int>& p) const noexcept {
        auto a = static_cast<uint32_t>(p.first);
        auto b = static_cast<uint32_t>(p.second);
        if (a > b) std::swap(a, b);
        uint64_t high = static_cast<uint64_t>(a) << 32;
        uint64_t low  = static_cast<uint64_t>(b);
        uint64_t packed = high | low;
        return static_cast<std::size_t>(packed);
    }
};
struct PairEqual {
    bool operator()(std::pair<int,int> x,
                    std::pair<int,int> y) const noexcept
    {
        if (x.first > x.second) std::swap(x.first,x.second);
        if (y.first > y.second) std::swap(y.first,y.second);
        return x.first==y.first && x.second==y.second;
    }
};
constexpr std::pair<int,int> INVALID_PATCH = { -1, -1 };

std::pair<std::vector<bool>, std::vector<int>>
computeValidPatch(int cellNum,
                  std::vector<double>& volInfo,
                  const arrangement::MatrixIr& cellData)
{
    std::vector<std::vector<bool>> adj(cellNum, std::vector<bool>(cellNum, false));
    std::vector<bool> valid(cellNum, false), merged(cellNum, false);
    for(int i = 0; i < cellNum; ++i) {
        if (volInfo[i] > threshold) valid[i] = true;
    }
    // build adjacency from cellData
    for( size_t e = 0; e < cellData.rows(); e++) {
        int u = cellData(e, 0), v = cellData(e, 1);
            adj[u][v] = adj[v][u] = true;
    }
    std::vector<int> owner(cellNum);
    for(int i = 0; i < cellNum; ++i) owner[i] = i;
    int cluster_size = cellNum;
    auto terminate = [&](){
        for(int i = 0; i < cluster_size ; i++)
            if (!(valid[i] || merged[i])) return false;
        return true;
    };

    while(!terminate()) {
        int mergeCell = 0;
        while(mergeCell < cluster_size && (valid[mergeCell] || merged[mergeCell]))
            ++mergeCell;
            // Do not need to consider cells that are either valid or merged.

        int target = -1;
        double target_vol = std::numeric_limits<double>::lowest();
        for (int clusterIt = 0; clusterIt < cluster_size; ++clusterIt) {
            if (!merged[clusterIt] && adj[mergeCell][clusterIt] && volInfo[clusterIt] > target_vol){
                target_vol = volInfo[clusterIt];
                target = clusterIt;
            }
        }
        assert(target < cluster_size);
        if (target >= cluster_size) break;    // nothing left to merge
        // merging two clusters
        cluster_size ++;
        // update ownership
        for(int j = 0; j < cellNum; ++j){
          if(owner[j] == target || owner[j] == mergeCell)
            owner[j] = cluster_size - 1;
        }
        merged[target] = merged[mergeCell] = true;
        merged.push_back(false);
        valid.push_back(valid[target] || valid[mergeCell]);
        volInfo.push_back(std::max(volInfo[target], volInfo[mergeCell]));
        std::vector<bool> new_cluster(cluster_size, false);
        for(int i = 0; i < cluster_size - 1; ++i) {
            new_cluster[i] = adj[target][i] || adj[mergeCell][i];
        }
        new_cluster[target] = new_cluster[mergeCell] = false;
        for (int i = 0; i < cluster_size - 1; ++i){
            adj[i].push_back(new_cluster[i]);
        }
        adj.push_back(new_cluster);
    }
    std::vector<bool> validPatch(cellData.rows(), false);
    for(int i = 0; i < cellData.rows(); ++i) {
        int u = cellData(i, 0),
            v = cellData(i, 1);
        if (owner[u] != owner[v]) {
            validPatch[i] = true;
        }
    }
    return {validPatch, owner};
}


/**
 * Compute the surface of the sweep volume.
 *
 * @param vertices Input vertices of the sweep generator.
 * @param faces Input faces of the sweep generator.
 * @param out_vertices Output vertices of the sweep volume surface.
 * @param out_faces Output faces of the sweep volume surface.
 */
void compute_sweep_volume(const arrangement::MatrixFr& vertices, const arrangement::MatrixIr& faces,
                          arrangement::MatrixFr& out_vertices, std::vector<arrangement::MatrixIr>& out_faces, std::string output_path) {
    // assume the input faces have already been oriented based on second order time derivatives
    // then the sweep volume consists of arrangement cells with positive winding number
    
    // Initialize face labels
    arrangement::VectorI face_labels = Eigen::VectorXi::LinSpaced(faces.rows(), 0, faces.rows() - 1);
    // create a mesh arrangement engine.
    auto engine = arrangement::Arrangement::create_mesh_arrangement(vertices, faces, face_labels);
    // Alternatively, create a fast mesh arrangement engine.
    // face_labels.setZero();  // fast arrangement cannot use too large face labels
    // auto engine = arrangement::Arrangement::create_fast_arrangement(vertices, faces, face_labels);
    // or, create a Geogram arrangement engine
    // the cell ids associated with patches can often be MAX_INT64.
    // auto engine = arrangement::Arrangement::create_geogram_arrangement(vertices, faces, face_labels);
    
// The following is based on James' arrangement code in Python: https://github.com/qnzhou/arrangement-benchmark/blob/main/python/arrangement/__main__.py
    engine->run();
    
    const auto& cell_data       = engine->get_cells();           // (#facets x 2) array
    const auto& patches         = engine->get_patches();         // list of facet indices
    const auto& parent_facets   = engine->get_out_face_labels();   // size = #facets
    const auto& winding_number  = engine->get_winding_number();  // (#facets x 2) array
    out_vertices = engine->get_vertices();
    const auto& arrangement_faces = engine->get_faces();
    int num_cells = engine->get_num_cells();
    int num_patches = engine->get_num_patches();
    int num_facets = arrangement_faces.rows();
    std::vector<int> wind_list(num_cells);
    std::vector<double> volInfo(num_cells);
    std::vector<size_t> cellIt;
    for (size_t i = 0; i < num_cells; ++i) {
        // collect all facets incident on cell i
        std::vector<int> active_facets;
        active_facets.reserve(num_patches);
        int active_facets_count = 0;
        bool active_orientations;
        double meshVol = 0.0;
        for (int facet_idx = 0; facet_idx < patches.size(); ++facet_idx) {
            int c0 = cell_data(patches[facet_idx], 0);
            int c1 = cell_data(patches[facet_idx], 1);
            if (c0 == i || c1 == i) {
                active_facets.emplace_back(facet_idx);
                active_facets_count++;
                int j = arrangement_faces(facet_idx, 0), k = arrangement_faces(facet_idx, 1), l = arrangement_faces(facet_idx, 2);
                Eigen::Vector3d vj = out_vertices.row(j);
                Eigen::Vector3d vk = out_vertices.row(k);
                Eigen::Vector3d vl = out_vertices.row(l);
                meshVol += vj.dot( vk.cross(vl) );
                if (active_facets_count == 1){
                    active_orientations = (c0 == i);
                }
            }
        }
        meshVol = std::abs(meshVol) / 6.0;
        volInfo[i] = meshVol;
        active_facets.resize(active_facets_count);
        if (active_facets.empty()) {
            continue;
        }
        
        // pick winding from the first active facet, based on its orientation
        int active_wind;
        if (!active_orientations) {
            active_wind = winding_number(active_facets[0], 1);
        } else {
            active_wind = winding_number(active_facets[0], 0);
        }
        wind_list[i] = active_wind;
    }
    std::vector<bool> valid(num_cells, false);
    int valid_num = 0;
    std::cout << "valid 0-winding cell iter: ";
    for(int i = 0; i < num_cells; ++i) {
        if (volInfo[i] > threshold) {
            valid[i] = true;
            valid_num++;
        }
        if (valid[i] && wind_list[i] == 0){
            cellIt.push_back(i);
            std::cout << i << " ";
        }
    }
    std::cout << std::endl;
    std::cout << "number of valid cells: " << valid_num << std::endl;
    
    // Pruning of error cells
    const auto& prunedInfo2 = computeValidPatch(num_cells, volInfo, cell_data);
    std::vector<bool> valid_patchInd = prunedInfo2.first;
    std::vector<int> cell_merge_map = prunedInfo2.second;
//    std::vector<bool> valid_patchInd (cell_data.size(), true);
//    std::vector<int> cell_merge_map (num_cells, -1);
//    for (int i = 0; i < num_cells; i++){
//        cell_merge_map[i] = i;
//    }
    
    std::vector<std::array<std::array<double, 3>, 2>> feature_lines;
    std::vector<std::array<double, 3>> corners;
    std::unordered_map<int, std::pair<int, int>> patchIdMap;
    for (int i = 0; i < valid_patchInd.size(); i++){
        if (valid_patchInd[i]){
            int c0 = cell_data(i, 0);
            int c1 = cell_data(i, 1);
            patchIdMap[i].first = cell_merge_map[c0];
            patchIdMap[i].second = cell_merge_map[c1];
        }
    }
    
    // New detection method to check if multiple patches meet (valences information is also needed after testing)
    std::unordered_map<std::pair<int, int>, std::pair<int, int>, PairHasher, PairEqual> edge_patch_map;
    std::unordered_map<std::pair<int, int>, int, PairHasher, PairEqual> edge_valence_map;
    std::unordered_map<int, int> vert_valence;
    auto push_vert_valance = [&](int v){
        vert_valence[v] ++;
    };
    auto push_edge_patch = [&](std::pair<int, int> edge, int patchId){
        int c0 = cell_data(patchId, 0), c1 = cell_data(patchId, 1);
        if (c0 > c1){
            std::swap(c0, c1);
        }
        edge_patch_map[edge] = std::pair<int, int>{c0, c1};
        edge_valence_map[edge]++;
    };
    auto check_edge_patch = [&](std::pair<int, int> edge, int patchId){
        edge_valence_map[edge]++;
        auto prev_patch = edge_patch_map[edge];
        if (prev_patch != INVALID_PATCH){
            PairHasher hasher;
            if (hasher(prev_patch) != hasher({cell_data(patchId, 0), cell_data(patchId, 1)})){
                edge_patch_map[edge] = INVALID_PATCH;
            }
        }
    };
    auto parse_edge_patch = [&] (std::pair<int, int> edge){
        if (edge_patch_map[edge] == INVALID_PATCH && edge_valence_map[edge] > 2){
            edge_valence_map[edge] = 2;
            push_vert_valance(edge.first);
            push_vert_valance(edge.second);
            feature_lines.push_back(std::array<std::array<double, 3>, 2>{std::array<double, 3>{out_vertices(edge.first, 0), out_vertices(edge.first, 1), out_vertices(edge.first, 2)}, std::array<double, 3>{out_vertices(edge.second, 0), out_vertices(edge.second, 1), out_vertices(edge.second, 2)}});
        }
    };

    for (size_t i = 0; i < arrangement_faces.rows(); i++){
        auto& p = patches(i);
        if (valid_patchInd[p]){
            const auto& face = arrangement_faces.row(i);
            std::array<std::pair<int, int>, 3> edge =
                {std::pair<int, int>{face[0], face[1]},
                std::pair<int, int>{face[0], face[2]},
                std::pair<int, int>{face[1], face[2]}};
            for (auto& e: edge){
                if (edge_patch_map.contains(e)){
                    check_edge_patch(e, p);
                }else{
                    push_edge_patch(e, p);
                }
            }
        }
    }
    
    std::vector<std::pair<size_t, bool>> patch_list;
    for (size_t i = 0; i < cellIt.size(); ++i){
        auto current_cell = cell_merge_map[cellIt[i]];
//        auto current_cell = cellIt[i];
        arrangement::MatrixIr patch_faces;
        patch_faces.resize(arrangement_faces.rows(), 3);
        size_t face_count = 0;
        for (size_t f = 0; f < num_facets; ++f){
            auto patch_id = patches(f);
            auto c0 = cell_merge_map[cell_data(patch_id, 0)], c1 = cell_merge_map[cell_data(patch_id, 1)];
//            auto c0 = cell_data(patch_id, 0), c1 = cell_data(patch_id, 1);
            if ((c0 == current_cell && c1 != current_cell) || (c1 == current_cell && c0 != current_cell)){
                if (c0 == current_cell){
                    patch_faces.row(face_count) = (arrangement_faces.row(f));
                }else{
                    patch_faces.row(face_count) = (arrangement_faces.row(f).reverse());
                }
                const auto& face = arrangement_faces.row(f);
                std::array<std::pair<int, int>, 3> edge =
                    {std::pair<int, int>{face[0], face[1]},
                    std::pair<int, int>{face[0], face[2]},
                    std::pair<int, int>{face[1], face[2]}};
                for (auto& e: edge){
                    parse_edge_patch(e);
                }
                face_count++;
            }
        }
        patch_faces.conservativeResize(face_count, 3);
        out_faces.emplace_back(patch_faces);
    }
    int count = 0;
    for (int i = 0; i < out_vertices.rows(); i++){
        if (vert_valence[i] != 2 && vert_valence[i] > 0){
            corners.push_back({out_vertices(i, 0), out_vertices(i, 1), out_vertices(i, 2)});
        }
    }
    std::cout << "# Verts: " << engine->get_vertices().rows() << std::endl;
    std::cout << "# Tris: " << arrangement_faces.rows() << std::endl;
    std::cout << "patches: " << num_patches << std::endl;
    std::cout << "cells: " << num_cells << std::endl;
    std::string feature_lines_file = output_path + "/features.json";
    if (std::filesystem::exists(feature_lines_file.c_str())){
        std::filesystem::remove(feature_lines_file.c_str());
    }
    using json = nlohmann::json;
    std::ofstream fout(feature_lines_file.c_str(),std::ios::app);
    json jOut;
    jOut["chains"] = feature_lines;
    jOut["joints"] = corners;
//    jOut.push_back(json(feature_lines));
//    jOut.push_back(json(corners));
    fout << jOut.dump(4, ' ', true, json::error_handler_t::replace) << std::endl;
    fout.close();
}
#endif /* post_processing_h */
