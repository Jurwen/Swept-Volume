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

struct pairHash {
    std::size_t operator()(const std::pair<int,int>& p) const noexcept {
        auto a = static_cast<uint32_t>(p.first);
        auto b = static_cast<uint32_t>(p.second);
        uint64_t high = static_cast<uint64_t>(a) << 32;
        uint64_t low  = static_cast<uint64_t>(b);
        uint64_t packed = high | low;
        return static_cast<std::size_t>(packed);
    }
};


// Compute the “validPatch” list as in your Mathematica code.
// - cellNum: number of cells
// - volInfo0: length‑cellNum vector of volumes
// - cellData: list of {u,v} zero‑based adjacency pairs
//
// Returns: a vector of indices i into cellData for which
//          hash[min(u,v)][max(u,v)] was never set to true.
std::pair<std::vector<bool>, std::vector<std::vector<bool>>>
computeValidPatch(size_t cellNum,
                  const std::vector<double>& volInfo0,
                  const arrangement::MatrixIr& cellData)
{
    std::vector<std::vector<bool>> adj(cellNum, std::vector<bool>(cellNum, false));
    std::vector<bool> valid(cellNum, false), merged(cellNum, true);
    std::vector<std::vector<bool>> hash(cellNum, std::vector<bool>(cellNum, false));

    //    We want, for each index i, a “rank” so that largest volInfo0 → rank 0,
    //    next largest → rank 1, etc.
    std::vector<int> order(cellNum);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b){
                  return volInfo0[a] > volInfo0[b];   // descending
              });
    std::vector<int> rankDesc(cellNum);
    for(int pos = 0; pos < cellNum; ++pos) {
        rankDesc[ order[pos] ] = pos;
    }
    const double threshold = 1e-5;
    for(int i = 0; i < cellNum; ++i) {
        if (volInfo0[i] > threshold) valid[i] = true;
    }
    // build adjacency from cellData
    for( size_t e = 0; e < cellData.rows(); e++) {
        int u = cellData(e, 0), v = cellData(e, 1);
        if (u>=0 && u<cellNum && v>=0 && v<cellNum) {
            adj[u][v] = adj[v][u] = true;
        }
    }

    // handle = Xnor/@Transpose[{valid,merged}]
    auto allHandleTrue = [&](){
        for(int i=0;i<cellNum;++i)
            if (valid[i] != merged[i]) return false;
        return true;
    };

    // While[Not@AllTrue@handle, … do merge steps …]
    while(!allHandleTrue()) {
        int target = 0;
        while(target < cellNum && valid[target] == merged[target])
            ++target;
        // find its neighbors: mergeCell = Flatten@Position[adj[[target]],True]
        std::vector<int> neigh;
        neigh.reserve(cellNum);
        for(int j=0; j<cellNum; ++j) {
            if (adj[target][j]) neigh.push_back(j);
        }
        if (neigh.empty()) break;  // nothing to merge with?

        auto it = std::min_element(
            neigh.begin(), neigh.end(),
            [&](int a,int b){ return rankDesc[a] < rankDesc[b]; }
        );
        int mergeCell = *it;

        // union adjacency rows & columns:
        for(int k=0;k<cellNum;++k) {
            adj[mergeCell][k] = adj[mergeCell][k] || adj[target][k];
            adj[k][mergeCell] = adj[k][mergeCell] || adj[k][target];
        }

        hash[target][mergeCell] = true;
        hash[mergeCell][target] = true;

        merged[target] = false;
    }

    // 9) Build validPatch: those edges {u,v} for which hash[min,max]==false
    std::vector<bool> validPatch(cellData.rows(), false);
    for(int i = 0; i < cellData.rows(); ++i) {
        int u = cellData(i, 0),
            v = cellData(i, 1);
        if (! hash[u][v]) {
            validPatch[i] = true;
        }
    }
    return {validPatch, hash};
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
    size_t num_cells = engine->get_num_cells();
    size_t num_patches = engine->get_num_patches();
    size_t num_facets = arrangement_faces.rows();
    std::cout << "Printing arrangement stats: " << std::endl;
    std::cout << "num_cells: " << num_cells << std::endl;
    std::cout << "num_patches: " << num_patches << std::endl;
    std::cout << "num_vertices: " << engine->get_vertices().rows() << std::endl;
    std::cout << "num_faces: " << arrangement_faces.rows() << std::endl;
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
    const double threshold = 1e-5;
    for(int i = 0; i < num_cells; ++i) {
        if (volInfo[i] > threshold) valid[i] = true;
        if (valid[i] && wind_list[i] == 0){
            cellIt.push_back(i);
            //std::cout << i << " ";
        }
    }
    //std::cout << std::endl;
    
    // Pruning of error cells
    const auto& prunedInfo = computeValidPatch(num_cells, volInfo, cell_data);
    const std::vector<std::vector<bool>>& hash = prunedInfo.second;
    const std::vector<bool>& valid_patchInd = prunedInfo.first;

    // Check edge valences for finding the feature lines and corners
    std::vector<std::array<std::array<double, 3>, 2>> feature_lines;
    std::vector<std::array<double, 3>> corners;
    std::unordered_map<std::pair<int, int>, int, pairHash> edge_valence;
    std::unordered_map<int, int> vert_valence;
    auto push_valance = [&](int v0, int v1){
        if (v0 < v1){
            edge_valence[{v0, v1}] ++;
        } else{
            edge_valence[{v1, v0}] ++;
        }
    };
    auto push_vert_valance = [&](int v){
        vert_valence[v] ++;
    };
    for (size_t i = 0; i < arrangement_faces.rows(); i++){
        if (valid_patchInd[patches(i, 0)]){
            const auto& face = arrangement_faces.row(i);
            push_valance(face[0], face[1]);
            push_valance(face[1], face[2]);
            push_valance(face[0], face[2]);
        }
    }
    auto check_valence = [&](int v0, int v1, int v2){
        std::vector<int> vert_list = {v0, v1, v2};
        std::sort(vert_list.begin(), vert_list.end());
        if (edge_valence[{vert_list[0], vert_list[1]}] > 2){
            edge_valence[{vert_list[0], vert_list[1]}] = 2;
            push_vert_valance(vert_list[0]);
            push_vert_valance(vert_list[1]);
            feature_lines.push_back(std::array<std::array<double, 3>, 2>{std::array<double, 3>{out_vertices(vert_list[0], 0), out_vertices(vert_list[0], 1), out_vertices(vert_list[0], 2)}, std::array<double, 3>{out_vertices(vert_list[1], 0), out_vertices(vert_list[1], 1), out_vertices(vert_list[1], 2)}});
        }
        if (edge_valence[{vert_list[0], vert_list[2]}] > 2){
            edge_valence[{vert_list[0], vert_list[2]}] = 2;
            push_vert_valance(vert_list[0]);
            push_vert_valance(vert_list[2]);
            feature_lines.push_back(std::array<std::array<double, 3>, 2>{std::array<double, 3>{out_vertices(vert_list[0], 0), out_vertices(vert_list[0], 1), out_vertices(vert_list[0], 2)}, std::array<double, 3>{out_vertices(vert_list[2], 0), out_vertices(vert_list[2], 1), out_vertices(vert_list[2], 2)}});
        }
        if (edge_valence[{vert_list[1], vert_list[2]}] > 2){
            edge_valence[{vert_list[1], vert_list[2]}] = 2;
            push_vert_valance(vert_list[1]);
            push_vert_valance(vert_list[2]);
            feature_lines.push_back(std::array<std::array<double, 3>, 2>{std::array<double, 3>{out_vertices(vert_list[1], 0), out_vertices(vert_list[1], 1), out_vertices(vert_list[1], 2)}, std::array<double, 3>{out_vertices(vert_list[2], 0), out_vertices(vert_list[2], 1), out_vertices(vert_list[2], 2)}});
        }
    };
    std::vector<std::pair<size_t, bool>> patch_list;
    for (size_t i = 0; i < cellIt.size(); ++i){
        auto cluster = hash[cellIt[i]];
        cluster[cellIt[i]] = true;
        arrangement::MatrixIr patch_faces;
        patch_faces.resize(arrangement_faces.rows(), 3);
        size_t face_count = 0;
        for (size_t f = 0; f < num_facets; ++f){
            auto patch_id = patches(f);
            if (cluster[cell_data(patch_id, 0)] != cluster[cell_data(patch_id, 1)]){
                if (cluster[cell_data(patch_id, 0)]){
                    patch_faces.row(face_count) = (arrangement_faces.row(f));
                }else{
                    patch_faces.row(face_count) = (arrangement_faces.row(f).reverse());
                }
                check_valence(patch_faces(face_count, 0), patch_faces(face_count, 1), patch_faces(face_count, 2));
                face_count++;
            }
        }
        patch_faces.conservativeResize(face_count, 3);
        out_faces.emplace_back(patch_faces);
    }
    for (int i = 0; i < out_vertices.rows(); i++){
        if (vert_valence[i] != 2 && vert_valence[i] == 1){
            corners.push_back({out_vertices(i, 0), out_vertices(i, 1), out_vertices(i, 2)});
        }
    }
    std::string feature_lines_file = output_path + "/features.json";
    if (std::filesystem::exists(feature_lines_file.c_str())){
        std::filesystem::remove(feature_lines_file.c_str());
    }
    using json = nlohmann::json;
    std::ofstream fout(feature_lines_file.c_str(),std::ios::app);
    json jOut;
    jOut.push_back(json(feature_lines));
    jOut.push_back(json(corners));
    fout << jOut.dump(4, ' ', true, json::error_handler_t::replace) << std::endl;
    fout.close();
}
#endif /* post_processing_h */
