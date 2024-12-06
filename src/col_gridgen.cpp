//
//  col_gridgen.cpp
//  adaptive_column_grid
//
//  Created by Yiwen Ju on 12/4/24.
//

#include "col_gridgen.h"
///   Lexicographical comparison for the extruding the 4D simplex from the base tetrahedra This first comparies the lowest time stamp, the second lowest time stamp, and the vertex id.
auto comp = [](const std::tuple<int, int, uint64_t>& e0, const std::tuple<int, int, uint64_t>& e1) {
    return e0 < e1; // Lexicographical comparison
};

/// Sample the list of 5-cells based on the base tetrahedra and 4 lists of time samples at its vertices. The extrusion/sampling is based on lowest time stamp, the second lowest time stamp, and the vertex id comparing the four incremental time stamps at each vertex.
/// @param[in] grid: the base tetrahedra grid in `mtet` structure
/// @param[in] tid: the tetrahedra id that is going to be extruded
/// @param[in] timeMap: the map from vertex to a list of time samples.
/// @return A list of cell5 elements.
std::vector<cell5> sampleCol(std::span<mtet::VertexId, 4> vs, vertexCol &timeMap){
    std::vector<vertex4d> ti = timeMap[value_of(vs[0])];
    std::vector<vertex4d> tj = timeMap[value_of(vs[1])];
    std::vector<vertex4d> tk = timeMap[value_of(vs[2])];
    std::vector<vertex4d> tl = timeMap[value_of(vs[3])];
    std::array<uint64_t, 4> quad = {value_of(vs[0]), value_of(vs[1]), value_of(vs[2]), value_of(vs[3])};
    vertex4d last;
    last.time = ti.back().time * 2;
    ti.push_back(last);
    tj.push_back(last);
    tk.push_back(last);
    tl.push_back(last);
    std::vector<cell5> cell5Col;
    cell5Col.reserve(ti.size()+tj.size()+tk.size()+tl.size() - 8);
    int i = 1, j = 1, k = 1, l = 1;
    while (i < ti.size() - 1 || j < tj.size() - 1 || k < tk.size() - 1 || l < tl.size() - 1) {
        std::array<std::tuple<int, int, uint64_t>, 4> candidates = {
            std::tuple<int, int, uint64_t>{ti[i - 1].time, ti[i].time, quad[0]},
            {tj[j - 1].time, tj[j].time, quad[1]},
            {tk[k - 1].time, tk[k].time, quad[2]},
            {tl[l - 1].time, tl[l].time, quad[3]}};
        // Find the index of the minimum tuple
        //auto minIt = std::min_element(candidates.begin(), candidates.end(), comp);
        size_t minInd = std::distance(candidates.begin(), std::min_element(candidates.begin(), candidates.end(), comp));
        cell5 simp;
        simp.level = 0; 
        switch (minInd) {
            case 0:
                ++i;
                simp.hash = {i - 1, j - 1, k - 1, l - 1, 0};
                break;
            case 1:
                ++j;
                simp.hash = {i - 1, j - 1, k - 1, l - 1, 1};
                break;
            case 2:
                ++k;
                simp.hash = {i - 1, j - 1, k - 1, l - 1, 2};
                break;
            case 3:
                ++l;
                simp.hash = {i - 1, j - 1, k - 1, l - 1, 3};
                break;
        }
        cell5Col.push_back(simp);
//        std::cout << simp.hash[0] << " " << simp.hash[1] << " " << simp.hash[2] << " " << simp.hash[3] << " " << simp.hash[4] << " " << std::endl;
    }
    return cell5Col;
}

/// @param[in] timeDep: depth of time intervals; exponential increase
/// @param[in] grid: base 3D grid. For each vertex, build a list of time stamps. For each tet, build a list of extruded 4D simplices
/// @param[in] func: the implicit function that represents the swept volume. The input of the function is the 4d coordinate, and the output is an size-4 vector with first entry as the value and the other three as the gradient. 
/// @param[in] maxTimeDep: maximum interger-valued time depth of the trajectory. Default: 1024
///
/// @param[out] timeList: a list of time stamps at this vertex
void init5CGrid(const int timeDep, mtet::MTetMesh grid, const std::function<Eigen::RowVector4d(std::span<const Scalar, 4>)> func, const int maxTimeDep, vertexCol &timeMap, tetCol &cell5Map){
    int timeLen = pow(2, timeDep - 1);
    int len = maxTimeDep / timeLen;
    std::vector<int> time3DList(timeLen + 1);
    for (int i = 0; i < timeLen+1; i++){
        int time = i * len;
        time3DList[i] = time;
    }
    grid.seq_foreach_vertex([&](mtet::VertexId vid, std::span<const mtet::Scalar, 3> data)
                            {
        std::vector<vertex4d> vertColList(timeLen + 1);
        for (int i = 0; i < timeLen + 1; i++){
            vertex4d vert;
            vert.time = time3DList[i];
            double time_fp = (double)vert.time / 1024;
            vert.coord = {data[0], data[1], data[2], time_fp};
//            std::cout << vert.coord[0] << " " << vert.coord[1] << " " << vert.coord[2] << " " << vert.coord[3] << std::endl;
            vert.valGradList = func(vert.coord);
            std::cout << vert.valGradList[0] << " " << vert.valGradList[1] << " " << vert.valGradList[2] << " " << vert.valGradList[3] << std::endl;
            vertColList[i] = vert;
        }
        timeMap[value_of(vid)] = vertColList;
        
    });
    grid.seq_foreach_tet([&](mtet::TetId tid, [[maybe_unused]] std::span<const mtet::VertexId, 4> data){
        std::span<VertexId, 4> vs = grid.get_tet(tid);
        std::vector<cell5> col = sampleCol(vs, timeMap);
        
        cell5Map[vs] = col;
    });
}

bool gridRefine(mtet::MTetMesh grid, vertexCol &timeMap, tetCol &cell5Map, const std::function<Eigen::RowVector4d(std::span<const Scalar, 4>)> func, const double threshold){
    init5CGrid(3, grid, func, 1024, timeMap, cell5Map);
    
    auto compTime = [](std::tuple<mtet::Scalar, cell5, mtet::VertexId, size_t> timeSub0,
                       std::tuple<mtet::Scalar, cell5, mtet::VertexId, size_t> timeSub1)
    { return std::get<0>(timeSub0) < std::get<0>(timeSub1); };
    std::vector<std::tuple<mtet::Scalar, cell5, mtet::VertexId, size_t>> timeQ;
    
    auto compSpace = [](std::tuple<mtet::Scalar, cell5, mtet::VertexId, size_t> timeSub0,
                       std::tuple<mtet::Scalar, cell5, mtet::VertexId, size_t> timeSub1)
    { return std::get<0>(timeSub0) < std::get<0>(timeSub1); };
    std::vector<std::tuple<mtet::Scalar, cell5, mtet::VertexId, size_t>> spaceQ;
    
    
    
    return true;
}
