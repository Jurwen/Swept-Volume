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
/// @param[in] vertexMap: the map from vertex to a list of time samples.
/// @return A list of cell5 elements.
simpCol::cell5_list sampleCol(std::span<mtet::VertexId, 4> vs, vertExtrude &vertexMap){
    vertexCol::vert4d_list ti = vertexMap[value_of(vs[0])].timeStamp;
    vertexCol::vert4d_list tj = vertexMap[value_of(vs[1])].timeStamp;
    vertexCol::vert4d_list tk = vertexMap[value_of(vs[2])].timeStamp;
    vertexCol::vert4d_list tl = vertexMap[value_of(vs[3])].timeStamp;
    std::array<uint64_t, 4> quad = {value_of(vs[0]), value_of(vs[1]), value_of(vs[2]), value_of(vs[3])};
    vertex4d last;
    last.time = ti.back().time * 2;
    ti.emplace_back(last);
    tj.emplace_back(last);
    tk.emplace_back(last);
    tl.emplace_back(last);
    simpCol::cell5_list cell5Col;
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
        cell5Col.emplace_back(simp);
    }
    return cell5Col;
}

/// @param[in] timeDep: depth of time intervals; exponential increase
/// @param[in] grid: base 3D grid. For each vertex, build a list of time stamps. For each tet, build a list of extruded 4D simplices
/// @param[in] func: the implicit function that represents the swept volume. The input of the function is the 4d coordinate, and the output is an size-4 vector with first entry as the value and the other three as the gradient. 
/// @param[in] maxTimeDep: maximum interger-valued time depth of the trajectory. Default: 1024
///
/// @param[out] timeList: a list of time stamps at this vertex
void init5CGrid(const int timeDep, mtet::MTetMesh grid, const std::function<std::pair<Scalar, Eigen::RowVector4d>(Eigen::RowVector4d)> func, const int maxTimeDep, vertExtrude &vertexMap, tetExtrude &cell5Map){
    int timeLen = pow(2, timeDep - 1);
    int len = maxTimeDep / timeLen;
    vertexCol::time_list time3DList(timeLen + 1);
    for (int i = 0; i < timeLen+1; i++){
        int time = i * len;
        time3DList[i] = time;
    }
    grid.seq_foreach_vertex([&](mtet::VertexId vid, std::span<const mtet::Scalar, 3> data)
                            {
        vertexCol col;
        vertexCol::vert4d_list vertColList(timeLen + 1);
        for (int i = 0; i < timeLen + 1; i++){
            vertex4d vert;
            vert.time = time3DList[i];
            double time_fp = (double)vert.time / 1024;
            vert.coord = {data[0], data[1], data[2], time_fp};
            vert.valGradList = func(vert.coord);
            vertColList[i] = vert;
        }
        col.timeStamp = vertColList;
        vertexMap[value_of(vid)] = col;
    });
    grid.seq_foreach_tet([&](mtet::TetId tid, [[maybe_unused]] std::span<const mtet::VertexId, 4> data){
        std::span<VertexId, 4> vs = grid.get_tet(tid);
        for (size_t i = 0; i < 4; i++){
            vertexMap[value_of(vs[i])].vertTetAssoc.push_back(tid);
        }
        simpCol col;
        col.cell5Col = sampleCol(vs, vertexMap);
        cell5Map[vs] = col;
    });
}

bool gridRefine(mtet::MTetMesh &grid, vertExtrude &vertexMap, tetExtrude &cell5Map, const std::function<std::pair<Scalar, Eigen::RowVector4d>(Eigen::RowVector4d)> func, const double threshold, const int max_splits){
    init5CGrid(3, grid, func, 1024, vertexMap, cell5Map);
    ///
    /// Initiate queue: timeQ and spaceQ
    auto compTime = [](std::tuple<mtet::Scalar, mtet::TetId, mtet::VertexId, int, size_t, int> timeSub0,
                       std::tuple<mtet::Scalar, mtet::TetId, mtet::VertexId, int, size_t, int> timeSub1)
    { return std::get<0>(timeSub0) <= std::get<0>(timeSub1); };
    std::vector<std::tuple<mtet::Scalar, mtet::TetId, mtet::VertexId, int, size_t, int>> timeQ;
    
    auto compSpace = [](std::tuple<mtet::Scalar, mtet::TetId, mtet::EdgeId, int> spaceSub0,
                       std::tuple<mtet::Scalar, mtet::TetId, mtet::EdgeId, int> spaceSub1)
    { return std::get<0>(spaceSub0) <= std::get<0>(spaceSub1); };
    std::vector<std::tuple<mtet::Scalar, mtet::TetId, mtet::EdgeId, int>> spaceQ;
    
    
    auto queryCell5 = [&](cell5 simp, std::array<vertexCol, 4> vs){
        
    };
    int splits = 0;
    ///
    /// Push Queue Function:
    auto push_one_col = [&](mtet::TetId tid)
    {
        std::span<VertexId, 4> vs = grid.get_tet(tid);
        simpCol colInfo = cell5Map[vs];
        simpCol::cell5_list cell5Col = colInfo.cell5Col;
        std::array<vertexCol, 4> baseVerts;
        for (size_t i = 0; i < 4; i++){
            baseVerts[i] = vertexMap[value_of(vs[i])];
        }
        /// Compute longest spatial edge
        mtet::EdgeId longest_edge;
        mtet::Scalar longest_edge_length = 0;
        bool baseSub = false;
        grid.foreach_edge_in_tet(tid, [&](mtet::EdgeId eid, mtet::VertexId v0, mtet::VertexId v1)
                                 {
            auto p0 = grid.get_vertex(v0);
            auto p1 = grid.get_vertex(v1);
            mtet::Scalar l = (p0[0] - p1[0]) * (p0[0] - p1[0]) + (p0[1] - p1[1]) * (p0[1] - p1[1]) +
            (p0[2] - p1[2]) * (p0[2] - p1[2]);
            if (l > longest_edge_length) {
                longest_edge_length = l;
                longest_edge = eid;
            } });
        for (size_t cell5It = 0; cell5It < cell5Col.size(); cell5It++){
            std::array<int, 5> cell5 = cell5Col[cell5It].hash;
            int lastInd = cell5[4];
            std::array<vertex4d, 5> verts;
            verts[0] = baseVerts[lastInd].timeStamp[cell5[lastInd]];
            size_t ind = 0;
            for (size_t i = 0; i < 4; i++){
                if (i != lastInd){
                    ind ++;
                    verts[ind] = baseVerts[i].timeStamp[cell5[i]];
                }
            }
            verts[4] = baseVerts[lastInd].timeStamp[cell5[lastInd] - 1];
            bool inside = false;
            bool ret = refine4D(verts, threshold, inside);
            if (ret){
                double timeLens =verts[0].coord[3] - verts[4].coord[3];
                if (timeLens * timeLens > longest_edge_length && (verts[0].time - verts[4].time) > 4){
                    timeQ.emplace_back(timeLens, tid, vs[lastInd], (verts[0].time + verts[4].time) / 2,(size_t)cell5[lastInd], colInfo.level);
                    std::push_heap(timeQ.begin(), timeQ.end(), compTime);
                }else{
                    if (!baseSub){
                        spaceQ.emplace_back(longest_edge_length, tid, longest_edge, colInfo.level);
                        std::push_heap(spaceQ.begin(), spaceQ.end(), compSpace);
                        baseSub = true;
                    }
                }
            }
        }
    };
    
    //alternative speed-up approach
//    auto push_one_col = [&](mtet::TetId tid)
//    {
//        std::span<VertexId, 4> vs = grid.get_tet(tid);
//        simpCol colInfo = cell5Map[vs];
//        simpCol::cell5_list cell5Col = colInfo.cell5Col;
//        std::array<vertexCol, 4> baseVerts;
//        for (size_t i = 0; i < 4; i++){
//            baseVerts[i] = vertexMap[value_of(vs[i])];
//        }
//        /// Compute longest spatial edge
//        mtet::EdgeId longest_edge;
//        mtet::Scalar longest_edge_length = 0;
////        mtet::VertexId subVertex0, subVertex1;
//        bool baseSub = false;
//        double longest_time = 0;
//        size_t longest_time_ind = -1;
//        grid.foreach_edge_in_tet(tid, [&](mtet::EdgeId eid, mtet::VertexId v0, mtet::VertexId v1)
//                                 {
//            auto p0 = grid.get_vertex(v0);
//            auto p1 = grid.get_vertex(v1);
//            mtet::Scalar l = (p0[0] - p1[0]) * (p0[0] - p1[0]) + (p0[1] - p1[1]) * (p0[1] - p1[1]) +
//            (p0[2] - p1[2]) * (p0[2] - p1[2]);
//            if (l > longest_edge_length) {
//                longest_edge_length = l;
//                longest_edge = eid;
////                subVertex0 = v0;
////                subVertex1 = v1;
//            } });
//        for (size_t cell5It = 0; cell5It < cell5Col.size(); cell5It++){
//            std::array<int, 5> cell5 = cell5Col[cell5It].hash;
//            int lastInd = cell5[4];
//            std::array<vertex4d, 5> verts;
//            verts[0] = baseVerts[lastInd].timeStamp[cell5[lastInd]];
//            size_t ind = 0;
//            for (size_t i = 0; i < 4; i++){
//                if (i != lastInd){
//                    ind ++;
//                    verts[ind] = baseVerts[i].timeStamp[cell5[i]];
//                }
//            }
//            verts[4] = baseVerts[lastInd].timeStamp[cell5[lastInd] - 1];
//            bool inside = false;
//            bool ret = refine4D(verts, 0.05, inside);
//            if (ret){
//                //std::cout << splits << std::endl;
//                double timeLens =verts[0].coord[3] - verts[4].coord[3];
//                if (timeLens > longest_edge_length && (verts[0].time - verts[4].time) > 4){
//                    if (timeLens > longest_time){
//                        longest_time = timeLens;
//                        longest_time_ind = cell5It;
//                    }
////                    timeQ.emplace_back(timeLens, tid, vs[lastInd], (verts[0].time + verts[4].time) / 2,(size_t)cell5[lastInd], colInfo.level);
////                    std::push_heap(timeQ.begin(), timeQ.end(), compTime);
//                }else{
//                    baseSub = true;
//                }
////                else{
////                    if (!baseSub){
////                        spaceQ.emplace_back(longest_edge_length, tid, longest_edge, colInfo.level);
////                        std::push_heap(spaceQ.begin(), spaceQ.end(), compSpace);
////                        baseSub = true;
////                    }
//                //                }
//            }
//        }
//        if (longest_time_ind != -1){
//            std::array<int, 5> cell5 = cell5Col[longest_time_ind].hash;
//            int lastInd = cell5[4];
//            std::array<vertex4d, 5> verts;
//            verts[0] = baseVerts[lastInd].timeStamp[cell5[lastInd]];
//            size_t ind = 0;
//            for (size_t i = 0; i < 4; i++){
//                if (i != lastInd){
//                    ind ++;
//                    verts[ind] = baseVerts[i].timeStamp[cell5[i]];
//                }
//            }
//            verts[4] = baseVerts[lastInd].timeStamp[cell5[lastInd] - 1];
//            timeQ.emplace_back(longest_edge_length, tid, vs[lastInd], (verts[0].time + verts[4].time) / 2,(size_t)cell5[lastInd], colInfo.level);
//            std::push_heap(timeQ.begin(), timeQ.end(), compTime);
//        }else if (baseSub){
//            spaceQ.emplace_back(longest_edge_length, tid, longest_edge, colInfo.level);
//            std::push_heap(spaceQ.begin(), spaceQ.end(), compSpace);
//        }
//    };
    
    grid.seq_foreach_tet([&](mtet::TetId tid, [[maybe_unused]] std::span<const mtet::VertexId, 4> vs)
                         { push_one_col(tid); });
    while ((!timeQ.empty() || !spaceQ.empty()) && splits < max_splits){
        if (!timeQ.empty()){
            /// temporal subdivision:
            //std::cout << "temporal sub" << std::endl;
            std::pop_heap(timeQ.begin(), timeQ.end(), compTime);
            auto [score, tid, vid, time, cind, level] = timeQ.back();
            timeQ.pop_back();
            
            if (!grid.has_tet(tid)){
                continue;
            }
            std::span<VertexId, 4> vs = grid.get_tet(tid);
            vertexCol timeList = vertexMap[value_of(vid)];
            simpCol cell5List = cell5Map[vs];
            if (level != cell5List.level || timeList.timeExist[time]){
                continue;
            }
            splits ++;
            timeList.timeExist[time] = true;
            vertexCol::vertToSimp newTets;
            newTets.reserve(timeList.vertTetAssoc.size());
            for (size_t i = 0; i < timeList.vertTetAssoc.size(); i++){
                auto assocTet = timeList.vertTetAssoc[i];
                if (grid.has_tet(assocTet)){
                    newTets.emplace_back(assocTet);
                }
            }
            timeList.vertTetAssoc = newTets;
            vertex4d newVert;
            newVert.time = time;
            newVert.coord = {timeList.timeStamp[0].coord[0], timeList.timeStamp[0].coord[1], timeList.timeStamp[0].coord[2], (double)time / 1024};
            newVert.valGradList = func(newVert.coord);
            timeList.timeStamp.insert(timeList.timeStamp.begin() + cind, newVert);
            vertexMap[value_of(vid)] = timeList;
            for (size_t i = 0; i < newTets.size(); i++){
                std::span<VertexId, 4> newVs = grid.get_tet(newTets[i]);
                cell5Map[newVs].cell5Col = sampleCol(newVs, vertexMap);
                cell5Map[newVs].level ++;
                push_one_col(newTets[i]);
            }
        }else{
            /// spatial subdivision:
            //std::cout << "spatial sub" << std::endl;
            std::pop_heap(spaceQ.begin(), spaceQ.end(), compSpace);
            auto [score, tid, eid, level] = spaceQ.back();
            spaceQ.pop_back();
            if (!grid.has_tet(tid)){
                continue;
            }
            std::span<VertexId, 4> vs = grid.get_tet(tid);
            simpCol cell5List = cell5Map[vs];
            if (level != cell5List.level){
                continue;
            }
            splits ++;
            vertexCol newVert;
            std::array<VertexId, 2> vs_old = grid.get_edge_vertices(eid);
            auto [vid, eid0, eid1] = grid.split_edge(eid);
            std::span<Scalar, 3> vidCoord = grid.get_vertex(vid);
            vertexCol v0Col = vertexMap[value_of(vs_old[0])];
            vertexCol v1Col = vertexMap[value_of(vs_old[1])];
            vertexCol::time_list v0Time = v0Col.getTimeList();
            vertexCol::time_list v1Time = v1Col.getTimeList();
            vertexCol::time_list tSamples;
            std::set_union(v0Time.begin(), v0Time.end(), v1Time.begin(), v1Time.end(), std::back_inserter(tSamples));
            vertexCol newVertCol;
            vertexCol::vert4d_list vertColList(tSamples.size());
            for (int i = 0; i < tSamples.size(); i++){
                vertex4d vert;
                vert.time = tSamples[i];
                double time_fp = (double)vert.time / 1024;
                vert.coord = {vidCoord[0], vidCoord[1], vidCoord[2], time_fp};
                vert.valGradList = func(vert.coord);
                vertColList[i] = vert;
            }
            newVertCol.timeStamp = vertColList;
            vertexMap[value_of(vid)] = newVertCol;
            grid.foreach_tet_around_edge(eid0, [&](mtet::TetId t0)
                                         {
                std::span<VertexId, 4> vs = grid.get_tet(t0);
                for (size_t i = 0; i < 4; i++){
                    vertexMap[value_of(vs[i])].vertTetAssoc.push_back(t0);
                }
                simpCol col;
                col.cell5Col = sampleCol(vs, vertexMap);
                cell5Map[vs] = col;
                push_one_col(t0);
            });
            grid.foreach_tet_around_edge(eid1, [&](mtet::TetId t1)
                                         {
                std::span<VertexId, 4> vs = grid.get_tet(t1);
                for (size_t i = 0; i < 4; i++){
                    vertexMap[value_of(vs[i])].vertTetAssoc.push_back(t1);
                }
                simpCol col;
                col.cell5Col = sampleCol(vs, vertexMap);
                cell5Map[vs] = col;
                push_one_col(t1);
            });
        }
    }
    
    return true;
}
