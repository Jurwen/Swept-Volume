//
//  col_gridgen.cpp
//  adaptive_column_grid
//
//  Created by Yiwen Ju on 12/4/24.
//

#include "col_gridgen.h"

namespace dr = drjit; // For nanothread

/// Sample the list of 5-cells based on the base tetrahedra and 4 lists of time samples at its vertices. The extrusion/sampling is based on lowest time stamp, the second lowest time stamp, and the vertex id comparing the four incremental time stamps at each vertex.
/// @param[in] grid: the base tetrahedra grid in `mtet` structure
/// @param[in] tid: the tetrahedra id that is going to be extruded
/// @param[in] vertexMap: the map from vertex to a list of time samples.
/// @return A list of cell5 elements.
simpCol::cell5_list sampleCol(std::span<mtet::VertexId, 4> vs, vertExtrude &vertexMap){
    vertexCol::vert4d_list ti = vertexMap[value_of(vs[0])].vert4dList;
    vertexCol::vert4d_list tj = vertexMap[value_of(vs[1])].vert4dList;
    vertexCol::vert4d_list tk = vertexMap[value_of(vs[2])].vert4dList;
    vertexCol::vert4d_list tl = vertexMap[value_of(vs[3])].vert4dList;
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
        std::array<std::pair<int, uint64_t>, 4> candidates = {
            std::pair<int, uint64_t>{ti[i].time, quad[0]},
            {tj[j].time, quad[1]},
            {tk[k].time, quad[2]},
            {tl[l].time, quad[3]}};
        // Find the index of the minimum tuple
        size_t minInd = std::distance(candidates.begin(), std::min_element(candidates.begin(), candidates.end()));
        cell5 simp;
        switch (minInd) {
            case 0:
                ++i;
                simp.hash = {i - 1, j - 1, k - 1, l - 1, 0};
                simp.time_list = {ti[i - 1].time, tj[j - 1].time, tk[k - 1].time, tl[l - 1].time, ti[i-2].time};
                break;
            case 1:
                ++j;
                simp.hash = {i - 1, j - 1, k - 1, l - 1, 1};
                simp.time_list = {ti[i - 1].time, tj[j - 1].time, tk[k - 1].time, tl[l - 1].time, tj[j-2].time};
                break;
            case 2:
                ++k;
                simp.hash = {i - 1, j - 1, k - 1, l - 1, 2};
                simp.time_list = {ti[i - 1].time, tj[j - 1].time, tk[k - 1].time, tl[l - 1].time, tk[k-2].time};
                break;
            case 3:
                ++l;
                simp.hash = {i - 1, j - 1, k - 1, l - 1, 3};
                simp.time_list = {ti[i - 1].time, tj[j - 1].time, tk[k - 1].time, tl[l - 1].time, tl[l-2].time};
                break;
        }
        cell5Col.emplace_back(simp);
    }
    return cell5Col;
}

/// Re-extrusion of simplex column after the temporal edge subdivision
/// @param[in] simpInfo:  The list of 4D simplex column
/// @param[in] time: The integer-valued time of the newly-inserted 4D vertex after the temporal edge split
/// @param[in] quad: The indices of four 3D vertices used by this simplex column as the base
llvm_vecsmall::SmallVector<size_t, 256> resampleTimeCol
(simpCol::cell5_list simpInfo, const int time, const size_t ind){
    llvm_vecsmall::SmallVector<size_t, 256> refineList;
    refineList.reserve(simpInfo.size());
    size_t cell5_num = 0;
    for (size_t i = 0; i < simpInfo.size(); i++){
        //std::cout << time << " " << simpInfo[i].bot(ind) << " " << simpInfo[i].top() << " " << simpInfo[i].hash[4] << " " << ind << std::endl;
        if (simpInfo[i].bot(ind) == time || (simpInfo[i].top() == time && simpInfo[i].hash[4] == ind)){
            cell5_num++;
            refineList.emplace_back(i);
        }
    }
    refineList.resize(cell5_num);
    return refineList;
}

/// @param[in] timeDep: depth of time intervals; exponential increase
/// @param[in] grid: base 3D grid. For each vertex, build a list of time stamps. For each tet, build a list of extruded 4D simplices
/// @param[in] func: the implicit function that represents the swept volume. The input of the function is the 4d coordinate, and the output is an size-4 vector with first entry as the value and the other three as the gradient.
/// @param[in] maxTimeDep: maximum interger-valued time depth of the trajectory. Default: 1024
///
/// @param[out] timeList: a list of time stamps at this vertex
void init5CGrid(const int timeDep, mtet::MTetMesh grid, const std::function<std::pair<Scalar, Eigen::RowVector4d>(Eigen::RowVector4d)> func, const int maxTimeDep, vertExtrude &vertexMap){
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
        col.vert4dList = vertColList;
        vertexMap[value_of(vid)] = col;
    });
    grid.seq_foreach_tet([&](mtet::TetId tid, [[maybe_unused]] std::span<const mtet::VertexId, 4> data){
        std::span<VertexId, 4> vs = grid.get_tet(tid);
        for (size_t i = 0; i < 4; i++){
            vertexMap[value_of(vs[i])].vertTetAssoc.push_back(tid);
        }
    });
}

bool sign(double x) {
    return x > 0;
}

bool gridRefine(mtet::MTetMesh &grid, vertExtrude &vertexMap, insidenessMap &insideMap, const std::function<std::pair<Scalar, Eigen::RowVector4d>(Eigen::RowVector4d)> func, const double threshold, const int max_splits, std::array<double, timer_amount>& profileTimer){
    init5CGrid(3, grid, func, 1024, vertexMap);
    ///
    /// Initiate queue: timeQ and spaceQ
    auto compTime = [](std::tuple<mtet::Scalar, mtet::TetId, mtet::VertexId, int> timeSub0,
                       std::tuple<mtet::Scalar, mtet::TetId, mtet::VertexId, int> timeSub1)
    { return std::get<0>(timeSub0) < std::get<0>(timeSub1); };
    std::vector<std::tuple<mtet::Scalar, mtet::TetId, mtet::VertexId, int>> timeQ;
    
    auto compSpace = [](std::tuple<mtet::Scalar, mtet::TetId, mtet::EdgeId> spaceSub0,
                        std::tuple<mtet::Scalar, mtet::TetId, mtet::EdgeId> spaceSub1)
    { return std::get<0>(spaceSub0) < std::get<0>(spaceSub1); };
    std::vector<std::tuple<mtet::Scalar, mtet::TetId, mtet::EdgeId>> spaceQ;
    
    
    auto queryCell5 = [&](cell5 simp, std::array<vertexCol, 4> vs){
        
    };
    int splits = 0, temporal_splits = 0, spatial_splits = 0;
    ///
    /// Push Queue Function:
    auto push_one_col = [&](mtet::TetId tid)
    {
        std::span<VertexId, 4> vs = grid.get_tet(tid);
        //        simpCol colInfo = cell5Map[vs];
        //        simpCol::cell5_list cell5Col = colInfo.cell5Col;
        simpCol::cell5_list cell5Col = sampleCol(vs, vertexMap);
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
        std::vector<int> timeLenList(cell5Col.size());
        std::vector<mtet::Scalar> timeList(cell5Col.size());
        std::vector<size_t> indList(cell5Col.size());
        std::vector<bool> subList(cell5Col.size(), false);
        std::vector<bool> choiceList(cell5Col.size());
        Timer push_col(eval_tet_col, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        // tried parallel evaluation for each 4d tet in the column. Failed due to the need for early termination. May try in the future.
        //        dr::parallel_for(
        //                         dr::blocked_range<size_t>(0, cell5Col.size(), cell5Col.size() / 12),
        //                         [&](dr::blocked_range<size_t> range) {
        //                             for (auto cell5It : range) {
        //                                 auto simp = cell5Col[cell5It];
        //                                 std::array<int, 5> cell5Index = simp.hash;
        //                                 int lastInd = cell5Index[4];
        //                                 std::array<vertex4d, 5> verts;
        //                                 verts[0] = baseVerts[lastInd].vert4dList[cell5Index[lastInd]];
        //                                 size_t ind = 0;
        //                                 for (size_t i = 0; i < 4; i++){
        //                                     if (i != lastInd){
        //                                         ind ++;
        //                                         verts[ind] = baseVerts[i].vert4dList[cell5Index[i]];
        //                                     }
        //                                 }
        //                                 verts[4] = baseVerts[lastInd].vert4dList[cell5Index[lastInd] - 1];
        //                                 bool inside = false;
        //                                 bool choice = false;
        //                                 bool ret = refineContour(verts, threshold, inside, choice, profileTimer, spatial_splits);
        //                                 if (inside) {
        //                                     cell5Map[vs].covered = true;
        //                                     push_col.Stop();
        //                                     return;
        //                                 }
        //                                 if (ret) {
        //                                     subList[cell5It] = true;
        //                                     timeLenList[cell5It] = verts[0].time - verts[4].time;
        //                                     timeList[cell5It] = (verts[0].time + verts[4].time) / 2;
        //                                     indList[cell5It] = lastInd;
        //                                     choiceList[cell5It] = choice;
        //                                 }
        //                             }
        //                         });
        //        push_col.Stop();
        //Timer push_col(eval_tet_col, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        for (size_t cell5It = 0; cell5It < cell5Col.size(); cell5It++){
            auto simp = cell5Col[cell5It];
            std::array<int, 5> cell5Index = simp.hash;
            int lastInd = cell5Index[4];
            std::array<vertex4d, 5> verts;
            verts[0] = baseVerts[lastInd].vert4dList[cell5Index[lastInd]];
            size_t ind = 0;
            for (size_t i = 0; i < 4; i++){
                if (i != lastInd){
                    ind ++;
                    verts[ind] = baseVerts[i].vert4dList[cell5Index[i]];
                }
            }
            verts[4] = baseVerts[lastInd].vert4dList[cell5Index[lastInd] - 1];
            bool inside = false;
            bool choice = false;
            //            Timer ref_crit_timer(ref_crit, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            bool ret = refineContour(verts, threshold, inside, choice, profileTimer);
            //            ref_crit_timer.Stop();
            if (inside) {
                insideMap[vs] = true;
                //push_col.Stop();
                return;
            }
            //choice = (verts[0].coord[3] - verts[4].coord[3]) > longest_edge_length;
            if (ret) {
                subList[cell5It] = true;
                timeLenList[cell5It] = verts[0].time - verts[4].time;
                timeList[cell5It] = (verts[0].time + verts[4].time) / 2;
                indList[cell5It] = lastInd;
                choiceList[cell5It] = choice;
            }
        }
        //push_col.Stop();
        for (size_t cell5It = 0; cell5It < cell5Col.size(); cell5It++){
            if (subList[cell5It]){
                if (choiceList[cell5It]){
                    if (timeLenList[cell5It] > 4){
                        timeQ.emplace_back(timeLenList[cell5It], tid, vs[indList[cell5It]], timeList[cell5It]);
                        std::push_heap(timeQ.begin(), timeQ.end(), compTime);
                    }
                }else{
                    if (!baseSub){
                        spaceQ.emplace_back(longest_edge_length, tid, longest_edge);
                        std::push_heap(spaceQ.begin(), spaceQ.end(), compSpace);
                        baseSub = true;
                    }
                }
                
            }
        }
        if (!baseSub){
            std::array<vertex4d, 4> verts;
            for (size_t i = 0; i < 4; i++){
                verts[i] = baseVerts[i].vert4dList.front();
            }
            bool zeroX = false;
            bool ret = refineCap(verts, threshold * 0.25, zeroX);
            if (ret){
                spaceQ.emplace_back(longest_edge_length, tid, longest_edge);
                std::push_heap(spaceQ.begin(), spaceQ.end(), compSpace);
                baseSub = true;
            }
        }
        if (!baseSub){
            std::array<vertex4d, 4> verts;
            for (size_t i = 0; i < 4; i++){
                verts[i] = baseVerts[i].vert4dList.back();
            }
            bool zeroX = false;
            bool ret = refineCap(verts, threshold * 0.25, zeroX);
            if (ret){
                spaceQ.emplace_back(longest_edge_length, tid, longest_edge);
                std::push_heap(spaceQ.begin(), spaceQ.end(), compSpace);
                baseSub = true;
            }
        }
    };
    
    auto push_simps = [&](mtet::TetId tid, simpCol::cell5_list cell5Col, llvm_vecsmall::SmallVector<size_t, 256> refineList)
    {
        std::span<VertexId, 4> vs = grid.get_tet(tid);
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
        std::vector<int> timeLenList(refineList.size());
        std::vector<mtet::Scalar> timeList(refineList.size());
        std::vector<size_t> indList(refineList.size());
        std::vector<bool> subList(refineList.size(), false);
        std::vector<bool> choiceList(refineList.size());
        for (size_t it = 0; it < refineList.size(); it++){
            size_t cell5It = refineList[it];
            auto simp = cell5Col[cell5It];
            std::array<int, 5> cell5Index = simp.hash;
            int lastInd = cell5Index[4];
            std::array<vertex4d, 5> verts;
            verts[0] = baseVerts[lastInd].vert4dList[cell5Index[lastInd]];
            size_t ind = 0;
            for (size_t i = 0; i < 4; i++){
                if (i != lastInd){
                    ind ++;
                    verts[ind] = baseVerts[i].vert4dList[cell5Index[i]];
                }
            }
            verts[4] = baseVerts[lastInd].vert4dList[cell5Index[lastInd] - 1];
            bool inside = false;
            bool choice = false;
            bool ret = refineContour(verts, threshold, inside, choice, profileTimer);
            if (inside) {
                insideMap[vs] = true;
                return;
            }
            if (ret) {
                subList[it] = true;
                timeLenList[it] = verts[0].time - verts[4].time;
                timeList[it] = (verts[0].time + verts[4].time) / 2;
                indList[it] = lastInd;
                choiceList[it] = choice;
            }
            
        }
        for (size_t it = 0; it < refineList.size(); it++){
            size_t cell5It = refineList[it];
            if (subList[it]){
                if (choiceList[it]){
                    if (timeLenList[it] > 4){
                        timeQ.emplace_back(timeLenList[it], tid, vs[indList[it]], timeList[it]);
                        std::push_heap(timeQ.begin(), timeQ.end(), compTime);
                    }
                }else{
                    if (!baseSub){
                        spaceQ.emplace_back(longest_edge_length, tid, longest_edge);
                        std::push_heap(spaceQ.begin(), spaceQ.end(), compSpace);
                        baseSub = true;
                    }
                }
            }
        }
    };
    
    grid.seq_foreach_tet([&](mtet::TetId tid, [[maybe_unused]] std::span<const mtet::VertexId, 4> vs)
                         { push_one_col(tid); });
    while ((!timeQ.empty() || !spaceQ.empty()) && splits < max_splits){
        if (!timeQ.empty()){
            /// temporal subdivision:
            Timer temporal_splits_timer(time_splits, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            std::pop_heap(timeQ.begin(), timeQ.end(), compTime);
            auto [_ , tid, vid, time] = timeQ.back();
            timeQ.pop_back();
            
            if (!grid.has_tet(tid)){
                continue;
            }
            std::span<VertexId, 4> vs = grid.get_tet(tid);
            vertexCol timeList = vertexMap[value_of(vid)];
            if (insideMap[vs] || timeList.timeExist[time]){
                continue;
            }
            splits ++;
            temporal_splits++;
            timeList.timeExist[time] = true;
            vertexCol::vertToSimp newTets;
            newTets.reserve(timeList.vertTetAssoc.size());
            Timer add_assoc(add_connectivity, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            for (size_t i = 0; i < timeList.vertTetAssoc.size(); i++){
                auto assocTet = timeList.vertTetAssoc[i];
                if (grid.has_tet(assocTet)){
                    newTets.emplace_back(assocTet);
                }
            }
            add_assoc.Stop();
            timeList.vertTetAssoc = newTets;
            vertex4d newVert;
            newVert.time = time;
            newVert.coord = {timeList.vert4dList[0].coord[0], timeList.vert4dList[0].coord[1], timeList.vert4dList[0].coord[2], (double)time / 1024};
            newVert.inherit = false;
            newVert.valGradList = func(newVert.coord);
            timeList.insertTime(newVert);
            vertexMap[value_of(vid)] = timeList;
            for (size_t i = 0; i < newTets.size(); i++){
                std::span<VertexId, 4> newVs = grid.get_tet(newTets[i]);
                std::array<uint64_t, 4> quad = {value_of(newVs[0]), value_of(newVs[1]), value_of(newVs[2]), value_of(newVs[3])};
                size_t ind = static_cast<int>(std::distance(quad.begin(), std::find(quad.begin(), quad.end(), value_of(vid))));
                simpCol::cell5_list cell5Col = sampleCol(newVs, vertexMap);
                llvm_vecsmall::SmallVector<size_t, 256> refineList = resampleTimeCol(cell5Col, time, ind);
                push_simps(newTets[i], cell5Col, refineList);
            }
            temporal_splits_timer.Stop();
        }else{
            /// spatial subdivision:
            Timer spatial_splits_timer(space_splits, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            std::pop_heap(spaceQ.begin(), spaceQ.end(), compSpace);
            auto [_ , tid, eid] = spaceQ.back();
            spaceQ.pop_back();
            if (!grid.has_tet(tid)){
                continue;
            }
            std::span<VertexId, 4> vs = grid.get_tet(tid);
            if (insideMap[vs]){
                continue;
            }
            Timer create_vert(init_vert_col, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            splits ++;
            spatial_splits++;
            vertexCol newVert;
            std::array<VertexId, 2> vs_old = grid.get_edge_vertices(eid);
            auto [vid, eid0, eid1] = grid.split_edge(eid);
            std::span<Scalar, 3> vidCoord = grid.get_vertex(vid);
            vertexCol v0Col = vertexMap[value_of(vs_old[0])];
            vertexCol v1Col = vertexMap[value_of(vs_old[1])];
            vertexCol newVertCol;
            vertexCol::time_list v0Time = v0Col.getTimeList();
            vertexCol::time_list v1Time = v1Col.getTimeList();
            vertexCol::time_list tSamples;
            std::set_union(v0Time.begin(), v0Time.end(), v1Time.begin(), v1Time.end(), std::back_inserter(tSamples));
            vertexCol::vert4d_list vertColList(tSamples.size());
            for (int i = 0; i < tSamples.size(); i++){
                vertex4d vert;
                vert.time = tSamples[i];
                double time_fp = (double)vert.time / 1024;
                vert.coord = {vidCoord[0], vidCoord[1], vidCoord[2], time_fp};
                vert.valGradList = func(vert.coord);
                vertColList[i] = vert;
            }
            newVertCol.vert4dList = vertColList;
            vertexMap[value_of(vid)] = newVertCol;
            create_vert.Stop();
            grid.foreach_tet_around_edge(eid0, [&](mtet::TetId t0)
                                         {
                std::span<VertexId, 4> vs = grid.get_tet(t0);
                Timer add_assoc(add_connectivity, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
                for (size_t i = 0; i < 4; i++){
                    vertexMap[value_of(vs[i])].vertTetAssoc.emplace_back(t0);
                }
                add_assoc.Stop();
                Timer add_col_timer(extrude_tet_col, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
                insideMap[vs] = false;
                add_col_timer.Stop();
                push_one_col(t0);
            });
            grid.foreach_tet_around_edge(eid1, [&](mtet::TetId t1)
                                         {
                std::span<VertexId, 4> vs = grid.get_tet(t1);
                Timer add_assoc(add_connectivity, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
                for (size_t i = 0; i < 4; i++){
                    vertexMap[value_of(vs[i])].vertTetAssoc.emplace_back(t1);
                }
                add_assoc.Stop();
                Timer add_col_timer(extrude_tet_col, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
                insideMap[vs] = false;
                add_col_timer.Stop();
                push_one_col(t1);
            });
            spatial_splits_timer.Stop();
        }
    }
    std::cout << "Temporal splits: " << temporal_splits << " Spatial splits " << spatial_splits << std::endl;
    return true;
}
