//
//  col_gridgen.cpp
//  adaptive_column_grid
//
//  Created by Yiwen Ju on 12/4/24.
//

#include "col_gridgen.h"

/// Sample the list of 5-cells based on the base tetrahedra and 4 lists of time samples at its vertices. The extrusion/sampling is based on lowest time stamp, the second lowest time stamp, and the vertex id comparing the four incremental time stamps at each vertex.
/// @param[in] grid: the base tetrahedra grid in `mtet` structure
/// @param[in] tid: the tetrahedra id that is going to be extruded
/// @param[in] vertexMap: the map from vertex to a list of time samples.
/// @return A list of cell5 elements.
simpCol sampleCol(std::span<mtet::VertexId, 4> vs, vertExtrude &vertexMap){
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
    simpCol cell5Col;
    cell5Col.vertex_active_tag = {llvm_vecsmall::SmallVector<int, 256>(ti.size() - 1, 0),
        llvm_vecsmall::SmallVector<int, 256>(tj.size() - 1, 0),
        llvm_vecsmall::SmallVector<int, 256>(tk.size() - 1, 0),
        llvm_vecsmall::SmallVector<int, 256>(tl.size() - 1, 0)};
    cell5Col.cell5Col.reserve(ti.size()+tj.size()+tk.size()+tl.size() - 8);
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
        cell5Col.cell5Col.emplace_back(std::make_shared<cell5>(simp));
        cell5Col.zeroX_cell.emplace_back(false);
        cell5Col.cell5_eval.emplace_back(false);
    }
    return cell5Col;
}

llvm_vecsmall::SmallVector<size_t, 256> resampleTimeCol(
//simpCol::cell5_list& simpInfo, const int level,
                                                        simpCol& simpCol, const int time, const std::array<uint64_t, 4> quad, const int ind){
    simpCol::cell5_list simpInfo = simpCol.cell5Col;
    int level = simpCol.level;
    llvm_vecsmall::SmallVector<size_t, 256> refineList;
    refineList.reserve(simpInfo.size());
    std::pair<int, uint64_t> insertTime = std::pair<int, uint64_t>{time, quad[ind]};
    size_t i = 0;
    size_t inserted_cell;
    while (i < simpInfo.size() && std::pair<int, uint64_t>{simpInfo[i]->top(), quad[simpInfo[i]->hash[4]]} < insertTime){
        simpInfo[i]->level = level;
        i++;
    }
    inserted_cell = i;
    simpInfo.insert(simpInfo.begin() + i, simpInfo[i]->rebuildCell5(time, ind));
    simpInfo[i]->level = level;
    refineList.emplace_back(i);
    i ++;
    while (i < simpInfo.size()){
        if( simpInfo[i]->bot(ind) < time){
            simpCol.active_tag_minus(i - 1);
            cell5 prevCell = simpInfo[i]->copyCell5();
            simpInfo[i] = std::make_shared<cell5>(prevCell);
            if (simpInfo[i]->hash[4] == ind){
                simpInfo[i]->time_list[4] = time;
            }else{
                simpInfo[i]->time_list[ind] = time;
            }
            refineList.emplace_back(i);
        }
        simpInfo[i]->hash[ind] ++;
        simpInfo[i]->level = level;
        i++;
    }
    simpCol.cell5Col = simpInfo;
    simpCol.zeroX_cell.insert(simpCol.zeroX_cell.begin()+inserted_cell, false);
    simpCol.cell5_eval.insert(simpCol.cell5_eval.begin()+inserted_cell, false);
    return refineList;
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
            vert.eval = true;
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
        cell5Map[vs] = sampleCol(vs, vertexMap);
        //cell5Map[vs] = col;
    });
}

bool gridRefine(mtet::MTetMesh &grid, vertExtrude &vertexMap, tetExtrude &cell5Map, const std::function<std::pair<Scalar, Eigen::RowVector4d>(Eigen::RowVector4d)> func, const double threshold, const int max_splits, std::array<double, timer_amount>& profileTimer){
    int eval = 0, not_eval = 0;
    init5CGrid(3, grid, func, 1024, vertexMap, cell5Map);
    ///
    /// Initiate queue: timeQ and spaceQ
    auto compTime = [](std::tuple<mtet::Scalar, mtet::TetId, mtet::VertexId, int, std::shared_ptr<cell5>> timeSub0,
                       std::tuple<mtet::Scalar, mtet::TetId, mtet::VertexId, int, std::shared_ptr<cell5>> timeSub1)
    { return std::get<0>(timeSub0) <= std::get<0>(timeSub1); };
    std::vector<std::tuple<mtet::Scalar, mtet::TetId, mtet::VertexId, int, std::shared_ptr<cell5>>> timeQ;
    
    auto compSpace = [](std::tuple<mtet::Scalar, mtet::TetId, mtet::EdgeId, std::shared_ptr<cell5>> spaceSub0,
                       std::tuple<mtet::Scalar, mtet::TetId, mtet::EdgeId, std::shared_ptr<cell5>> spaceSub1)
    { return std::get<0>(spaceSub0) <= std::get<0>(spaceSub1); };
    std::vector<std::tuple<mtet::Scalar, mtet::TetId, mtet::EdgeId, std::shared_ptr<cell5>>> spaceQ;
    
    
    auto queryCell5 = [&](cell5 simp, std::array<vertexCol, 4> vs){
        
    };
    int splits = 0;
    ///
    /// Push Queue Function:
    cell5 always_sub;
    auto always_sub_ptr = std::make_shared<cell5>(always_sub);
    
//    auto activate_vertex = [&](vertex4d& vert) -> vertex4d {
//        if (!vert.eval) {
//            vert.valGradList = func(vert.coord); // Modify the input
//            vert.eval = true;                   // Modify the input
//        }
//        return vert; // Return a copy of the modified object
//    };
    
    auto push_one_col = [&](mtet::TetId tid, const bool init)
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
            llvm_vecsmall::SmallVector<int, 256> timeLenList(cell5Col.size());
            llvm_vecsmall::SmallVector<mtet::Scalar, 256> timeList(cell5Col.size());
            llvm_vecsmall::SmallVector<size_t, 256> indList(cell5Col.size());
            llvm_vecsmall::SmallVector<bool, 256> subList(cell5Col.size(), false);
            llvm_vecsmall::SmallVector<bool, 256> choiceList(cell5Col.size());
            llvm_vecsmall::SmallVector<bool, 256> zeroXList = colInfo.cell5_eval;
            for (size_t cell5It = 0; cell5It < cell5Col.size(); cell5It++){
                if (zeroXList[cell5It] || init){
                    eval++;
                    auto simp = cell5Col[cell5It];
                    std::array<int, 5> cell5Index = simp->hash;
                    int lastInd = cell5Index[4];
                    std::array<vertex4d, 5> verts;
                    verts[0] = baseVerts[lastInd].timeStamp[cell5Index[lastInd]];
                    size_t ind = 0;
                    for (size_t i = 0; i < 4; i++){
                        if (i != lastInd){
                            ind ++;
                            verts[ind] = baseVerts[i].timeStamp[cell5Index[i]];
                            assert(verts[ind].eval);
                            assert(baseVerts[i].timeStamp[cell5Index[i]].eval);
                        }
                    }
                    verts[4] = baseVerts[lastInd].timeStamp[cell5Index[lastInd] - 1];
                    assert(verts[4].eval == true);
                    bool inside = false;
                    bool choice = false;
                    bool zeroX = false;
                    bool ret = refineContour(verts, threshold, inside, choice, zeroX);
                    if (inside) {
                        cell5Map[vs].covered = true;
                        return;
                    }
                    zeroXList[cell5It] = zeroX;
                    if (ret) {
                        subList[cell5It] = true;
                        timeLenList[cell5It] = verts[0].time - verts[4].time;
                        timeList[cell5It] = (verts[0].time + verts[4].time) / 2;
                        indList[cell5It] = lastInd;
                        choiceList[cell5It] = choice;
                    }
                }else{
                    not_eval++;
                }
            }
            for (size_t cell5It = 0; cell5It < cell5Col.size(); cell5It++){
                if (zeroXList[cell5It]){
                    //update active tag
                    cell5Map[vs].active_tag_plus(cell5It);
                }else{
                    cell5Map[vs].zeroX_cell[cell5It] = false;
                }
                if (subList[cell5It]){
                    assert(zeroXList[cell5It]);
                    if (choiceList[cell5It]){
                        if (timeLenList[cell5It] > 4){
                            timeQ.emplace_back(timeLenList[cell5It], tid, vs[indList[cell5It]], timeList[cell5It],
                                               cell5Col[cell5It]);//colInfo.level);
                            std::push_heap(timeQ.begin(), timeQ.end(), compTime);
                        }
                    }else{
                        if (!baseSub){
                            spaceQ.emplace_back(longest_edge_length, tid, longest_edge, cell5Col[cell5It]);
                            std::push_heap(spaceQ.begin(), spaceQ.end(), compSpace);
                            baseSub = true;
                        }
                    }
                    
                }
            }
            if (!baseSub){
                std::array<vertex4d, 4> verts;
                for (size_t i = 0; i < 4; i++){
                    verts[i] = baseVerts[i].timeStamp.front();
                }
                bool zeroX = false;
                bool ret = refineCap(verts, threshold/2, zeroX);
                if (ret){
                    spaceQ.emplace_back(longest_edge_length, tid, longest_edge, always_sub_ptr);
                    std::push_heap(spaceQ.begin(), spaceQ.end(), compSpace);
                    baseSub = true;
                }
            }
            if (!baseSub){
                std::array<vertex4d, 4> verts;
                for (size_t i = 0; i < 4; i++){
                    verts[i] = baseVerts[i].timeStamp.back();
                }
                bool zeroX = false;
                bool ret = refineCap(verts, threshold/2, zeroX);
                if (ret){
                    spaceQ.emplace_back(longest_edge_length, tid, longest_edge, always_sub_ptr);
                    std::push_heap(spaceQ.begin(), spaceQ.end(), compSpace);
                    baseSub = true;
                }
            }
            assert(cell5Map[vs].check_tag());
        };
        
    auto push_simps = [&](mtet::TetId tid, llvm_vecsmall::SmallVector<size_t, 256> refineList)
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
        llvm_vecsmall::SmallVector<int, 256> timeLenList(cell5Col.size());
        llvm_vecsmall::SmallVector<mtet::Scalar, 256> timeList(cell5Col.size());
        llvm_vecsmall::SmallVector<size_t, 256> indList(cell5Col.size());
        llvm_vecsmall::SmallVector<bool, 256> subList(cell5Col.size(), false);
        llvm_vecsmall::SmallVector<bool, 256> choiceList(cell5Col.size());
        llvm_vecsmall::SmallVector<bool, 256> zeroXList = colInfo.zeroX_cell;
        for (size_t it = 0; it < refineList.size(); it++){
            size_t cell5It = refineList[it];
            auto simp = cell5Col[cell5It];
            std::array<int, 5> cell5Index = simp->hash;
            int lastInd = cell5Index[4];
            std::array<vertex4d, 5> verts;
            verts[0] = baseVerts[lastInd].timeStamp[cell5Index[lastInd]];
            size_t ind = 0;
            for (size_t i = 0; i < 4; i++){
                if (i != lastInd){
                    ind ++;
                    verts[ind] = baseVerts[i].timeStamp[cell5Index[i]];
                    assert(verts[ind].eval);
                    assert(baseVerts[i].timeStamp[cell5Index[i]].eval);
                }
            }
            verts[4] = baseVerts[lastInd].timeStamp[cell5Index[lastInd] - 1];
            assert(zeroXList[cell5It] == false);
            bool inside = false;
            bool choice = false;
            bool zeroX = false;
            bool ret = refineContour(verts, threshold, inside, choice, zeroX);
            cell5Map[vs].cell5_eval[cell5It] = true;
            if (inside) {
                cell5Map[vs].covered = true;
                return;
            }
            zeroXList[cell5It] = zeroX;
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
            if (zeroXList[cell5It]){
                //update active tag
                cell5Map[vs].active_tag_plus(cell5It);
            }
            if (subList[it]){
                assert(zeroXList[cell5It]);
                if (choiceList[it]){
                    if (timeLenList[it] > 4){
                        timeQ.emplace_back(timeLenList[it], tid, vs[indList[it]], timeList[it],
                                           cell5Col[cell5It]);//colInfo.level);
                        std::push_heap(timeQ.begin(), timeQ.end(), compTime);
                    }
                }else{
                    if (!baseSub){
                        spaceQ.emplace_back(longest_edge_length, tid, longest_edge, cell5Col[cell5It]);
                        std::push_heap(spaceQ.begin(), spaceQ.end(), compSpace);
                        baseSub = true;
                    }
                }
            }
        }
        assert(cell5Map[vs].check_tag());
    };
    
    
    grid.seq_foreach_tet([&](mtet::TetId tid, [[maybe_unused]] std::span<const mtet::VertexId, 4> data)
                         {
        std::span<VertexId, 4> vs = grid.get_tet(tid);
        push_one_col(tid, true);
        for (size_t vertIt = 0; vertIt < 4; vertIt++){
            vertexMap[value_of(vs[vertIt])].addTetAssoc(cell5Map[vs].vertex_active_tag[vertIt]);
        }
    });
    while ((!timeQ.empty() || !spaceQ.empty()) && splits < max_splits){
        if (!timeQ.empty()){
            /// temporal subdivision:
            Timer temporal_splits_timer(time_splits, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            std::pop_heap(timeQ.begin(), timeQ.end(), compTime);
            auto [_ , tid, vid, time, simp] = timeQ.back();
            timeQ.pop_back();
            
            if (!grid.has_tet(tid)){
                continue;
            }
            std::span<VertexId, 4> vs = grid.get_tet(tid);
            vertexCol timeList = vertexMap[value_of(vid)];
            simpCol cell5List = cell5Map[vs];
            if (
                simp->level != cell5List.level || cell5List.covered || timeList.timeExist[time]){
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
            for (size_t i = 0; i < newTets.size(); i++){
                std::span<VertexId, 4> newVs = grid.get_tet(newTets[i]);
                if (cell5Map[newVs].cell5Col.size() == 0){
                    continue;
                }
                //udpate active tag on vertex
                for (size_t vertIt = 0; vertIt < 4; vertIt++){
                    vertexMap[value_of(newVs[vertIt])].delTetAssoc(cell5Map[newVs].vertex_active_tag[vertIt]);
                }
            }
            timeList.vertTetAssoc = newTets;
            vertex4d newVert;
            newVert.eval = true;
            newVert.time = time;
            newVert.coord = {timeList.timeStamp[0].coord[0], timeList.timeStamp[0].coord[1], timeList.timeStamp[0].coord[2], (double)time / 1024};
            newVert.valGradList = func(newVert.coord);
            int inserted_ind = timeList.insertTime(newVert);
            vertexMap[value_of(vid)] = timeList;

            for (size_t i = 0; i < newTets.size(); i++){
                std::span<VertexId, 4> newVs = grid.get_tet(newTets[i]);
                if (cell5Map[newVs].cell5Col.size() == 0){
                    continue;
                }
                cell5Map[newVs].level ++;
                std::array<uint64_t, 4> quad = {value_of(newVs[0]), value_of(newVs[1]), value_of(newVs[2]), value_of(newVs[3])};
                size_t ind = static_cast<int>(std::distance(quad.begin(), std::find(quad.begin(), quad.end(), value_of(vid))));
                llvm_vecsmall::SmallVector<size_t, 256> refineList =
                resampleTimeCol(cell5Map[newVs], time, quad, ind);
                cell5Map[newVs].insertTimeVert(ind, inserted_ind);
                push_simps(newTets[i], refineList);
                //udpate active tag on vertex
                for (size_t vertIt = 0; vertIt < 4; vertIt++){
                    vertexMap[value_of(newVs[vertIt])].addTetAssoc(cell5Map[newVs].vertex_active_tag[vertIt]);
                }
            }
            temporal_splits_timer.Stop();
        }else{
            /// spatial subdivision:
            Timer spatial_splits_timer(space_splits, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            std::pop_heap(spaceQ.begin(), spaceQ.end(), compSpace);
            auto [_ , tid, eid, simp] = spaceQ.back();
            spaceQ.pop_back();
            if (!grid.has_tet(tid)){
                continue;
            }
            std::span<VertexId, 4> vs = grid.get_tet(tid);
            simpCol cell5List = cell5Map[vs];
            if (cell5List.covered || (simp!= always_sub_ptr && simp->level != cell5List.level)){
                continue;
            }
            splits ++;
            vertexCol newVert;
            std::array<VertexId, 2> vs_old = grid.get_edge_vertices(eid);
            llvm_vecsmall::SmallVector<std::span<VertexId, 4>, 256> old_tets_buffer;
            grid.foreach_tet_around_edge(eid, [&](mtet::TetId oldTet)
                                         {
                std::span<VertexId, 4> vs = grid.get_tet(oldTet);
                //udpate active tag on vertex
                old_tets_buffer.push_back(vs);
            });
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
                vert.eval = true;
                double time_fp = (double)vert.time / 1024;
                vert.coord = {vidCoord[0], vidCoord[1], vidCoord[2], time_fp};
                vert.valGradList = func(vert.coord);
                vertColList[i] = vert;
            }
            newVertCol.timeStamp = vertColList;
            vertexMap[value_of(vid)] = newVertCol;
            llvm_vecsmall::SmallVector<bool, 256> evalList(tSamples.size(), false);
            
            auto activate_tet = [&](mtet::EdgeId sub_eid, mtet::VertexId vid) {
                grid.foreach_tet_around_edge(sub_eid, [&](mtet::TetId sub_tid)
                                             {
                    std::span<VertexId, 4> sub_vs = grid.get_tet(sub_tid);
                    std::array<vertexCol::vert4d_list, 4> vsCol;
                    for (size_t i = 0; i < 4; i++){
                        vertexMap[value_of(sub_vs[i])].vertTetAssoc.emplace_back(sub_tid);
                        vsCol[i] = vertexMap[value_of(sub_vs[i])].timeStamp;
                    }
                    simpCol col;
                    col = sampleCol(sub_vs, vertexMap);
                    cell5Map[sub_vs] = col;
                    std::array<uint64_t, 4> quad = {value_of(sub_vs[0]), value_of(sub_vs[1]), value_of(sub_vs[2]), value_of(sub_vs[3])};
                    size_t ind = std::distance(quad.begin(), std::find(quad.begin(), quad.end(), value_of(vid)));
                    for (size_t cell5It = 0; cell5It < col.cell5Col.size(); cell5It++){
                        bool refine = false;
                        std::array<int, 5> index = col.cell5Col[cell5It]->hash;
                        for (size_t vertIt = 0; vertIt < 4; vertIt++){
//                            if(!vsCol[vertIt][index[vertIt]].eval && !evalList[index[vertIt]]){
//                                evalList[index[vertIt]] = true;
//                            }
                            if (vsCol[vertIt][index[vertIt]].isActive()){
                                refine = true;
                                cell5Map[sub_vs].cell5_eval[cell5It] = true;
                                break;
                            }
                        }
//                        if(!vsCol[index[4]][index[index[4]] - 1].eval && !evalList[index[index[4]] - 1]){
//                            evalList[index[index[4]] - 1] = true;
//                        }
                        if (!refine){
                            if (vsCol[index[4]][index[index[4]] - 1].isActive()){
                                refine = true;
                                cell5Map[sub_vs].cell5_eval[cell5It] = true;
                            }
                        }
//                        if (refine){
//                            evalList[index[ind]] = true;
//                            if (index[4] == ind){
//                                evalList[index[ind] - 1] = true;
//                            }
//                        }
                    }
                });
            };
            activate_tet(eid0, vid);
            activate_tet(eid1, vid);
//            for (size_t i = 0; i < evalList.size(); i++){
//                if (evalList[i]){
//                    auto vert = vertexMap[value_of(vid)].timeStamp[i];
//                    vert.eval = true;
//                    vert.valGradList = func(vert.coord);
//                    vertexMap[value_of(vid)].timeStamp[i] = vert;
//                }
//            }
            grid.foreach_tet_around_edge(eid0, [&](mtet::TetId t0){
                push_one_col(t0, false);
                //udpate active tag on vertex
                for (size_t vertIt = 0; vertIt < 4; vertIt++){
                    vertexMap[value_of(vs[vertIt])].addTetAssoc(cell5Map[vs].vertex_active_tag[vertIt]);
                }
            });
            grid.foreach_tet_around_edge(eid1, [&](mtet::TetId t1){
                push_one_col(t1, false);
                //udpate active tag on vertex
                for (size_t vertIt = 0; vertIt < 4; vertIt++){
                    vertexMap[value_of(vs[vertIt])].addTetAssoc(cell5Map[vs].vertex_active_tag[vertIt]);
                }
            });
            for (size_t i = 0; i < old_tets_buffer.size(); i++){
                for (size_t vertIt = 0; vertIt < 4; vertIt++){
                    vertexMap[value_of(vs[vertIt])].delTetAssoc(cell5Map[vs].vertex_active_tag[vertIt]);
                }
            }
            spatial_splits_timer.Stop();
        }
    }
    std::cout << eval << " " << not_eval << " ";
    return true;
}
