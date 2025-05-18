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

struct spatialEq {
    bool operator()(std::span<const mtet::Scalar,3> a,
                    std::span<const mtet::Scalar,3> b) const noexcept
    {
        return a[0]==b[0]
        && a[1]==b[1]
        && a[2]==b[2];
    }
};

void parse_vertices(const mtetcol::Contour<4>& contour,
                    std::vector<double>& contour_time,
                    std::vector<int>& contour_index,
                    std::vector<Eigen::RowVector4d>& contour_pos,
                    std::vector<mtet::Scalar>& spatial_verts)
{
    spatialEq eq;
    int spatial_ind = 0;
    const int dim = 3;
    std::span<const mtet::Scalar, dim> spatial_it{ spatial_verts.data(), dim };
    auto num_vertices = contour.get_num_vertices();
    for (int i = 0 ; i < num_vertices; i++){
        std::span<const mtet::Scalar, 4> pos = contour.get_vertex(i);
        Eigen::Map<const Eigen::RowVector4d> pos_map(pos.data());
        contour_time.push_back(pos[3]);
        contour_pos.push_back(Eigen::RowVector4d{pos[0], pos[1], pos[2], pos[3]});
        if (!eq(pos.subspan<0 , dim>(), spatial_it)){
            spatial_ind++;
            spatial_it = std::span<const mtet::Scalar, 3>(spatial_verts.data() + dim * spatial_ind, dim);
        }
        contour_index.push_back(spatial_ind);
    }
}

void parse_polyhedron(const mtetcol::Contour<4>& contour,
                      mtetcol::Index poly_id,
                      std::vector<mtetcol::Index>& vert_id)
{
    vert_id.clear();
    std::vector<mtetcol::Index> vert_ind_tf(contour.get_num_vertices(), false);
    int vt = 0;
    auto poly = contour.get_polyhedron(poly_id);
    for (auto ci : poly) {
        auto cycle = contour.get_cycle(index(ci));
        for (auto si : cycle) {
            mtetcol::Index seg_id = index(si);
            auto seg = contour.get_segment(seg_id);
            mtetcol::Index v0 = seg[0];
            mtetcol::Index v1 = seg[1];
            if (!vert_ind_tf[v0]){
                vert_ind_tf[v0] = true;
                vert_id.push_back(v0);
            }
            if (!vert_ind_tf[v1]){
                vert_ind_tf[v1] = true;
                vert_id.push_back(v1);
            }
        }
    }
}

void compare_time(const double tet_time,
                  const double poly_time,
                  bool& intersect,
                  int& sign){
    if (tet_time == poly_time){
        intersect = true;
    }else
        if (tet_time > poly_time){
            if (sign == -1){
                intersect = true;
        }
        sign = 1;
    }else{
        if (sign == 1){
            intersect = true;
        }
        sign = -1;
    }
}

std::vector<uint32_t> one_column_simp = {0, 1, 2, 3};

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
        bool terminate = false;;
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
        std::vector<bool> zeroX_list(cell5Col.size());
        std::vector<std::array<int, 5>> cell5_index_list(cell5Col.size());
        for (size_t cell5It = 0; cell5It < cell5Col.size(); cell5It++){
            auto simp = cell5Col[cell5It];
            std::array<int, 5> cell5Index = simp.hash;
            cell5_index_list[cell5It] = cell5Index;
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
            bool zeroX = false;
            bool ret = refineFt(verts, 0.01, inside, choice, zeroX, profileTimer);
            zeroX_list[cell5It] = zeroX;
            if (inside) {
                insideMap[vs] = true;
                return;
            }
            if (ret) {
                subList[cell5It] = true;
                timeLenList[cell5It] = verts[0].time - verts[4].time;
                timeList[cell5It] = (verts[0].time + verts[4].time) / 2;
                indList[cell5It] = lastInd;
                choiceList[cell5It] = choice;
            }
        }
        for (size_t cell5It = 0; cell5It < cell5Col.size(); cell5It++){
            if (subList[cell5It]){
                terminate = true;
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
        if (terminate) return;
        
        std::vector<mtet::Scalar> spatial_verts;
        spatial_verts.reserve(12); //3 (dim) x 4 (vertices)tuple
        for (const auto& corners : vs){
            std::span<const Scalar, 3> data = grid.get_vertex(corners);
            spatial_verts.emplace_back(static_cast<double>(data[0]));
            spatial_verts.emplace_back(static_cast<double>(data[1]));
            spatial_verts.emplace_back(static_cast<double>(data[2]));
        }
        mtetcol::SimplicialColumn<4> column;
        std::array<vertexCol::time_list_f, 4> time = {baseVerts[0].getTimeList_f(),
            baseVerts[1].getTimeList_f(),
            baseVerts[2].getTimeList_f(),
            baseVerts[3].getTimeList_f()};
        std::function<std::span<double>(size_t)> time_func = [&](size_t index)->std::span<double>{
            return time[index];
        };
        std::array<vertexCol::value_list, 4> values = {baseVerts[0].getValueList(),
            baseVerts[1].getValueList(),
            baseVerts[2].getValueList(),
            baseVerts[3].getValueList()};
        std::function<std::span<double>(size_t)> values_func = [&](size_t index)->std::span<double>{
            return values[index];
        };
        column.set_vertices(spatial_verts);
        column.set_simplices(one_column_simp);
        column.set_time_samples(time_func, values_func);
        auto contour = column.extract_contour(0.0, false);
        auto num_polyhedra = contour.get_num_polyhedra();
        auto num_vertices = contour.get_num_vertices();
        
        std::vector<double> contour_time;
        contour_time.reserve(num_vertices);
        std::vector<int> contour_index;
        contour_index.reserve(num_vertices);
        std::vector<Eigen::RowVector4d> contour_pos;
        contour_pos.reserve(num_vertices);
        parse_vertices(contour, contour_time, contour_index, contour_pos, spatial_verts);
        for (int i = 0; i < num_polyhedra; i++){
            bool simple = contour.is_polyhedron_regular(i);
            std::vector<mtetcol::Index> vert_id;
            vert_id.reserve(num_vertices);
            parse_polyhedron(contour, i, vert_id);
            if (simple) assert(vert_id.size() == 4);
            bool caps = false;
            for (auto& vi : vert_id){
                auto poly_time = contour_time[vi];
                if (poly_time  == 1 || poly_time == 0){
                    caps = true;
                }
            }
            if (caps){
                if (simple){
                    std::array<Eigen::RowVector4d, 4> pts;
                    Eigen::RowVector4d vals;
                    std::array<Eigen::RowVector4d, 4> grads;
                    for (int vi = 0; vi < vert_id.size(); vi++){
                        pts[vi] = contour_pos[vert_id[vi]];
                        auto val_grad = func(pts[vi]);
                        vals[vi] = val_grad.first;
                        grads[vi] = val_grad.second;
                    }
                    if (refine3D(pts, vals, grads, threshold)){
                        spaceQ.emplace_back(longest_edge_length, tid, longest_edge);
                        std::push_heap(spaceQ.begin(), spaceQ.end(), compSpace);
                        terminate = true;
                    }
                }
            }
            if (terminate) return;
            // start traversing the active 5-cell.
            for (size_t cell5It = 0; cell5It < cell5Col.size(); cell5It++){
                if (!zeroX_list[cell5It]){
                    continue;
                }
                int sign = 0;
                bool intersect = false;
                auto& hash = cell5_index_list[cell5It];
                for (auto& vi : vert_id){
                    auto& poly_time = contour_time[vi];
                    auto& poly_ind = contour_index[vi];
                    auto tet_time = time[poly_ind][hash[poly_ind]];
                    compare_time(tet_time, poly_time, intersect, sign);
                    if (poly_ind == hash[4]){
                        tet_time = time[poly_ind][hash[poly_ind] - 1];
                        compare_time(tet_time, poly_time, intersect, sign);
                    }
                }
                if (intersect){
                    if (simple){
                        std::array<Eigen::RowVector4d, 4> pts;
                        Eigen::RowVector4d vals;
                        std::array<Eigen::RowVector4d, 4> grads;
                        for (int vi = 0; vi < vert_id.size(); vi++){
                            pts[vi] = contour_pos[vert_id[vi]];
                            auto val_grad = func(pts[vi]);
                            vals[vi] = val_grad.first;
                            grads[vi] = val_grad.second;
                        }
                        if (refine3D(pts, vals, grads, threshold)){
                            spaceQ.emplace_back(longest_edge_length, tid, longest_edge);
                            std::push_heap(spaceQ.begin(), spaceQ.end(), compSpace);
                            terminate = true;
                        }
                    }
                    else if (longest_edge_length > threshold * 0.1){
                        spaceQ.emplace_back(longest_edge_length, tid, longest_edge);
                        std::push_heap(spaceQ.begin(), spaceQ.end(), compSpace);
                        terminate = true;
                    }
                }
                if (terminate) return;
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
            bool zeroX = false;
            bool ret = refineFt(verts, 0.01, inside, choice, zeroX, profileTimer);
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
            for (size_t i = 0; i < timeList.vertTetAssoc.size(); i++){
                auto assocTet = timeList.vertTetAssoc[i];
                if (grid.has_tet(assocTet)){
                    newTets.emplace_back(assocTet);
                }
            }
            timeList.vertTetAssoc = newTets;
            vertex4d newVert;
            newVert.time = time;
            newVert.coord = {timeList.vert4dList[0].coord[0], timeList.vert4dList[0].coord[1], timeList.vert4dList[0].coord[2], (double)time / 1024};
            newVert.valGradList = func(newVert.coord);
            timeList.insertTime(newVert);
            vertexMap[value_of(vid)] = timeList;
            for (size_t i = 0; i < newTets.size(); i++){
                push_one_col(newTets[i]);
            }
        }else{
            /// spatial subdivision:
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
            grid.foreach_tet_around_edge(eid0, [&](mtet::TetId t0)
                                         {
                std::span<VertexId, 4> vs = grid.get_tet(t0);
                for (size_t i = 0; i < 4; i++){
                    vertexMap[value_of(vs[i])].vertTetAssoc.emplace_back(t0);
                }
                insideMap[vs] = false;
                push_one_col(t0);
            });
            grid.foreach_tet_around_edge(eid1, [&](mtet::TetId t1)
                                         {
                std::span<VertexId, 4> vs = grid.get_tet(t1);
                for (size_t i = 0; i < 4; i++){
                    vertexMap[value_of(vs[i])].vertTetAssoc.emplace_back(t1);
                }
                insideMap[vs] = false;
                push_one_col(t1);
            });
        }
    }
    std::cout << "Temporal splits: " << temporal_splits << " Spatial splits " << spatial_splits << std::endl;
    return true;
}
