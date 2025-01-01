//
//  adaptive_column_grid.h
//  adaptive_column_grid
//
//  Created by Yiwen Ju on 12/3/24.
//

#ifndef adaptive_column_grid_h
#define adaptive_column_grid_h
#include <iostream>
#include <nlohmann/json.hpp>
#include <ankerl/unordered_dense.h>
#include <SmallVector.h>
#include <fstream>
#include <mtet/mtet.h>
#include <mtet/io.h>
#include <math.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>

using namespace mtet;

class vertex4d{
public:
    int time; //int-valued hash; Default largest timestamp is 1024
    Eigen::RowVector4d coord;
    std::pair<Scalar, Eigen::RowVector4d> valGradList;
    int active_cells_num = 0;
    bool eval;
    vertex4d() = default;
    
    
    
    bool isActive(){
        return eval && ((active_cells_num > 0) || (valGradList.first <= 0)) ;
    }
    
    void printActive(){
        std::cout << active_cells_num << " " << valGradList.first << std::endl;
    }
};

auto compVertex = [](vertex4d v0, vertex4d v1){
    return v0.time < v1.time;};

class cell5{
public:
    std::array<int, 5> hash;
    std::array<int, 5> time_list;
    int level = 0;
    
    cell5() = default;
    
    int top(){
        return time_list[hash[4]];
    }
    
    int bot(int i){
        if (i == hash[4]){
            return time_list[4];
        }else{
            return time_list[i];
        }
    }
    
    std::shared_ptr<cell5> rebuildCell5(const int time, const int ind){
        cell5 simp;

        std::array<int, 5> botInd = {hash[0], hash[1], hash[2], hash[3], 0};
        int i = bot(0), j = bot(1), k = bot(2), l = bot(3);
        botInd[hash[4]] --;
        switch (ind) {
            case 0:
                botInd[0]++;
                botInd[4] = 0;
                simp.time_list = {time, j, k, l, i};
                break;
            case 1:
                botInd[1]++;
                botInd[4] = 1;
                simp.time_list = {i, time, k, l, j};
                break;
            case 2:
                botInd[2]++;
                botInd[4] = 2;
                simp.time_list = {i, j, time, l, k};
                break;
            case 3:
                botInd[3]++;
                botInd[4] = 3;
                simp.time_list = {i, j, k, time, l};
                break;
        }
        simp.hash = botInd;
        return std::make_shared<cell5>(simp);
    }
    
    cell5 copyCell5(){
        return *this;
    }
};

class vertexCol{
public:
    using vert4d_list = llvm_vecsmall::SmallVector<vertex4d, 256>;
    using time_list = llvm_vecsmall::SmallVector<int, 256>;
    using coord_list = llvm_vecsmall::SmallVector<Eigen::RowVector4d, 256>;
    using vertToSimp = llvm_vecsmall::SmallVector<mtet::TetId, 256>;
    
    vert4d_list timeStamp;
    ankerl::unordered_dense::map<int, bool> timeExist;
    vertToSimp vertTetAssoc;
    
    vertexCol() = default;
    
    int insertTime(vertex4d& newVert) {
        // Find the position to insert using binary search
        auto it = std::lower_bound(timeStamp.begin(), timeStamp.end(), newVert, compVertex);
        int ind = std::distance(timeStamp.begin(), it);
        timeStamp.insert(it, newVert);
        return ind;
    }
    
    time_list getTimeList(){
        time_list timeList(timeStamp.size());
        for (size_t i = 0; i < timeStamp.size(); i++){
            timeList[i] = timeStamp[i].time;
        }
        return timeList;
    }
    
    coord_list getCoordList(){
        coord_list coordList(timeStamp.size());
        for (size_t i = 0; i < timeStamp.size(); i++){
            coordList[i] = timeStamp[i].coord;
        }
        return coordList;
    }
    
    void addTetAssoc(llvm_vecsmall::SmallVector<int, 256> active_list){
        assert(active_list.size() == timeStamp.size());
        for (size_t i = 0; i < active_list.size(); i++){
            if (active_list[i] > 0){
                timeStamp[i].active_cells_num ++;
            }
        }
    }
    
    void delTetAssoc(llvm_vecsmall::SmallVector<int, 256> active_list){
        //std::cout << active_list.size() << " " << timeStamp.size() << std::endl;
        assert(active_list.size() == timeStamp.size());
        for (size_t i = 0; i < active_list.size(); i++){
            if (active_list[i] > 0){
                timeStamp[i].active_cells_num --;
            }
        }
    }
    
    void printActiveTag(){
        std::cout << "zeroX Tag: " << std::endl;
            for (size_t j = 0; j < timeStamp.size(); j++){
                std::cout << timeStamp[j].active_cells_num << " ";
            }
            std::cout << std::endl;
    }
};

class simpCol{
public:
    using cell5_list = llvm_vecsmall::SmallVector<std::shared_ptr<cell5>, 256>;
    using active_tag_list = llvm_vecsmall::SmallVector<int, 256>;
    using cell5_bool_list = llvm_vecsmall::SmallVector<bool, 256>;
    cell5_list cell5Col;
    int level = 0;// to prevent a subdivision of cell5 that's already been refined spatially
    bool covered = false;
    cell5_bool_list zeroX_cell;
    cell5_bool_list cell5_eval;
    std::array<active_tag_list, 4> vertex_active_tag; // to keep track of active-ness of each 4D vertex. This needs to be dynamically updated during every subdivision.
    simpCol() = default;
    
    void active_tag_plus(size_t cell5It){
        zeroX_cell[cell5It] = true;
        std::array<int, 5> cell5Index = cell5Col[cell5It]->hash;
        for (size_t i = 0; i < 4; i++){
            vertex_active_tag[i][cell5Index[i]]++;
        }
        int lastInd = cell5Index[4];
        vertex_active_tag[lastInd][cell5Index[lastInd] - 1]++;
    }
    
    void active_tag_minus(size_t cell5It){
        if (zeroX_cell[cell5It]){
            std::array<int, 5> cell5Index = cell5Col[cell5It]->hash;
            for (size_t i = 0; i < 4; i++){
                vertex_active_tag[i][cell5Index[i]]--;
            }
            int lastInd = cell5Index[4];
            vertex_active_tag[lastInd][cell5Index[lastInd] - 1]--;
        }
        zeroX_cell[cell5It] = false;
    }
    
    void insertTimeVert(size_t ind, int inserted_ind){
        auto front = vertex_active_tag[ind].begin();
        vertex_active_tag[ind].insert(front + inserted_ind, 0);
    }
    
    bool check_tag(){
        int active_cell_sum = 0;
        for (size_t i = 0; i < zeroX_cell.size(); i++){
            if (zeroX_cell[i]){
                active_cell_sum++;
            }
        }
        int active_vert_sum = 0;
        for (size_t i = 0; i < 4; i++){
            for (size_t j = 0; j < vertex_active_tag[i].size(); j++){
                active_vert_sum += vertex_active_tag[i][j];
            }
        }
        return active_cell_sum * 5 == active_vert_sum;
    }
    
    void printActiveTag(){
        for (size_t j = 0; j < zeroX_cell.size(); j++){
            std::cout << zeroX_cell[j] << " ";
        }
        std::cout << "zeroX cells" << std::endl;
        for (size_t i = 0; i < 4; i++){
            for (size_t j = 0; j < vertex_active_tag[i].size(); j++){
                std::cout << vertex_active_tag[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
    
    void printActiveTag(const size_t cell5It){
        std::array<int, 5> cell5Index = cell5Col[cell5It]->hash;
        for (size_t i = 0; i < 4; i++){
            std::cout << vertex_active_tag[i][cell5Index[i]] << " ";
        }
        int lastInd = cell5Index[4];
        std::cout << vertex_active_tag[lastInd][cell5Index[lastInd] - 1] << std::endl;
    }
};

using vertExtrude = ankerl::unordered_dense::map<uint64_t, vertexCol/*std::vector<vertex4d>*/>;

/// First, hash four tet vertices into a `uint64_t`
/// Since the tetid isn't const during the process, mount the boolean using vertexids of 4 corners.
struct TetHash
{
    using is_avalanching = void;
    using is_transparent = void;
    auto operator()(std::span<mtet::VertexId, 4> const& x) const noexcept -> uint64_t {
        ankerl::unordered_dense::hash<uint64_t> hash_fn;
        return ankerl::unordered_dense::detail::wyhash::hash(hash_fn(value_of(x[0])) + hash_fn(value_of(x[1])) + hash_fn(value_of(x[2])) + hash_fn(value_of(x[3])));
    }
};

/// Determine if a tet's hash is equal to another by comparing their vertices.
/// Two tet's vertex ids should be identical as each tet has a unique map to its tet vertices according to `mTet`.
struct TetEqual
{
    using is_transparent = void;
    bool operator()(std::span<mtet::VertexId, 4> const& lhs, std::span<mtet::VertexId, 4> const& rhs)
        const noexcept
    {
        return lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2] && lhs[3] == rhs[3];
    }
};

using tetExtrude = ankerl::unordered_dense::map<std::span<mtet::VertexId, 4>, simpCol, TetHash, TetEqual>;

#endif /* adaptive_column_grid_h */
