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

/// A 4D vertex, each is equipped with the following parameters:
/// @param time: an integer-valued time stamp. The default max of this value is 1024. The fourth coordinate of this vertex is a floating point of this value divided by 1024.
/// @param coord: The 4D coordinate of this vertex
/// @param valGradList: The value and gradients of this vertex given by the implicit function that represents the sweeping object.
/// @param active_cells_num: the number of 4D simplex uses this vertex
/// @param eval: If this vertex is evaluated
class vertex4d{
public:
    int time; //int-valued hash; Default largest vert4dList is 1024
    Eigen::RowVector4d coord;
    std::pair<Scalar, Eigen::RowVector4d> valGradList;
    bool inherit = true;
    vertex4d() = default;
};

auto compVertex = [](vertex4d v0, vertex4d v1){
    return v0.time < v1.time;};

/// A 4D simplex(5-cell) . Each is equipped iwth the following parameters:
/// @param hash: The index of the 5 vertices. This indexing system has the following structured 5-cell representation:
/// First four indices are the time stamp indices for each of the vertices {ti, tj, tk, tl}
/// The 5th index shows which vertex is being extruded 0/1/2/3 (The time stamp contains the high index for this vertex, i.e., need to retrieve the last vertex by finding the lower end of the extruded coordinate by the following form: hash[hash[4]] - 1)
/// @param time_list: The integer time at each index
/// @param level: The number of times this simplex is being visited in temporal subdivision. This needs to be aligned with the `level` of its simplex column in order to be subdivided
class cell5{
public:
    std::array<int, 5> hash;
    std::array<int, 5> time_list;
    int level = 0;
    
    cell5() = default;
    
    /// To obtain the time stamp of the top vertex of the extruded edge
    int top(){
        return time_list[hash[4]];
    }
    
    /// To obatin the time stamp at the given index i. For the index being extruded, it will pick the bottom time stamp.
    int bot(size_t i){
        if (i == hash[4]){
            return time_list[4];
        }else{
            return time_list[i];
        }
    }
    
    /// Build a simplex that will be inserted to the simplex column after the temporal edge subdivision.
    /// @param[in] time: The inserted time during the temporal split
    /// @param[in] ind: The index of the inserted time sample at this column of vertices
    /// @return: A new 5-cell/simplex
    cell5 rebuildCell5(const int time, const int ind){
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
        return simp;
    }
};

/// A column of 4D vertices equipped with the following members:
/// @param vert4dList: a list of 4D vertex in this column, each follows the data structure of `vertex4d`
/// @param timeExist: A hashed map to check if a time stamp exists in this column already
/// @param vertTetAssoc: A map from a vertex to a list of 3D tet indices
class vertexCol{
public:
    using vert4d_list = llvm_vecsmall::SmallVector<vertex4d, 256>;
    using time_list = llvm_vecsmall::SmallVector<int, 256>;
    using time_list_f = llvm_vecsmall::SmallVector<double, 256>;
    using value_list = llvm_vecsmall::SmallVector<double, 256>;
    using vertToSimp = llvm_vecsmall::SmallVector<mtet::TetId, 256>;
    
    vert4d_list vert4dList;
    ankerl::unordered_dense::map<int, bool> timeExist;
    vertToSimp vertTetAssoc;
    
    vertexCol() = default;

    /// Given a new 4D vertex, insert it to this column
    /// @param[in] newVert: A new 4D vertex
    /// @return The inserted index
    void insertTime(vertex4d& newVert) {
        // Find the position to insert using binary search
        auto it = std::lower_bound(vert4dList.begin(), vert4dList.end(), newVert, compVertex);
        vert4dList.insert(it, newVert);
    }
    
    /// find a list of time stamps of the 4D vertices in this column
    /// @return a list of integer-valued time stamp
    time_list getTimeList(){
        time_list timeList(vert4dList.size());
        for (size_t i = 0; i < vert4dList.size(); i++){
            timeList[i] = vert4dList[i].time;
        }
        return timeList;
    }
    
    time_list_f getTimeList_f(){
            time_list_f timeList(vert4dList.size());
            for (size_t i = 0; i < vert4dList.size(); i++){
                timeList[i] = vert4dList[i].coord(3);
            }
            return timeList;
        }
        
        value_list getValueList(){
            value_list valList(vert4dList.size());
            for (size_t i = 0; i < vert4dList.size(); i++){
                valList[i] = vert4dList[i].valGradList.second[3];
            }
            return valList;
        }
    
    void sortTime(){
        std::sort(vert4dList.begin(), vert4dList.end(), compVertex);
    }
};

/// A column of 4D simplices/5-cells equipped with the following members:
/// @param cell5_list: A list of shared pointers of 4D sipmlics, each follows the data structure of `cell5`
/// @param level: The level of this 4D column. The level increases with the number of times a temporal subdivision happened within the column. This keeps track of subdivisions of edges from a still valid simplex.
/// @param covered: Whether or not this column is covered; Covered means that there is a lifetd 3D tetrahedra face has all its bezier orndates below 0.
class simpCol{
public:
    using cell5_list = llvm_vecsmall::SmallVector<cell5, 256>;
    cell5_list cell5Col;
    bool covered = false;
    simpCol() = default;
};

/// A mount of a list of 4D vertex column to the 3D vertex
using vertExtrude = ankerl::unordered_dense::map<uint64_t, vertexCol>;

/// First, hash four tet vertices into a `uint64_t`
/// Since the tetid isn't const during the process, mount the boolean using vertexids of 4 corners.
struct TetHash
{
    using is_avalanching = void;
    using is_transparent = void;
    [[nodiscard]] auto operator()(std::span<mtet::VertexId, 4> const& x) const noexcept -> uint64_t {
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
        return value_of(lhs[0]) == value_of(rhs[0]) && value_of(lhs[1]) == value_of(rhs[1]) && value_of(lhs[2]) == value_of(rhs[2]) && value_of(lhs[3]) == value_of(rhs[3]);
    }
};

/// A mount of a list of 4D simplices column
using tetExtrude = ankerl::unordered_dense::map<std::span<mtet::VertexId, 4>, simpCol, TetHash, TetEqual>;

/// A mount of a boolean tag to every 3D tet to represent if the column is marked as "inside" of the sweep
using insidenessMap = ankerl::unordered_dense::map<std::span<mtet::VertexId, 4>, bool, TetHash, TetEqual>;

#endif /* adaptive_column_grid_h */
