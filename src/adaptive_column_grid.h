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

struct vertex4d{
    int time; //int-valued hash; Default largest timestamp is 1024, i.e., no 
    Eigen::RowVector4d coord;
    std::pair<Scalar, Eigen::RowVector4d> valGradList;
};

struct cell5{
    std::array<int, 5> hash; //same data structure as an element in cell5List(Mathematica)
    //std::pair<mtet::VertexId, int> vertices[5]; //not sure if its cost benefit
    //std::array<double, 4> coords[5]; //not sure if its cost benefit
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
};

struct simpCol{
    using cell5_list = llvm_vecsmall::SmallVector<cell5, 256>;
    cell5_list cell5Col;
    int level = 0;// to prevent a subdivision of cell5 that's already been refined spatially
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
