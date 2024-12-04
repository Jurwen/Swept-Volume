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
#include <fstream>
#include <mtet/mtet.h>
#include <mtet/io.h>
#include <math.h>

using namespace mtet;

struct vertex4d{
    int time;
    //double coord; //not sure if its cost benefit
};

struct cell5{
    std::array<int, 5> hash; //same data structure as an element in cell5List(Mathematica)
    //std::pair<mtet::VertexId, int> vertices[5]; //not sure if its cost benefit
    //std::array<double, 4> coords[5]; //not sure if its cost benefit
};

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

using vertexCol = ankerl::unordered_dense::map<uint64_t, std::vector<vertex4d>>;
using tetCol = ankerl::unordered_dense::map<std::span<mtet::VertexId, 4>, std::vector<cell5>, TetHash, TetEqual>;

#endif /* adaptive_column_grid_h */
