//
//  post_processing.h
//  general_sweep
//
//  Created by Yiwen Ju on 5/13/25.
//

#ifndef post_processing_h
#define post_processing_h
#include <ankerl/unordered_dense.h>
#include <arrangement/Arrangement.h>
#include <lagrange/SurfaceMesh.h>
#include <lagrange/mesh_cleanup/remove_isolated_vertices.h>
#include <lagrange/utils/SmallVector.h>
#include <lagrange/views.h>

#include <algorithm>

template <typename Scalar, typename Index>
lagrange::SurfaceMesh<Scalar, Index> isocontour_to_mesh(mtetcol::Contour<4>& isocontour)
{
    lagrange::SurfaceMesh<Scalar, Index> envelope;

    // Port the isocontour into lagrange mesh
    size_t num_vertices = isocontour.get_num_vertices();
    size_t num_cycles = isocontour.get_num_cycles();

    // Add vertices and time
    envelope.add_vertices(num_vertices);
    envelope.template create_attribute<double>(
        "time",
        lagrange::AttributeElement::Vertex,
        lagrange::AttributeUsage::Scalar,
        1);
    auto time_values = attribute_vector_ref<double>(envelope, "time");

    for (size_t i = 0; i < num_vertices; i++) {
        auto xyzt = isocontour.get_vertex(i);
        auto pos = envelope.ref_position(i);
        pos[0] = xyzt[0];
        pos[1] = xyzt[1];
        pos[2] = xyzt[2];
        time_values[i] = xyzt[3];
    }

    // Add polygons
    lagrange::SmallVector<uint32_t, 16> polygon;
    for (size_t i = 0; i < num_cycles; i++) {
        auto cycle = isocontour.get_cycle(i);
        size_t cycle_size = cycle.size();
        polygon.clear();
        polygon.resize(cycle_size);

        size_t ind = 0;
        for (auto si : cycle) {
            mtetcol::Index seg_id = index(si);
            bool seg_ori = mtetcol::orientation(si);
            auto seg = isocontour.get_segment(seg_id);
            polygon[ind] = (seg_ori ? seg[0] : seg[1]);
            std::pair<mtetcol::Index, mtetcol::Index> edge_key = {
                std::min(seg[0], seg[1]),
                std::max(seg[0], seg[1])};

            ind++;
        }
        envelope.add_polygon({polygon.data(), polygon.size()});
    }

    // Add regular attribute
    envelope.template create_attribute<uint8_t>(
        "regular",
        lagrange::AttributeElement::Facet,
        lagrange::AttributeUsage::Scalar,
        1);
    auto regular_values = attribute_vector_ref<uint8_t>(envelope, "regular");
    for (size_t i = 0; i < num_cycles; i++) {
        regular_values[i] = isocontour.is_cycle_regular(i) ? 1 : 0;
    }
    
    // Add polygon (triangle) orientation label from isocontour cycle labels
    envelope.template create_attribute<int8_t>(
        "orientation",
        lagrange::AttributeElement::Facet,
        lagrange::AttributeUsage::Scalar,
        1);
    auto orientation_values = attribute_vector_ref<int8_t>(envelope, "orientation");
    const auto& cycle_labels = isocontour.get_cycle_labels();
    if (cycle_labels.empty()) {
        orientation_values.setZero();
    } else {
        for (size_t i = 0; i < num_cycles; i++) {
            orientation_values[i] = isocontour.get_cycle_label(i);
        }
    }
    return envelope;
}

template <typename Scalar, typename Index>
lagrange::SurfaceMesh<Scalar, Index> compute_envelope_arrangement(
    const lagrange::SurfaceMesh<Scalar, Index>& envelope,
    Scalar vol_threshold,
    size_t face_count_threshold,
    bool cyclic)
{
    using Point = Eigen::Matrix<Scalar, 3, 1>;

    // Compute arrangement
    auto V = vertex_view(envelope).template cast<double>();
    auto F = facet_view(envelope).template cast<int>();
    auto T = attribute_vector_view<Scalar>(envelope, "time");
    auto O = attribute_vector_view<int8_t>(envelope, "orientation");
    arrangement::VectorI face_labels = Eigen::VectorXi::LinSpaced(F.rows(), 0, F.rows() - 1);
    auto engine = arrangement::Arrangement::create_mesh_arrangement(V, F, face_labels);
    engine->run();

    lagrange::SurfaceMesh<Scalar, Index> sweep_arrangement;

    const auto& cell_data = engine->get_cells(); // (#facets x 2) array
    const auto& patches = engine->get_patches(); // list of facet indices
    const auto& parent_facets = engine->get_out_face_labels(); // size = #facets
    const auto& winding_number = engine->get_winding_number(); // (#facets x 2) array

    const auto& out_vertices = engine->get_vertices();
    const auto& arrangement_faces = engine->get_faces();
    size_t num_cells = engine->get_num_cells();
    size_t num_patches = engine->get_num_patches();
    size_t num_facets = arrangement_faces.rows();
    assert(patches.size() == num_facets);
    assert(cell_data.rows() == num_patches);

    sweep_arrangement.add_vertices(
        static_cast<Index>(out_vertices.rows()),
        {out_vertices.data(), static_cast<size_t>(out_vertices.size())});
    sweep_arrangement.add_triangles(static_cast<Index>(num_facets));
    {
        auto facets = facet_ref(sweep_arrangement);
        std::copy_n(arrangement_faces.data(), arrangement_faces.size(), facets.data());
    }
    sweep_arrangement.template create_attribute<Index>(
        "envelope_facet_id",
        lagrange::AttributeElement::Facet,
        lagrange::AttributeUsage::Scalar,
        1);
    auto envelope_facet_id = attribute_vector_ref<Index>(sweep_arrangement, "envelope_facet_id");
    std::copy_n(parent_facets.data(), parent_facets.size(), envelope_facet_id.data());
    sweep_arrangement.initialize_edges();

    // Build cell adjacency graph and compute cell volumes
    constexpr int32_t invalid_winding_number = std::numeric_limits<int32_t>::max();
    std::vector<ankerl::unordered_dense::set<Index>> cell_graph(num_cells);
    std::vector<Scalar> cell_volumes(num_cells, 0);
    std::vector<size_t> cell_face_counts(num_cells, 0);
    std::vector<int32_t> cell_winding_numbers(num_cells, invalid_winding_number);
    std::vector<int> cell_label_sums(num_cells, 0);
    for (size_t fid = 0; fid < num_facets; fid++) {
        Index c0 = static_cast<Index>(cell_data(patches[fid], 0)); // Cell on the positive side
        Index c1 = static_cast<Index>(cell_data(patches[fid], 1)); // Cell on the negative side
        int w0 = winding_number(fid, 0); // Winding number on the positive side
        int w1 = winding_number(fid, 1); // Winding number on the negative side

        if (cell_winding_numbers[c0] == invalid_winding_number) {
            cell_winding_numbers[c0] = w0;
        } else {
            if (cell_winding_numbers[c0] != w0) {
                // This should never happen
                throw std::runtime_error("Inconsistent winding numbers detected!");
            }
        }
        if (cell_winding_numbers[c1] == invalid_winding_number) {
            cell_winding_numbers[c1] = w1;
        } else {
            if (cell_winding_numbers[c1] != w1) {
                // This should never happen
                throw std::runtime_error("Inconsistent winding numbers detected!");
            }
        }

        cell_graph[c0].insert(c1);
        cell_graph[c1].insert(c0);

        Eigen::Matrix<int, 1, 3> f = arrangement_faces.row(fid);
        Point p0 = out_vertices.row(f[0]).template cast<Scalar>();
        Point p1 = out_vertices.row(f[1]).template cast<Scalar>();
        Point p2 = out_vertices.row(f[2]).template cast<Scalar>();
        Scalar vol = p0.dot((p1).cross(p2));

        cell_volumes[c0] -= vol / 6;
        cell_volumes[c1] += vol / 6;
        cell_face_counts[c0]++;
        cell_face_counts[c1]++;

        // Accumulate orientation labels per cell
        int8_t ori_label = O[parent_facets(fid)];
        cell_label_sums[c0] += ori_label;
        cell_label_sums[c1] += ori_label;
    }
    auto cell_is_small = [&](Index cid) {
        return (
            std::abs(cell_volumes[cid]) < vol_threshold ||
            cell_face_counts[cid] < face_count_threshold);
    };

    std::vector<int> parent_cell(num_cells);
    std::iota(parent_cell.begin(), parent_cell.end(), 0);
    auto get_parent = [&](auto&& self, int cid) -> int {
        if (parent_cell[cid] != cid) {
            parent_cell[cid] = self(self, parent_cell[cid]);
        }
        return parent_cell[cid];
    };
    for (size_t cid = 0; cid < num_cells; cid++) {
        if (cell_is_small(cid)) {
            // Union small cell with one of its neighbors
            int parent = parent_cell[cid];
            Scalar max_vol = std::abs(cell_volumes[cid]);
            for (auto adj_cid : cell_graph[cid]) {
                if (std::abs(cell_volumes[adj_cid]) > max_vol) {
                    max_vol = std::abs(cell_volumes[adj_cid]);
                    parent = adj_cid;
                }
            }
            parent_cell[cid] = get_parent(get_parent, parent);
        }
    }

    // Compute sweep surface facets
    sweep_arrangement.template create_attribute<int8_t>(
        "valid",
        lagrange::AttributeElement::Facet,
        lagrange::AttributeUsage::Scalar,
        1);
    auto is_valid = attribute_vector_ref<int8_t>(sweep_arrangement, "valid");
    is_valid.setZero();
    size_t num_valid_facets = 0;
    for (size_t fid = 0; fid < num_facets; fid++) {
        int c0 = cell_data(patches[fid], 0); // Cell on the positive side
        int c1 = cell_data(patches[fid], 1); // Cell on the negative side
        c0 = get_parent(get_parent, c0); // Find the representative parent cell
        c1 = get_parent(get_parent, c1); // Find the representative parent cell
        int w0 = cell_winding_numbers[c0]; // Winding number on the positive side
        int w1 = cell_winding_numbers[c1]; // Winding number on the negative side
        int ls0 = cell_label_sums[c0]; // Label sum on the positive side
        int ls1 = cell_label_sums[c1]; // Label sum on the negative side
        
        // A cell is "outside" if its winding number is 0 AND its label sum >= 0.
        // A cell is "inside" (negation) if its winding number != 0 OR label sum < 0.
        bool c0_outside, c1_outside;
        if(cyclic){
            c0_outside = (w0 == 0 && ls0 >= 0);
            c1_outside = (w1 == 0 && ls1 >= 0);
        }else{
            c0_outside = (w0 == 0);
            c1_outside = (w1 == 0);
        }
        if (c0_outside && !c1_outside) {
            is_valid[fid] = 1;
            num_valid_facets++;
        } else if (c1_outside && !c0_outside) {
            is_valid[fid] = -1;
            num_valid_facets++;
        }
    }

    // Compute per-corner time attribute
    auto interpolate_time = [&](Index fid, Point p) {
        Eigen::Matrix<int, 1, 3> f = F.row(fid);
        Point p0 = V.row(f[0]);
        Point p1 = V.row(f[1]);
        Point p2 = V.row(f[2]);
        auto t0 = T[f[0]];
        auto t1 = T[f[1]];
        auto t2 = T[f[2]];

        Scalar b0 = ((p1 - p).cross(p2 - p)).norm();
        Scalar b1 = ((p2 - p).cross(p0 - p)).norm();
        Scalar b2 = ((p0 - p).cross(p1 - p)).norm();
        Scalar denom = b0 + b1 + b2;
        if (std::abs(denom) > std::numeric_limits<Scalar>::epsilon()) {
            return (t0 * b0 + t1 * b1 + t2 * b2) / denom;
        } else {
            // Simple average to avoid numerical issues
            return (t0 + t1 + t2) / Scalar(3);
        }
    };

    sweep_arrangement.template create_attribute<Scalar>(
        "time",
        lagrange::AttributeElement::Corner,
        lagrange::AttributeUsage::Scalar,
        1);
    auto time_values = attribute_vector_ref<Scalar>(sweep_arrangement, "time");
    Index count = 0;
    // TODO: parallelize this loop
    for (size_t fid = 0; fid < num_facets; fid++) {
        auto idx = count++;
        auto parent_fid = parent_facets(fid);
        Index v0 = static_cast<Index>(arrangement_faces(fid, 0));
        Index v1 = static_cast<Index>(arrangement_faces(fid, 1));
        Index v2 = static_cast<Index>(arrangement_faces(fid, 2));
        Point p0(out_vertices.row(v0).template cast<Scalar>());
        Point p1(out_vertices.row(v1).template cast<Scalar>());
        Point p2(out_vertices.row(v2).template cast<Scalar>());
        Scalar t0 = interpolate_time(parent_fid, p0);
        Scalar t1 = interpolate_time(parent_fid, p1);
        Scalar t2 = interpolate_time(parent_fid, p2);
        time_values[3 * idx + 0] = t0;
        time_values[3 * idx + 1] = t1;
        time_values[3 * idx + 2] = t2;
    }

    // Extract feature edges
    sweep_arrangement.template create_attribute<int8_t>(
        "is_feature",
        lagrange::AttributeElement::Edge,
        lagrange::AttributeUsage::Scalar,
        1);
    auto is_feature = attribute_vector_ref<int8_t>(sweep_arrangement, "is_feature");
    is_feature.setZero();
    Index num_edges = sweep_arrangement.get_num_edges();
    for (Index eid = 0; eid < num_edges; eid++) {
        Index edge_valence = sweep_arrangement.count_num_corners_around_edge(eid);
        if (edge_valence != 2) {
            bool feature_edge = false;
            sweep_arrangement.foreach_facet_around_edge(eid, [&](Index fid) {
                if (is_valid[fid] != 0) {
                    feature_edge = true;
                }
            });
            if (feature_edge) is_feature[eid] = 1;
        }
    }

    return sweep_arrangement;
}

template <typename Scalar, typename Index>
lagrange::SurfaceMesh<Scalar, Index> extract_sweep_surface_from_arrangement(
    lagrange::SurfaceMesh<Scalar, Index>& sweep_arrangement)
{
    Index num_arrangement_facets = sweep_arrangement.get_num_facets();
    auto V = vertex_view(sweep_arrangement);
    auto F = facet_view(sweep_arrangement);
    auto is_valid = attribute_vector_view<int8_t>(sweep_arrangement, "valid");
    auto time_values = attribute_vector_view<Scalar>(sweep_arrangement, "time");

    lagrange::SurfaceMesh<Scalar, Index> sweep_surface;
    sweep_surface.add_vertices(
        static_cast<Index>(V.rows()),
        {V.data(), static_cast<size_t>(V.size())});

    Index num_valid_facets = 0;
    for (Index fid = 0; fid < num_arrangement_facets; fid++) {
        if (is_valid[fid] != 0) num_valid_facets++;
    }
    sweep_surface.add_triangles(num_valid_facets);
    auto sweep_F = facet_ref(sweep_surface);

    Index count = 0;
    for (Index fid = 0; fid < num_arrangement_facets; fid++) {
        if (is_valid[fid] == 0) {
            continue;
        } else if (is_valid[fid] == 1) {
            sweep_F.row(count) = F.row(fid);
        } else {
            sweep_F.row(count) = F.row(fid).reverse();
        }
        count++;
    }

    sweep_surface.template create_attribute<Scalar>(
        "time",
        lagrange::AttributeElement::Corner,
        lagrange::AttributeUsage::Scalar,
        1);
    auto sweep_time_values = attribute_vector_ref<Scalar>(sweep_surface, "time");

    count = 0;
    for (Index fid = 0; fid < num_arrangement_facets; fid++) {
        if (is_valid[fid] == 0) {
            continue;
        }
        const Index sweep_c_begin = sweep_surface.get_facet_corner_begin(count);
        const Index arrang_c_begin = sweep_arrangement.get_facet_corner_begin(fid);
        if (is_valid[fid] == 1) {
            sweep_time_values[sweep_c_begin + 0] = time_values[arrang_c_begin + 0];
            sweep_time_values[sweep_c_begin + 1] = time_values[arrang_c_begin + 1];
            sweep_time_values[sweep_c_begin + 2] = time_values[arrang_c_begin + 2];
        } else {
            sweep_time_values[sweep_c_begin + 0] = time_values[arrang_c_begin + 2];
            sweep_time_values[sweep_c_begin + 1] = time_values[arrang_c_begin + 1];
            sweep_time_values[sweep_c_begin + 2] = time_values[arrang_c_begin + 0];
        }
        count++;
    }

    lagrange::remove_isolated_vertices(sweep_surface);
    return sweep_surface;
}

#endif /* post_processing_h */
