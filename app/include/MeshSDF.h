#pragma once

#include <stf/primitives/implicit_function.h>
#include <igl/AABB.h>
#include <igl/fast_winding_number.h>
#include <igl/read_triangle_mesh.h>

#include <array>
#include <iostream>
#include <stdexcept>

class MeshSDF : public stf::ImplicitFunction<3>
{
public:
    using Scalar  = stf::Scalar;

    MeshSDF(const std::string& filename, std::array<Scalar, 3> center, Scalar radius, Scalar offset=0)
        : m_offset(offset)
    {
        igl::read_triangle_mesh(filename,m_V,m_F);

        if (m_V.rows() == 0 || m_F.rows() == 0)
        {
            throw std::runtime_error("Failed to read mesh from file: " + filename);
        }

        // Normalize the mesh to fit within the sphere of given radius and center.
        Eigen::RowVector3d c = Eigen::Map<Eigen::RowVector3d>(center.data());
        Eigen::RowVector3d bbox_min = m_V.colwise().minCoeff();
        Eigen::RowVector3d bbox_max = m_V.colwise().maxCoeff();
        Eigen::RowVector3d bbox_center = (bbox_min + bbox_max) / 2.0;
        Scalar bbox_diag = (bbox_max - bbox_center).norm();
        m_V = ((m_V.rowwise() - bbox_center) * (radius / bbox_diag)).rowwise() + c;

        m_tree.init(m_V,m_F);
        int order = 2;
        igl::fast_winding_number(m_V, m_F, order, m_fwn_bvh);\
    }

    Scalar value(std::array<Scalar, 3> pos) const override
    {
        Eigen::RowVector3d p(pos[0], pos[1], pos[2]);
        int i;
        Eigen::RowVector3d c;
        auto sq_dist = m_tree.squared_distance(m_V, m_F, p, i, c);
        Scalar dist = std::sqrt(sq_dist);
        auto w = igl::fast_winding_number(m_fwn_bvh,2.0,p);
        return dist * (1 - 2 * w) + m_offset;
    }

    std::array<Scalar, 3> gradient(std::array<Scalar, 3> pos) const override
    {
        Eigen::RowVector3d p(pos[0], pos[1], pos[2]);
        int i;
        Eigen::RowVector3d c;
        auto sq_dist = m_tree.squared_distance(m_V, m_F, p, i, c);
        auto w = igl::fast_winding_number(m_fwn_bvh,2.0,p);

        Eigen::RowVector3d v = p - c;
        v.stableNormalize();
        v = v * (1 - 2 * w);
        return {v(0), v(1), v(2)};
    }


private:
    Eigen::MatrixXd m_V;
    Eigen::MatrixXi m_F;
    igl::AABB<Eigen::MatrixXd, 3> m_tree;
    igl::FastWindingNumberBVH m_fwn_bvh;
    Scalar m_offset;
};
