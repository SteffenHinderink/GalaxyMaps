#pragma once

#include <Eigen/Dense>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>

class Mesh : public OpenVolumeMesh::TetrahedralGeometryKernel<Eigen::Vector3d> {
public:
    OpenVolumeMesh::HalfEdgeHandle halfedge_opposite_halfedge(OpenVolumeMesh::HalfEdgeHandle e, OpenVolumeMesh::CellHandle c);

    OpenVolumeMesh::HalfFaceHandle vertex_opposite_halfface(OpenVolumeMesh::VertexHandle v, OpenVolumeMesh::CellHandle c);

    void remove_tet(OpenVolumeMesh::CellHandle c);

    OpenVolumeMesh::VertexHandle split_tet(OpenVolumeMesh::CellHandle c, Eigen::Vector3d p);

    std::vector<OpenVolumeMesh::VertexHandle> boundary_vertices();

    std::vector<OpenVolumeMesh::EdgeHandle> boundary_edges();

    std::vector<OpenVolumeMesh::HalfFaceHandle> boundary_halffaces();

    OpenVolumeMesh::HalfFaceHandle boundary_halfface(OpenVolumeMesh::HalfEdgeHandle e);

    OpenVolumeMesh::HalfEdgeHandle boundary_next(OpenVolumeMesh::HalfEdgeHandle e);

    OpenVolumeMesh::HalfEdgeHandle boundary_out(OpenVolumeMesh::VertexHandle v);
};
