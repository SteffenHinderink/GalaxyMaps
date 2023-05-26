#include "mesh.h"

OpenVolumeMesh::HalfEdgeHandle Mesh::halfedge_opposite_halfedge(OpenVolumeMesh::HalfEdgeHandle e, OpenVolumeMesh::CellHandle c) {
    OpenVolumeMesh::VertexHandle a = from_vertex_handle(e);
    OpenVolumeMesh::VertexHandle b = to_vertex_handle(e);
    std::vector<OpenVolumeMesh::VertexHandle> opposite_vertices;
    for (auto cv : tet_vertices(c)) {
        if (cv != a && cv != b) {
            opposite_vertices.push_back(cv);
        }
    }
    if (halfface_opposite_vertex(halfface({a, b, opposite_vertices[0]})) != opposite_vertices[1]) {
        return halfedge(opposite_vertices[1], opposite_vertices[0]);
    }
    return halfedge(opposite_vertices[0], opposite_vertices[1]);
}

OpenVolumeMesh::HalfFaceHandle Mesh::vertex_opposite_halfface(OpenVolumeMesh::VertexHandle v, OpenVolumeMesh::CellHandle c) {
    for (auto ch : cell_halffaces(c)) {
        bool has_v = false;
        for (auto chv : halfface_vertices(ch)) {
            if (chv == v) {
                has_v = true;
                break;
            }
        }
        if (!has_v) {
            return ch;
        }
    }
    return InvalidHalfFaceHandle;
}

void Mesh::remove_tet(OpenVolumeMesh::CellHandle c) {
    std::vector<OpenVolumeMesh::VertexHandle> vertices;
    for (auto cv : cell_vertices(c)) {
        vertices.push_back(cv);
    }
    std::vector<OpenVolumeMesh::EdgeHandle> edges;
    for (auto ce : cell_edges(c)) {
        edges.push_back(ce);
    }
    std::vector<OpenVolumeMesh::FaceHandle> faces;
    for (auto cf : cell_faces(c)) {
        faces.push_back(cf);
    }
    delete_cell(c);
    for (auto f : faces) {
        auto [c0, c1] = face_cells(f);
        if (!c0.is_valid() && !c1.is_valid()) {
            delete_face(f);
        }
    }
    for (auto e : edges) {
        if (!ec_iter(e).is_valid()) {
            delete_edge(e);
        }
    }
    for (auto v : vertices) {
        if (!vc_iter(v).is_valid()) {
            delete_vertex(v);
        }
    }
}

OpenVolumeMesh::VertexHandle Mesh::split_tet(OpenVolumeMesh::CellHandle c, Eigen::Vector3d p) {
    OpenVolumeMesh::VertexHandle v = add_vertex(p);
    std::vector<std::vector<OpenVolumeMesh::VertexHandle>> tets;
    for (auto ch : cell_halffaces(c)) {
        std::vector<OpenVolumeMesh::VertexHandle> tet;
        for (auto chv : halfface_vertices(ch)) {
            tet.push_back(chv);
        }
        tet.push_back(v);
        tets.push_back(tet);
    }
    delete_cell(c);
    for (auto tet : tets) {
        add_cell(tet, true);
    }
    return v;
}

std::vector<OpenVolumeMesh::VertexHandle> Mesh::boundary_vertices() {
    std::vector<OpenVolumeMesh::VertexHandle> boundary_vertices;
    for (auto v : vertices()) {
        if (is_boundary(v)) {
            boundary_vertices.push_back(v);
        }
    }
    return boundary_vertices;
}

std::vector<OpenVolumeMesh::EdgeHandle> Mesh::boundary_edges() {
    std::vector<OpenVolumeMesh::EdgeHandle> boundary_edges;
    for (auto e : edges()) {
        if (is_boundary(e)) {
            boundary_edges.push_back(e);
        }
    }
    return boundary_edges;
}

std::vector<OpenVolumeMesh::HalfFaceHandle> Mesh::boundary_halffaces() {
    std::vector<OpenVolumeMesh::HalfFaceHandle> boundary_halffaces;
    for (auto h : halffaces()) {
        if (is_boundary(h)) {
            boundary_halffaces.push_back(h);
        }
    }
    return boundary_halffaces;
}

OpenVolumeMesh::HalfFaceHandle Mesh::boundary_halfface(OpenVolumeMesh::HalfEdgeHandle e) {
    for (auto eh : halfedge_halffaces(e)) {
        if (is_boundary(eh)) {
            return eh;
        }
    }
    return InvalidHalfFaceHandle;
}

OpenVolumeMesh::HalfEdgeHandle Mesh::boundary_next(OpenVolumeMesh::HalfEdgeHandle e) {
    return next_halfedge_in_halfface(e, boundary_halfface(e));
}

OpenVolumeMesh::HalfEdgeHandle Mesh::boundary_out(OpenVolumeMesh::VertexHandle v) {
    for (auto ve : vertex_edges(v)) {
        if (is_boundary(ve)) {
            OpenVolumeMesh::HalfEdgeHandle e = halfedge_handle(ve, 0);
            if (from_vertex_handle(e) == v) {
                return e;
            } else {
                return opposite_halfedge_handle(e);
            }
        }
    }
    return InvalidHalfEdgeHandle;
}
