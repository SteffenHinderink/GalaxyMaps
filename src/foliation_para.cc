#include "foliation_para.h"

FoliationParameterization::FoliationParameterization(Mesh& mesh) : mesh_(mesh),
                                                                   surface_parameterized_(false),
                                                                   shelled_(false),
                                                                   d_(mesh.request_cell_property<std::vector<int>>("_")),
                                                                   h_(mesh.request_cell_property<mpq_class>("_")),
                                                                   surface_parameter_(mesh.request_vertex_property<Vector3q>("_")) {}

bool FoliationParameterization::set_surface_parameterization(std::string property) {
    if (!mesh_.vertex_property_exists<Vector3q>(property)) {
        surface_parameterized_ = false;
        return surface_parameterized_;
    }
    surface_parameter_ = mesh_.request_vertex_property<Vector3q>(property);

    for (auto h : mesh_.halffaces()) {
        if (mesh_.is_boundary(h)) {
            std::vector<Vector3q> triangle;
            for (auto hv : mesh_.halfface_vertices(h)) {
                triangle.push_back(surface_parameter_[hv]);
            }
            Vector3q n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
            if (n.dot(triangle[0]) <= 0) {
                surface_parameterized_ = false;
                return surface_parameterized_;
            }
        }
    }
    surface_parameterized_ = true;
    return surface_parameterized_;
}

bool FoliationParameterization::is_center(OpenVolumeMesh::CellHandle c) {
    std::vector<int> d_c = d_[c];
    return d_c[0] == 0 && d_c[1] == 0 && d_c[2] == 0 && d_c[3] == 0;
}

int& FoliationParameterization::direction_coordinate(OpenVolumeMesh::CellHandle c, OpenVolumeMesh::VertexHandle v) {
    int i = 0;
    for (auto cv : mesh_.cell_vertices(c)) {
        if (cv == v) {
            return d_[c][i];
        }
        i++;
    }
    std::cout << "Error: Vertex not part of tet" << std::endl;
    return d_[c][4];
}

mpq_class FoliationParameterization::direction_length(OpenVolumeMesh::CellHandle c) {
    if (is_center(c)) {
        std::vector<Eigen::Vector3d> vertices;
        for (auto cv : mesh_.tet_vertices(c)) {
            vertices.push_back(mesh_.vertex(cv));
        }
        Eigen::Matrix3d m;
        for (int i = 0; i < 3; i++) {
            m.col(i) = vertices[i + 1] - vertices[0];
        }
        double volume = std::abs(m.determinant() / 6);
        std::vector<double> face_areas;
        for (auto cf : mesh_.cell_faces(c)) {
            std::vector<Eigen::Vector3d> triangle;
            for (auto cfv : mesh_.face_vertices(cf)) {
                triangle.push_back(mesh_.vertex(cfv));
            }
            face_areas.push_back((triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]).norm() / 2);
        }
        double face_area_sum = 0;
        for (int i = 0; i < 4; i++) {
            face_area_sum += face_areas[i];
        }
        mpq_class h = 3 * volume / face_area_sum;
        if (h <= 0) {
            std::cout << "Error: Nonpositive inradius" << std::endl;
        }
        return h;
    } else {
        std::vector<mpq_class> a = {0.25, 0.25, 0.25, 0.25};
        std::vector<mpq_class> b(4);
        std::vector<int> d = d_[c];
        for (int i = 0; i < 4; i++) {
            b[i] = a[i] + d[i];
        }
        Eigen::Vector3d ac = Eigen::Vector3d::Zero();
        Eigen::Vector3d bc = Eigen::Vector3d::Zero();
        int i = 0;
        for (auto cv : mesh_.cell_vertices(c)) {
            Eigen::Vector3d pos = mesh_.vertex(cv);
            ac += a[i].get_d() * pos;
            bc += b[i].get_d() * pos;
            i++;
        }
        mpq_class h = (bc - ac).norm();
        if (h <= 0) {
            std::cout << "Error: Nonpositive direction length" << std::endl;
        }
        return h;
    }
}

bool FoliationParameterization::generate_directions() {
    d_ = mesh_.request_cell_property<std::vector<int>>("d");
    h_ = mesh_.request_cell_property<mpq_class>("h");

    auto marked = mesh_.request_face_property<bool>("marked");

    // c must have at least one marked face
    auto free = [this, marked](OpenVolumeMesh::CellHandle c) {
        int n_marked_faces = 0;
        for (auto cf : mesh_.cell_faces(c)) {
            if (marked[cf]) {
                n_marked_faces++;
            }
        }
        if (n_marked_faces > 2) {
            return true;
        }
        // 1 <= n_marked_faces <= 2

        int n_marked_edges = 0;
        for (auto ce : mesh_.cell_edges(c)) {
            bool edge_marked = false;
            for (auto cef : mesh_.edge_faces(ce)) {
                if (marked[cef]) {
                    edge_marked = true;
                }
            }
            if (edge_marked) {
                n_marked_edges++;
            }
        }
        if (n_marked_faces == 2) {
            return n_marked_edges == 5;
        }
        if (n_marked_edges > 3) {
            return false;
        }
        // n_marked_faces = 1, n_marked_edges = 3

        int n_marked_vertices = 0;
        for (auto cv : mesh_.cell_vertices(c)) {
            bool vertex_marked = false;
            for (auto cvf : mesh_.vertex_faces(cv)) {
                if (marked[cvf]) {
                    vertex_marked = true;
                }
            }
            if (vertex_marked) {
                n_marked_vertices++;
            }
        }
        return n_marked_vertices == 3;
    };

    enum State {
        BLOCKED,
        FREE,
        SHELLED
    };

    // The source is the entire boundary
    for (auto f : mesh_.faces()) {
        if (mesh_.is_boundary(f)) {
            marked[f] = true;
        } else {
            marked[f] = false;
        }
    }

    // Shelling order heuristic
    auto dist = mesh_.request_cell_property<int>("dist");
    auto visited = mesh_.request_cell_property<bool>("visited");
    auto cmp_greater = [dist](OpenVolumeMesh::CellHandle a, OpenVolumeMesh::CellHandle b) {
        return dist[a] > dist[b];
    };
    std::priority_queue<OpenVolumeMesh::CellHandle, std::vector<OpenVolumeMesh::CellHandle>, decltype(cmp_greater)> q_dijkstra(cmp_greater);
    // In
    for (auto c : mesh_.cells()) {
        visited[c] = false;
        if (mesh_.is_boundary(c)) {
            dist[c] = 0;
        } else {
            dist[c] = mesh_.n_cells();
        }
        q_dijkstra.push(c);
    }
    while (!q_dijkstra.empty()) {
        OpenVolumeMesh::CellHandle c = q_dijkstra.top();
        q_dijkstra.pop();
        if (!visited[c]) {
            visited[c] = true;
            for (auto cc : mesh_.cell_cells(c)) {
                if (dist[c] + 1 < dist[cc]) {
                    dist[cc] = dist[c] + 1;
                    q_dijkstra.push(cc);
                }
            }
        }
    }
    OpenVolumeMesh::CellHandle source = *mesh_.cells_begin();
    for (auto c : mesh_.cells()) {
        if (dist[c] > dist[source]) {
            source = c;
        }
    }
    // Out
    for (auto c : mesh_.cells()) {
        visited[c] = false;
        dist[c] = c == source ? 0 : mesh_.n_cells();
        q_dijkstra.push(c);
    }
    while (!q_dijkstra.empty()) {
        OpenVolumeMesh::CellHandle c = q_dijkstra.top();
        q_dijkstra.pop();
        if (!visited[c]) {
            visited[c] = true;
            for (auto cc : mesh_.cell_cells(c)) {
                if (dist[c] + 1 < dist[cc]) {
                    dist[cc] = dist[c] + 1;
                    q_dijkstra.push(cc);
                }
            }
        }
    }

    // Initialize cells
    auto cmp_less = [dist](OpenVolumeMesh::CellHandle a, OpenVolumeMesh::CellHandle b) {
        return dist[a] < dist[b];
    };
    std::priority_queue<OpenVolumeMesh::CellHandle, std::vector<OpenVolumeMesh::CellHandle>, decltype(cmp_less)> q(cmp_less);
    auto state = mesh_.request_cell_property<int>("state");
    for (auto c : mesh_.cells()) {
        if (mesh_.is_boundary(c) && free(c)) {
            state[c] = FREE;
            q.push(c);
        } else {
            state[c] = BLOCKED;
        }
    }

    int i = 0;
    while (!q.empty()) {
        OpenVolumeMesh::CellHandle c = q.top();
        q.pop();
        if (state[c] == FREE) {
            if (i == mesh_.n_cells() - 1) {
                // Sink
                d_[c] = {0, 0, 0, 0};
                state[c] = SHELLED;
            } else {
                // d_c
                int n_cell = 0;
                std::vector<OpenVolumeMesh::FaceHandle> incident;
                for (auto cf : mesh_.cell_faces(c)) {
                    incident.push_back(cf);
                    if (marked[cf]) {
                        n_cell++;
                    }
                }
                std::vector<int> d_c(4);
                int j = 0;
                for (auto cv : mesh_.cell_vertices(c)) {
                    int n_vertex = 0;
                    for (auto cvf : mesh_.vertex_faces(cv)) {
                        if (std::find(incident.begin(), incident.end(), cvf) != incident.end() && marked[cvf]) {
                            n_vertex++;
                        }
                    }
                    switch (n_vertex) {
                    case 0:
                        d_c[j] = 3;
                        break;
                    case 1:
                        d_c[j] = n_cell == 2 ? 1 : -1;
                        break;
                    case 2:
                        d_c[j] = n_cell == 2 ? -1 : 1;
                        break;
                    case 3:
                        d_c[j] = -3;
                        break;
                    }
                    j++;
                }
                d_[c] = d_c;

                // Remove c and mark revealed faces
                state[c] = SHELLED;
                for (auto cf : mesh_.cell_faces(c)) {
                    marked[cf] = true;
                }

                // Check if neighboring cells are still free
                for (auto cv : mesh_.cell_vertices(c)) {
                    for (auto cvc : mesh_.vertex_cells(cv)) {
                        if (state[cvc] == FREE && !free(cvc)) {
                            state[cvc] = BLOCKED;
                        }
                    }
                }

                // New free cells
                for (auto cc : mesh_.cell_cells(c)) {
                    if (state[cc] == BLOCKED && free(cc)) {
                        state[cc] = FREE;
                        q.push(cc);
                    }
                }
            }
            h_[c] = direction_length(c);
            i++;
        }
    }

    shelled_ = i == mesh_.n_cells();
    if (shelled_) {
        mesh_.set_persistent(d_);
        align_directions();
    }
    return shelled_;
}

void FoliationParameterization::align_directions() {
    // Align to edges
    for (auto e : mesh_.halfedges()) {
        if (mesh_.is_boundary(e)) {
            continue;
        }

        OpenVolumeMesh::VertexHandle u = mesh_.from_vertex_handle(e);
        OpenVolumeMesh::VertexHandle v = mesh_.to_vertex_handle(e);
        bool alignable = true;
        for (auto ec : mesh_.halfedge_cells(e)) {
            if (direction_coordinate(ec, u) >= 0) {
                alignable = false;
                break;
            }
            if (direction_coordinate(ec, v) <= 0) {
                alignable = false;
                break;
            }
        }
        if (!alignable) {
            continue;
        }

        for (auto ec : mesh_.halfedge_cells(e)) {
            direction_coordinate(ec, u) = -1;
            direction_coordinate(ec, v) = 1;
            OpenVolumeMesh::HalfEdgeHandle e2 = mesh_.halfedge_opposite_halfedge(e, ec);
            direction_coordinate(ec, mesh_.from_vertex_handle(e2)) = 0;
            direction_coordinate(ec, mesh_.to_vertex_handle(e2)) = 0;
        }

        mpq_class h = direction_length(*mesh_.hec_iter(e));
        for (auto ec : mesh_.halfedge_cells(e)) {
            h_[ec] = h;
        }
    }

    // Align to faces
    for (auto f : mesh_.faces()) {
        if (mesh_.is_boundary(f)) {
            continue;
        }

        std::vector<OpenVolumeMesh::VertexHandle> face_vertices;
        for (auto fv : mesh_.face_vertices(f)) {
            face_vertices.push_back(fv);
        }
        std::vector<OpenVolumeMesh::HalfFaceHandle> face_halffaces;
        std::vector<OpenVolumeMesh::CellHandle> face_cells;
        for (auto fh : mesh_.face_halffaces(f)) {
            face_halffaces.push_back(fh);
            face_cells.push_back(mesh_.incident_cell(fh));
        }
        bool alignable = true;
        std::vector<bool> positive(3);
        int n_positive = 0;
        for (int i = 0; i < 3; i++) {
            int dc = direction_coordinate(face_cells[0], face_vertices[i]);
            if (dc == 0) {
                alignable = false;
                break;
            }
            positive[i] = dc > 0;
            if (positive[i]) {
                n_positive++;
            }
        }
        if (n_positive == 0 || n_positive == 3) {
            alignable = false;
        }
        if (!alignable) {
            continue;
        }
        for (int i = 0; i < 3; i++) {
            int dc = direction_coordinate(face_cells[1], face_vertices[i]);
            if (dc == 0) {
                alignable = false;
                break;
            }
            if (positive[i] != (dc > 0)) {
                alignable = false;
                break;
            }
        }
        if (!alignable) {
            continue;
        }

        for (int i = 0; i < 2; i++) {
            direction_coordinate(face_cells[i], mesh_.halfface_opposite_vertex(face_halffaces[i])) = 0;
            for (int j = 0; j < 3; j++) {
                if (positive[j]) {
                    direction_coordinate(face_cells[i], face_vertices[j]) = n_positive == 2 ? 1 : 2;
                } else {
                    direction_coordinate(face_cells[i], face_vertices[j]) = n_positive == 2 ? -2 : -1;
                }
            }
        }

        mpq_class h = direction_length(face_cells[0]);
        for (auto fc : face_cells) {
            h_[fc] = h;
        }
    }
}

std::pair<std::vector<mpq_class>, mpq_class> FoliationParameterization::project(MeshPoint p, bool in) {
    std::vector<int> d_c = d_[p.c];
    if (d_c[0] == 0 && d_c[1] == 0 && d_c[2] == 0 && d_c[3] == 0) {
        // Projection in sink tetrahedron
        mpq_class min(1, 4);
        for (int i = 0; i < 4; i++) {
            if (p.bary[i] < min) {
                min = p.bary[i];
            }
        }
        if (min == mpq_class(1, 4)) {
            if (in) {
                return {{mpq_class(1, 4), mpq_class(1, 4), mpq_class(1, 4), mpq_class(1, 4)}, 0};
            }
            return {{1, 0, 0, 0}, h_[p.c]};
        }
        mpq_class lambda = min / (mpq_class(1, 4) - min);
        if (in) {
            return {{mpq_class(1, 4), mpq_class(1, 4), mpq_class(1, 4), mpq_class(1, 4)}, 1 / (1 + lambda) * h_[p.c]};
        }
        std::vector<mpq_class> l(4);
        for (int j = 0; j < 4; j++) {
            l[j] = p.bary[j] - lambda * (mpq_class(1, 4) - p.bary[j]);
        }
        return {l, lambda / (1 + lambda) * h_[p.c]};
    }
    std::vector<mpq_class> projection = p.bary;
    mpq_class length = 0;
    for (int i = 0; i < 4; i++) {
        if (p.bary[i] > 0 && d_c[i] != 0) {
            std::vector<mpq_class> l(4);
            mpq_class lambda = p.bary[i] / d_c[i];
            bool inside = true;
            for (int j = 0; j < 4 && inside; j++) {
                l[j] = p.bary[j] - lambda * d_c[j];
                if (l[j] < 0) {
                    inside = false;
                }
            }
            if (inside && ((in && lambda < length) || (!in && lambda > length))) {
                projection = l;
                length = lambda;
            }
        }
    }
    if (in) {
        length *= -1;
    }
    return {projection, length * h_[p.c]};
}

MeshPoint FoliationParameterization::next(MeshPoint p, bool in) {
    MeshPoint next;
    auto idx = mesh_.request_vertex_property<int>("idx");
    for (auto v : mesh_.vertices()) {
        idx[v] = -1;
    }
    int i = 0;
    for (auto cv : mesh_.cell_vertices(p.c)) {
        idx[cv] = i;
        i++;
    }
    int n_zeros = 0;
    for (i = 0; i < 4; i++) {
        if (p.bary[i] == 0) {
            n_zeros++;
        }
    }
    switch (n_zeros) {
    case 1:
        // Through face
        for (auto cf : mesh_.cell_faces(p.c)) {
            bool out = true;
            for (auto cfv : mesh_.face_vertices(cf)) {
                if (p.bary[idx[cfv]] == 0) {
                    out = false;
                    break;
                }
            }
            if (out) {
                for (auto cfc : mesh_.face_cells(cf)) {
                    if (cfc.is_valid()) {
                        if (is_center(cfc)) {
                            if (in) {
                                next.c = cfc;
                                break;
                            } else {
                                continue;
                            }
                        }
                        bool outgoing = true;
                        std::vector<int> d_cfc = d_[cfc];
                        i = 0;
                        for (auto cfcv : mesh_.cell_vertices(cfc)) {
                            bool part_of_face = false;
                            for (auto cfv : mesh_.face_vertices(cf)) {
                                if (cfcv == cfv) {
                                    part_of_face = true;
                                    break;
                                }
                            }
                            if (!part_of_face && ((in && d_cfc[i] < 0) || (!in && d_cfc[i] > 0))) {
                                outgoing = false;
                                break;
                            }
                            i++;
                        }
                        if (outgoing) {
                            next.c = cfc;
                            break;
                        }
                    }
                }
                break;
            }
        }
        break;
    case 2:
        // Through edge
        for (auto ce : mesh_.cell_edges(p.c)) {
            bool out = true;
            for (auto cev : mesh_.edge_vertices(ce)) {
                if (p.bary[idx[cev]] == 0) {
                    out = false;
                    break;
                }
            }
            if (out) {
                for (auto cec : mesh_.edge_cells(ce)) {
                    if (is_center(cec)) {
                        if (in) {
                            next.c = cec;
                            break;
                        } else {
                            continue;
                        }
                    }
                    bool outgoing = true;
                    std::vector<int> d_cec = d_[cec];
                    i = 0;
                    for (auto cecv : mesh_.cell_vertices(cec)) {
                        bool part_of_edge = false;
                        for (auto cev : mesh_.edge_vertices(ce)) {
                            if (cecv == cev) {
                                part_of_edge = true;
                                break;
                            }
                        }
                        if (!part_of_edge && ((in && d_cec[i] < 0) || (!in && d_cec[i] > 0))) {
                            outgoing = false;
                            break;
                        }
                        i++;
                    }
                    if (outgoing) {
                        next.c = cec;
                        break;
                    }
                }
                break;
            }
        }
        break;
    case 3:
        // Through vertex
        for (auto cv : mesh_.cell_vertices(p.c)) {
            if (p.bary[idx[cv]] != 0) {
                for (auto cvc : mesh_.vertex_cells(cv)) {
                    if (is_center(cvc)) {
                        if (in) {
                            next.c = cvc;
                            break;
                        } else {
                            continue;
                        }
                    }
                    bool outgoing = true;
                    std::vector<int> d_cvc = d_[cvc];
                    i = 0;
                    for (auto cvcv : mesh_.cell_vertices(cvc)) {
                        if (cvcv != cv && ((in && d_cvc[i] < 0) || (!in && d_cvc[i] > 0))) {
                            outgoing = false;
                            break;
                        }
                        i++;
                    }
                    if (outgoing) {
                        next.c = cvc;
                        break;
                    }
                }
                break;
            }
        }
        break;
    }
    if (next.c.is_valid()) {
        for (auto cv : mesh_.cell_vertices(next.c)) {
            if (idx[cv] < 0) {
                next.bary.push_back(0);
            } else {
                next.bary.push_back(p.bary[idx[cv]]);
            }
        }
    }
    return next;
}

Vector3q FoliationParameterization::psi(MeshPoint p) {
    if (!surface_parameterized_ || !shelled_) {
        throw FoliationParameterizationException();
    }

    // Trace leaf
    // In
    mpq_class length, length_in = 0;
    MeshPoint tmp = p;
    while (tmp.c.is_valid()) {
        std::tie(tmp.bary, length) = project(tmp, true);
        length_in += length;
        tmp = next(tmp, true);
    }
    // Out
    mpq_class length_sum = length_in;
    tmp = p;
    while (tmp.c.is_valid()) {
        p = tmp;
        std::tie(p.bary, length) = project(p, false);
        length_sum += length;
        tmp = next(p, false);
    }

    // Get position in star
    Vector3q q = Vector3q::Zero();
    int i = 0;
    for (auto cv : mesh_.cell_vertices(p.c)) {
        q += surface_parameter_[cv] * p.bary[i];
        i++;
    }
    q *= length_in / length_sum;

    return q;
}

MeshPoint FoliationParameterization::psi_inverse(Vector3q q) {
    if (!surface_parameterized_ || !shelled_) {
        throw FoliationParameterizationException();
    }

    if (q == Vector3q::Zero()) {
        MeshPoint p;
        for (auto c : mesh_.cells()) {
            std::vector<int> d_c = d_[c];
            if (d_c[0] == 0 && d_c[1] == 0 && d_c[2] == 0 && d_c[3] == 0) {
                p.c = c;
                break;
            }
        }
        p.bary = {mpq_class(1, 4), mpq_class(1, 4), mpq_class(1, 4), mpq_class(1, 4)};
        return p;
    }

    MeshPoint p;
    mpq_class lambda;
    for (auto f : mesh_.faces()) {
        if (mesh_.is_boundary(f)) {
            std::vector<Vector3q> triangle;
            for (auto fv : mesh_.face_vertices(f)) {
                triangle.push_back(surface_parameter_[fv]);
            }
            Vector3q n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
            if (n.dot(q) == 0) {
                // Parallel
                continue;
            }
            lambda = n.dot(triangle[0]) / n.dot(q);
            if (lambda < 1) {
                // q behind s
                continue;
            }
            Vector3q s = q * lambda;

            auto area = [n](Vector3q a, Vector3q b, Vector3q c) {
                Matrix3q m;
                m << b - a, c - a, n;
                return m.determinant();
            };

            mpq_class a = area(triangle[0], triangle[1], triangle[2]);
            std::vector<mpq_class> triangle_bary(3);
            bool inside = true;
            for (int i = 0; i < 3 && inside; i++) {
                triangle_bary[i] = area(s, triangle[(i + 1) % 3], triangle[(i + 2) % 3]) / a;
                if (triangle_bary[i] < 0 || triangle_bary[i] > 1) {
                    inside = false;
                }
            }
            if (inside) {
                auto idx = mesh_.request_vertex_property<int>("idx");
                for (auto v : mesh_.vertices()) {
                    idx[v] = -1;
                }
                int i = 0;
                for (auto fv : mesh_.face_vertices(f)) {
                    idx[fv] = i;
                    i++;
                }
                for (auto fc : mesh_.face_cells(f)) {
                    if (fc.is_valid()) {
                        p.c = fc;
                        for (auto fcv : mesh_.cell_vertices(fc)) {
                            p.bary.push_back(idx[fcv] == -1 ? 0 : triangle_bary[idx[fcv]]);
                        }
                        break;
                    }
                }
                break;
            }
        }
    }
    if (!p.c.is_valid()) {
        return p;
    }

    // Trace leaf
    mpq_class length, length_sum = 0;
    MeshPoint tmp = p;
    while (tmp.c.is_valid()) {
        std::tie(tmp.bary, length) = project(tmp);
        length_sum += length;
        tmp = next(tmp);
    }

    // Trace leaf up to depth
    mpq_class distance = length_sum * (1 - 1 / lambda);
    while (p.c.is_valid()) {
        tmp = p;
        std::tie(tmp.bary, length) = project(p);
        if (length > 0 && length >= distance) {
            mpq_class alpha = distance / length;
            for (int i = 0; i < 4; i++) {
                p.bary[i] = alpha * tmp.bary[i] + (1 - alpha) * p.bary[i];
            }
            break;
        }
        distance -= length;
        p = next(tmp);
    }

    return p;
}

Mesh FoliationParameterization::refine(bool verbose) {
    if (!surface_parameterized_ || !shelled_) {
        throw FoliationParameterizationException();
    }
    if (verbose) {
        std::cout << "Refining mesh ..." << std::endl;
    }

    auto pattern = mesh_.request_face_property<CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>>>("pattern");
    auto idx = mesh_.request_vertex_property<int>("idx");
    auto aligned = mesh_.request_face_property<bool>("aligned");
    for (auto f : mesh_.faces()) {
        OpenVolumeMesh::HalfFaceHandle h = mesh_.face_halffaces(f)[0];
        if (mesh_.is_boundary(h)) {
            h = mesh_.face_halffaces(f)[1];
        }
        OpenVolumeMesh::CellHandle c = mesh_.incident_cell(h);
        if (is_center(c)) {
            continue;
        }
        OpenVolumeMesh::VertexHandle v = mesh_.halfface_opposite_vertex(h);
        if (direction_coordinate(c, v) == 0) {
            aligned[f] = true;
        }
    }

    auto slice = [this, &pattern, &idx, &aligned](std::vector<MeshPoint> edge) {
        for (bool in : {false, true}) {
            std::queue<std::vector<MeshPoint>> blades;
            blades.push(edge);
            while (!blades.empty()) {
                std::vector<MeshPoint> blade = blades.front();
                blades.pop();
                for (auto v : mesh_.vertices()) {
                    idx[v] = -1;
                }
                MeshPoint center = {blade[0].c, {0, 0, 0, 0}};
                int ii = 0;
                for (auto cv : mesh_.cell_vertices(blade[0].c)) {
                    center.bary[ii] = (blade[0].bary[ii] + blade[1].bary[ii]) / 2;
                    idx[cv] = ii;
                    ii++;
                }

                // Insert blade into all patterns that it touches
                for (auto f : mesh_.faces()) {
                    if (aligned[f]) {
                        continue;
                    }
                    std::vector<std::vector<mpq_class>> face_blade(2);
                    std::vector<mpq_class> sum = {0, 0};
                    for (auto fv : mesh_.face_vertices(f)) {
                        if (idx[fv] == -1) {
                            face_blade[0].push_back(0);
                            face_blade[1].push_back(0);
                        } else {
                            face_blade[0].push_back(blade[0].bary[idx[fv]]);
                            sum[0] += blade[0].bary[idx[fv]];
                            face_blade[1].push_back(blade[1].bary[idx[fv]]);
                            sum[1] += blade[1].bary[idx[fv]];
                        }
                    }
                    if (sum[0] == 1 || sum[1] == 1) {
                        std::vector<CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>::Point_2> v = {
                            {face_blade[0][0], face_blade[0][1]},
                            {face_blade[1][0], face_blade[1][1]}};
                        if (sum[0] == 1 && sum[1] == 1) {
                            std::vector<CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>::X_monotone_curve_2> e = {{v[0], v[1]}};
                            CGAL::insert(pattern[f], e.begin(), e.end());
                        } else if (sum[0] == 1) {
                            CGAL::insert_point(pattern[f], v[0]);
                        } else {
                            CGAL::insert_point(pattern[f], v[1]);
                        }
                    }
                }

                // Blade endpoints a and b in next tetrahedron c with direction d
                OpenVolumeMesh::CellHandle c = next(center, in).c;
                if (c.is_valid()) {
                    std::vector<int> d = d_[c];
                    if (!(d[0] == 0 && d[1] == 0 && d[2] == 0 && d[3] == 0)) {
                        std::vector<mpq_class> a(4);
                        std::vector<mpq_class> b(4);
                        ii = 0;
                        for (auto cv : mesh_.cell_vertices(c)) {
                            a[ii] = idx[cv] < 0 ? 0 : blade[0].bary[idx[cv]];
                            b[ii] = idx[cv] < 0 ? 0 : blade[1].bary[idx[cv]];
                            ii++;
                        }

                        // Project blade
                        std::vector<std::pair<std::vector<mpq_class>, mpq_class>> projections;
                        projections.push_back({project({c, a}, in).first, 1});
                        projections.push_back({project({c, b}, in).first, 0});
                        std::vector<int> uv_indices;
                        for (int i = 0; i < 4; i++) {
                            if ((in && d[i] < 0) || (!in && d[i] > 0)) {
                                uv_indices.push_back(i);
                            }
                        }
                        for (int i = 0; i < uv_indices.size(); i++) {
                            for (int j = i + 1; j < uv_indices.size(); j++) {
                                int u = uv_indices[i];
                                int v = uv_indices[j];
                                // Prevent division by 0 when blade lies on edge
                                if (b[u] - a[u] - d[u] / d[v] * (b[v] - a[v]) != 0) {
                                    mpq_class alpha = (b[u] - d[u] / d[v] * b[v]) / (b[u] - a[u] - d[u] / d[v] * (b[v] - a[v]));
                                    if (alpha > 0 && alpha < 1) {
                                        mpq_class lambda = (alpha * a[u] + (1 - alpha) * b[u]) / d[u];
                                        std::vector<mpq_class> projection(4);
                                        bool inside = true;
                                        for (int k = 0; k < 4 && inside; k++) {
                                            projection[k] = alpha * a[k] + (1 - alpha) * b[k] - lambda * d[k];
                                            if (projection[k] < 0) {
                                                inside = false;
                                            }
                                        }
                                        if (inside) {
                                            projections.push_back({projection, alpha});
                                        }
                                    }
                                }
                            }
                        }
                        std::sort(projections.begin(), projections.end(), [](auto a, auto b) {
                            return a.second < b.second;
                        });

                        // Enqueue projection of blade as new blades
                        mpq_class last_alpha = 0;
                        for (int i = 1; i < projections.size(); i++) {
                            // Blade must have a length > 0
                            if (projections[i].second > last_alpha) {
                                int n_common_zeros = 0;
                                for (int j = 0; j < 4; j++) {
                                    if (projections[i - 1].first[j] == 0 && projections[i].first[j] == 0) {
                                        n_common_zeros++;
                                    }
                                }
                                // Blade must not lie on an edge
                                if (n_common_zeros < 2) {
                                    blades.push({{c, projections[i - 1].first}, {c, projections[i].first}});
                                }
                                last_alpha = projections[i].second;
                            }
                        }
                    }
                }
            }
        }
    };

    if (verbose) {
        std::cout << "Creating cut patterns" << std::endl;
    }
    // Fill patterns by slicing through the mesh with every edge
    for (auto e : mesh_.edges()) {
        if (verbose) {
            std::cout << "Edge " << e.idx() + 1 << " / " << mesh_.n_edges() << std::endl;
        }
        // Get e as two mesh points so that it can be used as a blade for slicing
        OpenVolumeMesh::CellHandle ec = *mesh_.edge_cells(e).first;
        std::vector<MeshPoint> edge = {{ec, {0, 0, 0, 0}}, {ec, {0, 0, 0, 0}}};
        int i = 0;
        for (auto ecv : mesh_.cell_vertices(ec)) {
            idx[ecv] = i;
            i++;
        }
        i = 0;
        for (auto ev : mesh_.edge_vertices(e)) {
            edge[i].bary[idx[ev]] = 1;
            i++;
        }

        slice(edge);
    }

    if (verbose) {
        std::cout << "Triangulating cut patterns" << std::endl;
    }
    // Triangulate patterns
    for (auto f : mesh_.faces()) {
        if (verbose) {
            std::cout << "Face " << f.idx() + 1 << " / " << mesh_.n_faces() << std::endl;
        }
        if (aligned[f]) {
            continue;
        }
        CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>> arr = pattern[f];
        for (CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>>::Face_const_iterator f_it = arr.faces_begin(); f_it != arr.faces_end(); f_it++) {
            if (!f_it->is_unbounded()) {
                std::vector<CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>::Point_2> polygon;
                CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>>::Ccb_halfedge_const_circulator circ = f_it->outer_ccb();
                CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>>::Ccb_halfedge_const_circulator curr = circ;
                do {
                    polygon.push_back(curr->target()->point());
                } while (++curr != circ);
                if (polygon.size() > 3) {

                    auto linear_dependent = [](CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>::Vector_2 u, CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>::Vector_2 v) {
                        if (v.x() == 0) {
                            return u.x() == 0;
                        }
                        mpq_class lambda = u.x() / v.x();
                        return u.y() == lambda * v.y();
                    };

                    int offset = -1;
                    for (int i = 0; i < polygon.size() && offset == -1; i++) {
                        offset = i;
                        CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>::Point_2 a = polygon[i];
                        CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>::Point_2 b = polygon[(i + 1) % polygon.size()];
                        CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>::Point_2 c = polygon[(i + 2) % polygon.size()];
                        if (linear_dependent(b - a, c - a)) {
                            offset = -1;
                        }
                        b = polygon[(polygon.size() + i - 1) % polygon.size()];
                        c = polygon[(polygon.size() + i - 2) % polygon.size()];
                        if (linear_dependent(b - a, c - a)) {
                            offset = -1;
                        }
                    }

                    CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>::Point_2 a = polygon[offset];
                    std::vector<mpq_class> face_a(3);
                    face_a[0] = a.x();
                    face_a[1] = a.y();
                    face_a[2] = 1 - a.x() - a.y();
                    for (int i = 2; i < polygon.size() - 1; i++) {
                        CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>::Point_2 b = polygon[(i + offset) % polygon.size()];
                        std::vector<mpq_class> face_b(3);
                        face_b[0] = b.x();
                        face_b[1] = b.y();
                        face_b[2] = 1 - b.x() - b.y();
                        for (auto fc : mesh_.face_cells(f)) {
                            if (fc.is_valid()) {
                                std::vector<MeshPoint> edge = {{fc, {0, 0, 0, 0}}, {fc, {0, 0, 0, 0}}};
                                int j = 0;
                                for (auto fcv : mesh_.cell_vertices(fc)) {
                                    idx[fcv] = j;
                                    j++;
                                }
                                j = 0;
                                for (auto fv : mesh_.face_vertices(f)) {
                                    edge[0].bary[idx[fv]] = face_a[j];
                                    edge[1].bary[idx[fv]] = face_b[j];
                                    j++;
                                }
                                slice(edge);
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    if (verbose) {
        std::cout << "Cutting tetrahedra" << std::endl;
    }
    // Refined mesh
    Mesh refined;
    // Vertices of the refined mesh in sparse barycentric coordinates
    // Each vertex is sorted by the indices of its spanning vertices
    std::vector<std::vector<std::pair<int, mpq_class>>> vertices;
    // Corresponding vertex handles of the refined mesh
    std::vector<OpenVolumeMesh::VertexHandle> handles;
    // Map from vertices to their indices
    std::map<std::vector<std::pair<int, mpq_class>>, int> vertex_map;

    auto refined_vertex = refined.request_vertex_property<std::vector<std::pair<int, mpq_class>>>("vertex");
    refined.set_persistent(refined_vertex);
    auto refined_position = refined.request_vertex_property<Vector3q>("position");
    refined.set_persistent(refined_position);

    auto arrpoint_to_vertex = [this](OpenVolumeMesh::FaceHandle f, CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>::Point_2 p) {
        std::vector<mpq_class> face_bary(3);
        face_bary[0] = p.x();
        face_bary[1] = p.y();
        face_bary[2] = 1 - p.x() - p.y();
        std::vector<std::pair<int, mpq_class>> vertex;
        int i = 0;
        for (auto fv : mesh_.face_vertices(f)) {
            if (face_bary[i] > 0) {
                vertex.push_back({fv.idx(), face_bary[i]});
            }
            i++;
        }
        std::sort(vertex.begin(), vertex.end(), [](auto a, auto b) {
            return a.first < b.first;
        });
        return vertex;
    };

    auto meshpoint_to_vertex = [this](MeshPoint p) {
        std::vector<std::pair<int, mpq_class>> vertex;
        int i = 0;
        for (auto cv : mesh_.cell_vertices(p.c)) {
            if (p.bary[i] > 0) {
                vertex.push_back({cv.idx(), p.bary[i]});
            }
            i++;
        }
        std::sort(vertex.begin(), vertex.end(), [](auto a, auto b) {
            return a.first < b.first;
        });
        return vertex;
    };

    auto vertex_to_meshpoint = [this](OpenVolumeMesh::CellHandle c, std::vector<std::pair<int, mpq_class>> vertex) {
        MeshPoint p = {c, {0, 0, 0, 0}};
        int i = 0;
        for (auto cv : mesh_.cell_vertices(c)) {
            for (int j = 0; j < vertex.size(); j++) {
                if (vertex[j].first == cv.idx()) {
                    p.bary[i] = vertex[j].second;
                }
            }
            i++;
        }
        return p;
    };

    auto add_vertex = [this, &refined, &refined_position, &refined_vertex](std::vector<std::pair<int, mpq_class>> vertex) {
        Vector3q position = Vector3q::Zero();
        for (int i = 0; i < vertex.size(); i++) {
            Eigen::Vector3d pos = mesh_.vertex(*(mesh_.vertices_begin() + vertex[i].first));
            position += pos.cast<mpq_class>() * vertex[i].second;
        }
        OpenVolumeMesh::VertexHandle v = refined.add_vertex(Eigen::Vector3d(position[0].get_d(), position[1].get_d(), position[2].get_d()));
        refined_position[v] = position;
        refined_vertex[v] = vertex;
        return v;
    };

    auto add_tetrahedron = [&refined, &handles](std::vector<int> tet) {
        refined.add_cell(handles[tet[1]], handles[tet[0]], handles[tet[2]], handles[tet[3]], true);
    };

    // Create vertices on faces
    for (auto f : mesh_.faces()) {
        if (aligned[f]) {
            continue;
        }
        CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>> arr = pattern[f];
        for (CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>>::Vertex_const_iterator v_it = arr.vertices_begin(); v_it != arr.vertices_end(); v_it++) {
            std::vector<std::pair<int, mpq_class>> vertex = arrpoint_to_vertex(f, v_it->point());
            if (vertex_map.find(vertex) == vertex_map.end()) {
                vertices.push_back(vertex);
                handles.push_back(add_vertex(vertex));
                vertex_map.insert({vertex, vertices.size() - 1});
            }
        }
    }

    std::vector<OpenVolumeMesh::VertexHandle> spiral_prism_centers;

    // Cut every tetrahedron along the patterns on its faces
    for (auto cell : mesh_.cells()) {
        if (verbose) {
            std::cout << "Cell " << cell.idx() + 1 << " / " << mesh_.n_cells() << std::endl;
        }
        // Cut in the foliation direction
        // The incoming faces are opposite of the vertices with a positive coordinate
        std::vector<int> d_c = d_[cell];
        int sink = -1;
        if (d_c[0] == 0 && d_c[1] == 0 && d_c[2] == 0 && d_c[3] == 0) {
            std::vector<std::pair<int, mpq_class>> vertex;
            for (auto cv : mesh_.cell_vertices(cell)) {
                vertex.push_back({cv.idx(), mpq_class(1, 4)});
            }
            vertices.push_back(vertex);
            handles.push_back(add_vertex(vertex));
            sink = vertices.size() - 1;
        }
        int side = 0;
        for (auto cv : mesh_.cell_vertices(cell)) {
            if (d_c[side] > 0 || is_center(cell)) {
                for (auto cf : mesh_.cell_faces(cell)) {
                    bool opposite = true;
                    for (auto cfv : mesh_.face_vertices(cf)) {
                        if (cfv == cv) {
                            opposite = false;
                            break;
                        }
                    }
                    if (opposite) {

                        // Get triangles of the face in face orientation
                        std::vector<std::vector<int>> face_triangles;
                        CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>> arr = pattern[cf];
                        for (CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>>::Face_const_iterator f_it = arr.faces_begin(); f_it != arr.faces_end(); f_it++) {
                            if (!f_it->is_unbounded()) {
                                std::vector<int> face_triangle;
                                CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>>::Ccb_halfedge_const_circulator circ = f_it->outer_ccb();
                                CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<CGAL::Cartesian<mpq_class>>>::Ccb_halfedge_const_circulator curr = circ;
                                do {
                                    face_triangle.push_back(vertex_map[arrpoint_to_vertex(cf, curr->target()->point())]);
                                } while (++curr != circ);
                                face_triangles.push_back(face_triangle);
                            }
                        }

                        // Flip triangles if the face is flipped
                        if (mesh_.incident_cell(mesh_.halfface_handle(cf, 0)) == cell) {
                            for (int i = 0; i < face_triangles.size(); i++) {
                                int tmp = face_triangles[i][1];
                                face_triangles[i][1] = face_triangles[i][2];
                                face_triangles[i][2] = tmp;
                            }
                        }

                        // Cut from every face triangle
                        for (auto t : face_triangles) {
                            if (sink != -1) {
                                // Refinement in sink
                                add_tetrahedron({t[0], t[1], t[2], sink});
                            } else {
                                std::vector<int> proj(3);
                                std::vector<int> indices;
                                for (int i = 0; i < 3; i++) {
                                    proj[i] = vertex_map[meshpoint_to_vertex({cell, project(vertex_to_meshpoint(cell, vertices[t[i]])).first})];
                                    if (proj[i] != t[i]) {
                                        indices.push_back(i);
                                    }
                                }
                                switch (indices.size()) {
                                case 1:
                                    // Tetrahedron
                                    add_tetrahedron({t[0], t[1], t[2], proj[indices[0]]});
                                    break;
                                case 2: {
                                    // Pyramid made of two tetrahedra
                                    int a = 3 - indices[0] - indices[1];
                                    int b = (a + 1) % 3;
                                    int c = (b + 1) % 3;
                                    OpenVolumeMesh::VertexHandle vc = handles[t[c]];
                                    OpenVolumeMesh::VertexHandle vpb = handles[proj[b]];
                                    bool cw = false;
                                    for (auto vv : refined.vertex_vertices(vc)) {
                                        if (vv == vpb) {
                                            cw = true;
                                            break;
                                        }
                                    }
                                    if (cw) {
                                        add_tetrahedron({t[a], t[b], t[c], proj[b]});
                                        add_tetrahedron({proj[c], proj[b], t[a], t[c]});
                                    } else {
                                        add_tetrahedron({t[a], t[b], t[c], proj[c]});
                                        add_tetrahedron({proj[c], proj[b], t[a], t[b]});
                                    }
                                } break;
                                case 3:
                                    // Prism made of three tetrahedra
                                    // Get directions of the diagonals of the quad sides of the prism
                                    std::vector<int> diag(3, 0);
                                    int sum = 0;
                                    for (int i = 0; i < 3; i++) {
                                        OpenVolumeMesh::VertexHandle viplus = handles[t[(i + 1) % 3]];
                                        OpenVolumeMesh::VertexHandle vpi = handles[proj[i]];
                                        bool cw = false;
                                        for (auto vv : refined.vertex_vertices(viplus)) {
                                            if (vv == vpi) {
                                                cw = true;
                                                break;
                                            }
                                        }
                                        if (cw) {
                                            diag[i] = -1;
                                            sum--;
                                        } else {
                                            OpenVolumeMesh::VertexHandle vi = handles[t[i]];
                                            OpenVolumeMesh::VertexHandle vpiplus = handles[proj[(i + 1) % 3]];
                                            bool ccw = false;
                                            for (auto vv : refined.vertex_vertices(vi)) {
                                                if (vv == vpiplus) {
                                                    ccw = true;
                                                    break;
                                                }
                                            }
                                            if (ccw) {
                                                diag[i] = 1;
                                                sum++;
                                            }
                                        }
                                    }
                                    if (std::abs(sum) == 3) {
                                        // Spiral prism
                                        MeshPoint center = {cell, {0, 0, 0, 0}};
                                        for (int i = 0; i < 3; i++) {
                                            MeshPoint p = vertex_to_meshpoint(cell, vertices[t[i]]);
                                            MeshPoint q = vertex_to_meshpoint(cell, vertices[proj[i]]);
                                            for (int j = 0; j < 4; j++) {
                                                center.bary[j] += p.bary[j] + q.bary[j];
                                            }
                                        }
                                        for (int i = 0; i < 4; i++) {
                                            center.bary[i] /= 6;
                                        }
                                        std::vector<std::pair<int, mpq_class>> vertex = meshpoint_to_vertex(center);
                                        vertices.push_back(vertex);
                                        OpenVolumeMesh::VertexHandle spiral_prism_center = add_vertex(vertex);
                                        handles.push_back(spiral_prism_center);
                                        spiral_prism_centers.push_back(spiral_prism_center);
                                        int x = vertices.size() - 1;
                                        add_tetrahedron({t[0], t[1], t[2], x});
                                        for (int i = 0; i < 3; i++) {
                                            if (diag[0] == 1) {
                                                add_tetrahedron({t[i], proj[(i + 1) % 3], t[(i + 1) % 3], x});
                                                add_tetrahedron({t[i], proj[i], proj[(i + 1) % 3], x});
                                            } else {
                                                add_tetrahedron({t[i], proj[i], t[(i + 1) % 3], x});
                                                add_tetrahedron({t[(i + 1) % 3], proj[i], proj[(i + 1) % 3], x});
                                            }
                                        }
                                        add_tetrahedron({proj[2], proj[1], proj[0], x});
                                    } else {
                                        int next = sum > 0 ? -1 : 1;
                                        for (int i = 0; i < 3; i++) {
                                            if (diag[i] == 0) {
                                                diag[i] = next;
                                                next *= -1;
                                            }
                                        }
                                        // Construct prism with tetrahedra
                                        for (int i = 0; i < 3; i++) {
                                            int a = i;
                                            int b = (i + 1) % 3;
                                            int c = (i + 2) % 3;
                                            if (diag[a] == 1 && diag[b] == -1) {
                                                // Front of the prism
                                                add_tetrahedron({t[a], t[b], t[c], proj[b]});
                                            } else if (diag[a] == -1 && diag[b] == 1) {
                                                // Back of the prism
                                                add_tetrahedron({proj[c], proj[b], proj[a], t[b]});
                                            } else if (diag[a] == 1 && diag[b] == 1) {
                                                // Middle of the prism with positive diagonal directions
                                                add_tetrahedron({t[a], t[b], proj[c], proj[b]});
                                            } else {
                                                // Middle of the prism with negative diagonal directions
                                                add_tetrahedron({t[c], proj[b], proj[a], t[b]});
                                            }
                                        }
                                    }
                                    break;
                                }
                            }
                        }

                        break;
                    }
                }
            }
            side++;
        }
    }

    if (verbose) {
        std::cout << "Parameterizing refined mesh" << std::endl;
    }
    // Parameterization of refined mesh
    auto refined_psi = refined.request_vertex_property<Vector3q>("qarameter");
    refined.set_persistent(refined_psi);
    for (int i = 0; i < handles.size(); i++) {
        for (auto c : mesh_.cells()) {
            int cnt = 0;
            for (auto cv : mesh_.cell_vertices(c)) {
                for (auto x : vertices[i]) {
                    if (x.first == cv.idx()) {
                        cnt++;
                    }
                }
            }
            if (cnt == vertices[i].size()) {
                refined_psi[handles[i]] = psi(vertex_to_meshpoint(c, vertices[i]));
                break;
            }
        }
    }
    for (auto v : spiral_prism_centers) {
        Vector3q center = Vector3q::Zero();
        for (auto vv : refined.vertex_vertices(v)) {
            center += refined_psi[vv];
        }
        center /= 6;
        refined_psi[v] = center;
    }

    if (verbose) {
        std::cout << "... Refinement finished" << std::endl;
    }

    return refined;
}
