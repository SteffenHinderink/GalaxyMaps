#include "galaxy.h"

bool inverted(Mesh& mesh, OpenVolumeMesh::CellHandle c, bool parameter_space) {
    OpenVolumeMesh::VertexPropertyT<Vector3q> qarameter = mesh.request_vertex_property<Vector3q>("qarameter");
    std::vector<Vector3q> vertices;
    for (auto cv : mesh.tet_vertices(c)) {
        vertices.push_back(parameter_space ? qarameter[cv] : mesh.vertex(cv).cast<mpq_class>());
    }
    Matrix3q m;
    for (int i = 0; i < 3; i++) {
        m.col(i) = vertices[i + 1] - vertices[0];
    }
    return m.determinant() <= 0;
}

Eigen::Vector3d kernel_chebyshev_center(Mesh& mesh, std::vector<OpenVolumeMesh::HalfFaceHandle> polyhedron) {
    CGAL::Quadratic_program<double> lp(CGAL::SMALLER, false);
    lp.set_c(3, -1);
    OpenVolumeMesh::VertexPropertyT<Vector3q> qarameter = mesh.request_vertex_property<Vector3q>("qarameter");
    for (int i = 0; i < polyhedron.size(); i++) {
        OpenVolumeMesh::HalfFaceHandle h = polyhedron[i];
        std::vector<Vector3q> triangle;
        for (auto hv : mesh.halfface_vertices(h)) {
            triangle.push_back(qarameter[hv]);
        }
        Vector3q n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
        Eigen::Vector3d n_normalized = n.unaryExpr([](mpq_class x) { return x.get_d(); }).normalized(); // (0, 0, 0) stays (0, 0, 0)
        // n_0 * x + n_1 * y + n_2 * z + r <= -d
        for (int j = 0; j < 3; j++) {
            lp.set_a(j, i, n_normalized(j));
        }
        lp.set_a(3, i, 1);
        lp.set_b(i, n_normalized.cast<mpq_class>().dot(triangle[0]).get_d());
    }
    CGAL::Quadratic_program_solution<double> sol = CGAL::solve_linear_program(lp, 21.0);
    Eigen::Vector3d x;
    auto it = sol.variable_values_begin();
    for (int i = 0; i < 3; i++) {
        x(i) = CGAL::quotient_truncation(*it);
        it++;
    }
    if (!sol.is_optimal() || std::isnan(x(0)) || std::isnan(x(1)) || std::isnan(x(2))) {
        Vector3q center = Vector3q::Zero();
        int cnt = 0;
        for (auto h : polyhedron) {
            for (auto hv : mesh.halfface_vertices(h)) {
                center += qarameter[hv];
                cnt++;
            }
        }
        center /= cnt;
        return center.unaryExpr([](mpq_class x) { return x.get_d(); });
    }
    return x;
}

void decimate_refined(Mesh& mesh) {
    OpenVolumeMesh::VertexPropertyT<Vector3q> qarameter = mesh.request_vertex_property<Vector3q>("qarameter");
    OpenVolumeMesh::VertexPropertyT<std::vector<std::pair<int, mpq_class>>> bary = mesh.request_vertex_property<std::vector<std::pair<int, mpq_class>>>("vertex");
    std::vector<OpenVolumeMesh::HalfEdgeHandle> halfedges;
    for (auto e : mesh.halfedges()) {
        halfedges.push_back(e);
    }
    std::default_random_engine rng;
    rng.seed(6364);
    std::shuffle(halfedges.begin(), halfedges.end(), rng);
    std::queue<OpenVolumeMesh::HalfEdgeHandle> q;
    for (auto e : halfedges) {
        q.push(e);
    }
    while (!q.empty()) {
        OpenVolumeMesh::HalfEdgeHandle e = q.front();
        q.pop();
        if (mesh.is_deleted(e)) {
            continue;
        }
        OpenVolumeMesh::VertexHandle u = mesh.from_vertex_handle(e);
        OpenVolumeMesh::VertexHandle v = mesh.to_vertex_handle(e);

        // Barycentric coordinate entries of v must be a subset of those of u
        bool collapse_ok = true;
        for (int i = 0; i < bary[v].size() && collapse_ok; i++) {
            int entry = bary[v][i].first;
            bool found = false;
            for (int j = 0; j < bary[u].size() && !found; j++) {
                if (bary[u][j].first == entry) {
                    found = true;
                }
            }
            if (!found) {
                collapse_ok = false;
            }
        }
        if (!collapse_ok) {
            continue;
        }

        // Simulate halfedge collapse, tets in the back must not be inverted
        Eigen::Vector3d position_before = mesh.vertex(u);
        Vector3q qarameter_before = qarameter[u];
        mesh.set_vertex(u, mesh.vertex(v));
        qarameter[u] = qarameter[v];
        for (auto uc : mesh.vertex_cells(u)) {
            bool back = true;
            for (auto ec : mesh.halfedge_cells(e)) {
                if (ec == uc) {
                    back = false;
                    break;
                }
            }
            if (back && (inverted(mesh, uc, true) || inverted(mesh, uc, false))) {
                collapse_ok = false;
                break;
            }
        }
        mesh.set_vertex(u, position_before);
        qarameter[u] = qarameter_before;
        if (!collapse_ok) {
            continue;
        }

        // Update outgoing halfedges of vertices of changed tets
        for (auto uv : mesh.vertex_vertices(u)) {
            for (auto uve : mesh.outgoing_halfedges(uv)) {
                q.push(uve);
            }
        }
        // Halfedge collapse
        mesh.collapse_edge(e);
    }
    mesh.collect_garbage();
}

bool galaxy_map(Mesh& mesh) {
    OpenVolumeMesh::CellPropertyT<Eigen::Vector3d> color_c = mesh.request_cell_property<Eigen::Vector3d>("color_c", Eigen::Vector3d{1, 1, 1});
    mesh.set_persistent(color_c);
    OpenVolumeMesh::FacePropertyT<Eigen::Vector3d> color_f = mesh.request_face_property<Eigen::Vector3d>("color_f", Eigen::Vector3d{1, 1, 1});
    mesh.set_persistent(color_f);
    OpenVolumeMesh::EdgePropertyT<Eigen::Vector3d> color_e = mesh.request_edge_property<Eigen::Vector3d>("color_e", Eigen::Vector3d{1, 1, 1});
    mesh.set_persistent(color_e);
    OpenVolumeMesh::VertexPropertyT<Eigen::Vector3d> color_v = mesh.request_vertex_property<Eigen::Vector3d>("color_v", Eigen::Vector3d{1, 1, 1});
    mesh.set_persistent(color_v);

    bool parameter_exists = mesh.vertex_property_exists<Eigen::Vector3d>("parameter");
    bool qarameter_exists = mesh.vertex_property_exists<Vector3q>("qarameter");
    if (!parameter_exists && !qarameter_exists) {
        std::cout << "Parameterization missing" << std::endl;
        return false;
    }
    OpenVolumeMesh::VertexPropertyT<Eigen::Vector3d> parameter = mesh.request_vertex_property<Eigen::Vector3d>("parameter");
    OpenVolumeMesh::VertexPropertyT<Vector3q> qarameter = mesh.request_vertex_property<Vector3q>("qarameter");
    if (!parameter_exists) {
        mesh.set_persistent(parameter);
        for (auto v : mesh.vertices()) {
            parameter[v] = qarameter[v].unaryExpr([](mpq_class x) { return x.get_d(); });
        }
    }
    if (!qarameter_exists) {
        mesh.set_persistent(qarameter);
        for (auto v : mesh.vertices()) {
            qarameter[v] = parameter[v].cast<mpq_class>();
        }
    }

    for (auto c : mesh.cells()) {
        if (inverted(mesh, c)) {
            color_c[c] = {1, 0, 0.2};
        }
    }

    // Grow stars simultaneously
    std::vector<std::vector<OpenVolumeMesh::CellHandle>> stars;
    OpenVolumeMesh::CellPropertyT<int> star_idx = mesh.request_cell_property<int>("star_idx", -1);
    std::vector<Vector3q> centers;
    std::vector<bool> alive;
    int i = 0;

    for (auto cell : mesh.cells()) {
        if (inverted(mesh, cell) && star_idx[cell] == -1) {
            // Grow star around inverted tet
            std::vector<OpenVolumeMesh::CellHandle> star;

            std::function<void(OpenVolumeMesh::CellHandle)> add_tet_to_star = [&](OpenVolumeMesh::CellHandle c) {
                star.push_back(c);
                star_idx[c] = stars.size();
                color_c[c] = {0.8, 0, 1};
            };

            Vector3q p;
            OpenVolumeMesh::CellHandle next = cell;
            bool star_shaped = false;
            while (!star_shaped) {
                // Add next tet and necessary neighbors
                std::queue<OpenVolumeMesh::CellHandle> cc_queue;
                OpenVolumeMesh::CellPropertyT<bool> in_cc_queue = mesh.request_cell_property<bool>("in_cc_queue");
                in_cc_queue[next] = true;
                cc_queue.push(next);
                while (!cc_queue.empty()) {
                    OpenVolumeMesh::CellHandle c = cc_queue.front();
                    cc_queue.pop();
                    int star_idx_c = star_idx[c];
                    if (star_idx_c != stars.size()) {
                        if (star_idx_c == -1) {
                            add_tet_to_star(c);
                        } else {
                            // Gobble up colliding star
                            for (auto collision_c : stars[star_idx_c]) {
                                add_tet_to_star(collision_c);
                            }
                            alive[star_idx_c] = false;
                        }
                        // Add inverted elements next to the boundary to prevent self intersections
                        for (auto cc : mesh.cell_cells(c)) {
                            if (star_idx[cc] != stars.size() && !in_cc_queue[cc] && inverted(mesh, cc)) {
                                in_cc_queue[cc] = true;
                                cc_queue.push(cc);
                            }
                        }
                    }
                }

                if (star.size() == 1) {
                    // First select inverted-encompassing tet if it exists (impossible with Tutte)
                    star_shaped = false;
                    OpenVolumeMesh::CellHandle c = star[0];
                    next = *(mesh.cc_iter(c));
                    for (auto ch : mesh.cell_halffaces(c)) {
                        OpenVolumeMesh::VertexHandle cv = mesh.halfface_opposite_vertex(ch);
                        OpenVolumeMesh::HalfFaceHandle oh = mesh.opposite_halfface_handle(ch);
                        if (mesh.is_boundary(oh)) {
                            continue;
                        }
                        OpenVolumeMesh::VertexHandle ov = mesh.halfface_opposite_vertex(oh);
                        bool encompassing = true;
                        for (auto cvh : mesh.cell_halffaces(c)) {
                            if (cvh == ch) {
                                continue;
                            }
                            std::vector<Vector3q> triangle;
                            for (auto cvhv : mesh.halfface_vertices(cvh)) {
                                triangle.push_back(qarameter[cvhv]);
                            }
                            Vector3q n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
                            if (n.dot(qarameter[ov]) - n.dot(triangle[0]) <= 0) {
                                encompassing = false;
                            }
                        }
                        if (encompassing) {
                            next = mesh.incident_cell(oh);
                            break;
                        }
                    }
                    continue;
                }
                // Get next tet by heuristic
                std::vector<OpenVolumeMesh::HalfFaceHandle> star_boundary;
                for (auto c : star) {
                    for (auto ch : mesh.cell_halffaces(c)) {
                        OpenVolumeMesh::HalfFaceHandle h = mesh.opposite_halfface_handle(ch);
                        if (mesh.is_boundary(h) || star_idx[mesh.incident_cell(h)] != stars.size()) {
                            star_boundary.push_back(h);
                        }
                    }
                }
                Eigen::Vector3d pd = kernel_chebyshev_center(mesh, star_boundary);
                p = pd.cast<mpq_class>();
                star_shaped = true;
                double max = std::numeric_limits<double>::lowest();
                OpenVolumeMesh::HalfFaceHandle max_h = Mesh::InvalidHalfFaceHandle;
                for (auto c : star) {
                    for (auto ch : mesh.cell_halffaces(c)) {
                        OpenVolumeMesh::HalfFaceHandle h = mesh.opposite_halfface_handle(ch);
                        if (mesh.is_boundary(h) || star_idx[mesh.incident_cell(h)] != stars.size()) {
                            std::vector<Vector3q> triangle;
                            for (auto hv : mesh.halfface_vertices(h)) {
                                triangle.push_back(qarameter[hv]);
                            }
                            Vector3q n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
                            if (n.dot(p) - n.dot(triangle[0]) >= 0) {
                                star_shaped = false;
                            }
                            Eigen::Vector3d n_normalized = n.unaryExpr([](mpq_class x) { return x.get_d(); }).normalized();
                            double distance = n_normalized.dot(pd) - n_normalized.dot(parameter[*mesh.hfv_iter(h)]);
                            if (distance > max && !mesh.is_boundary(h)) {
                                max = distance;
                                max_h = h;
                            }
                        }
                    }
                }
                if (!star_shaped) {
                    if (!max_h.is_valid()) {
                        std::cout << "Star cannot grow (domain not star-shaped)" << std::endl;
                        return false;
                    }
                    next = mesh.incident_cell(max_h);
                }
            }
            std::cout << "Grown star with " << star.size() << " tetrahedra" << std::endl;

            stars.push_back(star);
            centers.push_back(p);
            alive.push_back(true);
            i++;
        }
    }

    std::cout << "Created galaxy" << std::endl;

    // Reparameterize stars using foliations
    for (i = 0; i < stars.size(); i++) {
        if (alive[i]) {
            // Cut out star
            Mesh cutout;
            OpenVolumeMesh::VertexPropertyT<Vector3q> cutout_qarameter = cutout.request_vertex_property<Vector3q>("qarameter");
            cutout.set_persistent(cutout_qarameter);
            OpenVolumeMesh::VertexPropertyT<OpenVolumeMesh::VertexHandle> v_c2m = cutout.request_vertex_property<OpenVolumeMesh::VertexHandle>("v_c2m");
            OpenVolumeMesh::HalfFacePropertyT<OpenVolumeMesh::HalfFaceHandle> h_c2m = cutout.request_halfface_property<OpenVolumeMesh::HalfFaceHandle>("h_c2m");
            std::map<OpenVolumeMesh::VertexHandle, OpenVolumeMesh::VertexHandle> v_m2c;
            std::map<OpenVolumeMesh::HalfFaceHandle, OpenVolumeMesh::HalfFaceHandle> h_m2c;

            // Adds a tet from mesh to cutout maintaining the maps
            std::function<bool(OpenVolumeMesh::CellHandle)> add_tet_to_cutout = [&](OpenVolumeMesh::CellHandle c) {
                for (auto cv : mesh.cell_vertices(c)) {
                    if (v_m2c.count(cv) == 0) {
                        OpenVolumeMesh::VertexHandle cutout_v = cutout.add_vertex(mesh.vertex(cv));
                        cutout_qarameter[cutout_v] = qarameter[cv];
                        v_c2m[cutout_v] = cv;
                        v_m2c[cv] = cutout_v;
                    }
                }
                std::vector<OpenVolumeMesh::HalfFaceHandle> halffaces;
                for (auto ch : mesh.cell_halffaces(c)) {
                    std::vector<OpenVolumeMesh::VertexHandle> vertices;
                    for (auto chv : mesh.halfface_vertices(ch)) {
                        vertices.push_back(v_m2c[chv]);
                    }
                    OpenVolumeMesh::HalfFaceHandle cutout_h = cutout.add_halfface(vertices[0], vertices[1], vertices[2], true);
                    h_c2m[cutout_h] = ch;
                    h_m2c[ch] = cutout_h;
                    h_c2m[cutout.opposite_halfface_handle(cutout_h)] = mesh.opposite_halfface_handle(ch);
                    h_m2c[mesh.opposite_halfface_handle(ch)] = cutout.opposite_halfface_handle(cutout_h);
                    halffaces.push_back(cutout_h);
                }
                OpenVolumeMesh::CellHandle cutout_c = cutout.add_cell(halffaces, true);
                return cutout_c.is_valid();
            };

            std::vector<OpenVolumeMesh::CellHandle> star;
            for (auto c : stars[i]) {
                if (!mesh.is_deleted(c)) {
                    star.push_back(c);
                }
            }
            for (auto c : star) {
                add_tet_to_cutout(c);
            }

            std::cout << "Processing star with " << star.size() << " tetrahedra" << std::endl;
            if (star.size() > 1000) {
                std::cout << "Long runtime expected" << std::endl;
            }

            // Parameterize cutout using a foliation
            FoliationParameterization foliation_para(cutout);
            for (auto cutout_v : cutout.boundary_vertices()) {
                cutout_qarameter[cutout_v] = cutout_qarameter[cutout_v] - centers[i];
            }
            if (!foliation_para.set_surface_parameterization()) {
                std::cout << "\nSurface parameterization for foliation parameterization not star-shaped" << std::endl;
                return false;
            }
            if (!foliation_para.generate_directions()) {
                std::cout << "\nFoliation parameterization directions not generated" << std::endl;
                return false;
            }

            for (auto c : star) {
                mesh.remove_tet(c);
            }

            // Refine cutout
            Mesh refined = foliation_para.refine();
            OpenVolumeMesh::VertexPropertyT<Vector3q> refined_qarameter = refined.request_vertex_property<Vector3q>("qarameter");
            OpenVolumeMesh::VertexPropertyT<std::vector<std::pair<int, mpq_class>>> refined_bary = refined.request_vertex_property<std::vector<std::pair<int, mpq_class>>>("vertex");
            for (auto refined_v : refined.vertices()) {
                refined_qarameter[refined_v] += centers[i];
            }
            decimate_refined(refined);

            // Split tets incident to cutout such that the link can be filled
            std::set<OpenVolumeMesh::CellHandle> link_cells;
            for (auto cutout_v : cutout.boundary_vertices()) {
                OpenVolumeMesh::VertexHandle v = v_c2m[cutout_v];
                for (auto vc : mesh.vertex_cells(v)) {
                    link_cells.insert(vc);
                }
            }
            OpenVolumeMesh::CellPropertyT<std::vector<OpenVolumeMesh::CellHandle>> children = mesh.request_cell_property<std::vector<OpenVolumeMesh::CellHandle>>("children");

            // Recursively splits a tet until there are no two cutout edges remaining and saves the children
            std::function<void(OpenVolumeMesh::CellHandle)> split = [&](OpenVolumeMesh::CellHandle c) {
                int n_1_cutout_edge_faces = 0;
                for (auto ch : mesh.cell_halffaces(c)) {
                    if (h_m2c.count(ch) == 0) {
                        int n_cutout_edges = 0;
                        for (auto che : mesh.halfface_edges(ch)) {
                            bool cutout_edge = false;
                            for (auto cheh : mesh.edge_halffaces(che)) {
                                if (h_m2c.count(cheh) == 1) {
                                    cutout_edge = true;
                                    break;
                                }
                            }
                            if (cutout_edge) {
                                n_cutout_edges++;
                            }
                        }
                        if (n_cutout_edges > 1) {
                            // Split non-cutout face with >1 cutout edges
                            OpenVolumeMesh::CellHandle c2 = mesh.incident_cell(mesh.opposite_halfface_handle(ch));
                            OpenVolumeMesh::VertexHandle h = mesh.halfface_opposite_vertex(ch);
                            OpenVolumeMesh::FaceHandle split_f = mesh.face_handle(ch);
                            Eigen::Vector3d v_position = Eigen::Vector3d::Zero();
                            Vector3q v_qarameter = Vector3q::Zero();
                            for (auto fv : mesh.face_vertices(split_f)) {
                                v_position += mesh.vertex(fv);
                                v_qarameter += qarameter[fv];
                            }
                            v_position /= 3;
                            v_qarameter /= 3;
                            OpenVolumeMesh::VertexHandle v = mesh.split_face(split_f, v_position); // Cell properties are copied
                            parameter[v] = v_qarameter.unaryExpr([](mpq_class x) { return x.get_d(); });
                            qarameter[v] = v_qarameter;
                            for (auto vc : mesh.vertex_cells(v)) {
                                bool has_h = false;
                                for (auto vcv : mesh.cell_vertices(vc)) {
                                    if (vcv == h) {
                                        has_h = true;
                                        break;
                                    }
                                }
                                if (has_h) {
                                    children[c].push_back(vc);
                                } else if (c2.is_valid()) {
                                    children[c2].push_back(vc);
                                }
                                if (star_idx[vc] != -1) {
                                    stars[star_idx[vc]].push_back(c);
                                }
                            }
                            // Split children recursively
                            for (auto child : children[c]) {
                                split(child);
                            }
                            return;
                        }
                        if (n_cutout_edges == 1) {
                            n_1_cutout_edge_faces++;
                        }
                    }
                }
                if (n_1_cutout_edge_faces == 4) {
                    // Split tet with two opposite cutout edges (Case 19)
                    Eigen::Vector3d v_position = Eigen::Vector3d::Zero();
                    Vector3q v_qarameter = Vector3q::Zero();
                    for (auto cv : mesh.cell_vertices(c)) {
                        v_position += mesh.vertex(cv);
                        v_qarameter += qarameter[cv];
                    }
                    v_position /= 4;
                    v_qarameter /= 4;
                    int star_idx_c = star_idx[c];
                    OpenVolumeMesh::VertexHandle v = mesh.split_tet(c, v_position);
                    parameter[v] = v_qarameter.unaryExpr([](mpq_class x) { return x.get_d(); });
                    qarameter[v] = v_qarameter;
                    for (auto vc : mesh.vertex_cells(v)) {
                        children[c].push_back(vc);
                        if (star_idx_c != -1) {
                            star_idx[vc] = star_idx_c;
                            stars[star_idx_c].push_back(c);
                        }
                    }
                }
            };

            // Calls the split method with a tet or its children
            std::function<void(OpenVolumeMesh::CellHandle)> split_init = [&](OpenVolumeMesh::CellHandle c) {
                if (children[c].size() == 0) {
                    split(c);
                } else {
                    for (auto child : children[c]) {
                        split_init(child);
                    }
                }
            };

            for (auto c : link_cells) {
                split_init(c);
            }

            // Add refined vertices
            std::map<OpenVolumeMesh::VertexHandle, OpenVolumeMesh::VertexHandle> v_c2r;
            OpenVolumeMesh::VertexPropertyT<OpenVolumeMesh::VertexHandle> v_r2n = refined.request_vertex_property<OpenVolumeMesh::VertexHandle>("v_r2n");
            for (auto refined_v : refined.vertices()) {
                std::vector<std::pair<int, mpq_class>> bary = refined_bary[refined_v];
                if (bary.size() == 1) {
                    v_c2r[OpenVolumeMesh::VertexHandle(bary[0].first)] = refined_v;
                }
                OpenVolumeMesh::VertexHandle v = mesh.add_vertex(refined.vertex(refined_v));
                parameter[v] = refined_qarameter[refined_v].unaryExpr([](mpq_class x) { return x.get_d(); });
                qarameter[v] = refined_qarameter[refined_v];
                v_r2n[refined_v] = v;
            }

            // Add refined tets
            for (auto refined_c : refined.cells()) {
                std::vector<OpenVolumeMesh::VertexHandle> tet;
                for (auto refined_cv : refined.tet_vertices(refined_c)) {
                    tet.push_back(v_r2n[refined_cv]);
                }
                OpenVolumeMesh::CellHandle c = mesh.add_cell(tet, true);
                color_c[c] = {0, 0.2, 1};
                if (inverted(mesh, c)) {
                    std::cout << "\nRefined tet inverted" << std::endl;
                    return false;
                }
            }

            // Fill link between mesh and refined
            // Add link tets from faces
            for (auto refined_h : refined.boundary_halffaces()) {
                std::vector<OpenVolumeMesh::VertexHandle> tet;
                std::vector<int> bary_indices;
                for (auto refined_hv : refined.halfface_vertices(refined_h)) {
                    tet.push_back(v_r2n[refined_hv]);
                    for (std::pair<int, mpq_class> bary : refined_bary[refined_hv]) {
                        if (std::find(bary_indices.begin(), bary_indices.end(), bary.first) == bary_indices.end()) {
                            bary_indices.push_back(bary.first);
                        }
                    }
                }
                std::vector<OpenVolumeMesh::VertexHandle> cutout_vertices(3);
                for (int j = 0; j < 3; j++) {
                    cutout_vertices[j] = OpenVolumeMesh::VertexHandle(bary_indices[j]);
                }
                OpenVolumeMesh::HalfFaceHandle cutout_h = cutout.halfface(cutout_vertices);
                if (!cutout.is_boundary(cutout_h)) {
                    cutout_h = cutout.opposite_halfface_handle(cutout_h);
                }
                OpenVolumeMesh::HalfFaceHandle h = h_c2m[cutout_h];
                if (!mesh.is_boundary(h)) {
                    OpenVolumeMesh::VertexHandle v = mesh.halfface_opposite_vertex(h);
                    if (v_m2c.count(v) == 1) {
                        v = v_r2n[v_c2r[v_m2c[v]]];
                    }
                    tet.push_back(v);
                    OpenVolumeMesh::CellHandle c = mesh.add_cell(tet, true);
                    OpenVolumeMesh::CellHandle hc = mesh.incident_cell(h);
                    if (star_idx[hc] != -1) {
                        star_idx[c] = star_idx[hc];
                        stars[star_idx[hc]].push_back(c);
                        color_c[c] = {0.8, 0, 1};
                    }
                }
            }
            // Add link tets from edges
            for (auto refined_e : refined.boundary_edges()) {
                std::vector<OpenVolumeMesh::VertexHandle> tet_ab;
                std::vector<int> bary_indices;
                std::vector<mpq_class> alpha(2);
                int j = 0;
                for (auto refined_ev : refined.edge_vertices(refined_e)) {
                    tet_ab.push_back(v_r2n[refined_ev]);
                    alpha[j] = 0;
                    for (std::pair<int, mpq_class> bary : refined_bary[refined_ev]) {
                        if (std::find(bary_indices.begin(), bary_indices.end(), bary.first) == bary_indices.end()) {
                            bary_indices.push_back(bary.first);
                        }
                        if (bary.first == bary_indices[0]) {
                            alpha[j] = bary.second;
                        }
                    }
                    j++;
                }
                if (bary_indices.size() == 2) {
                    std::vector<OpenVolumeMesh::VertexHandle> cutout_vertices(2);
                    for (j = 0; j < 2; j++) {
                        cutout_vertices[j] = OpenVolumeMesh::VertexHandle(bary_indices[alpha[0] > alpha[1] ? j : 1 - j]);
                    }
                    OpenVolumeMesh::HalfEdgeHandle e = mesh.halfedge(v_c2m[cutout_vertices[0]], v_c2m[cutout_vertices[1]]);
                    for (auto ec : mesh.halfedge_cells(e)) {
                        bool face_connected = false;
                        for (auto ech : mesh.cell_halffaces(ec)) {
                            if (h_m2c.count(ech) == 1) {
                                face_connected = true;
                                break;
                            }
                        }
                        if (!face_connected) {
                            OpenVolumeMesh::HalfEdgeHandle e2 = mesh.halfedge_opposite_halfedge(e, ec);
                            OpenVolumeMesh::VertexHandle v0 = mesh.from_vertex_handle(e2);
                            if (v_m2c.count(v0) == 1) {
                                v0 = v_r2n[v_c2r[v_m2c[v0]]];
                            }
                            OpenVolumeMesh::VertexHandle v1 = mesh.to_vertex_handle(e2);
                            if (v_m2c.count(v1) == 1) {
                                v1 = v_r2n[v_c2r[v_m2c[v1]]];
                            }
                            std::vector<OpenVolumeMesh::VertexHandle> tet = tet_ab;
                            tet.push_back(v0);
                            tet.push_back(v1);
                            OpenVolumeMesh::CellHandle c = mesh.add_cell(tet, true);
                            if (c.is_valid()) {
                                if (star_idx[ec] != -1) {
                                    star_idx[c] = star_idx[ec];
                                    stars[star_idx[ec]].push_back(c);
                                    color_c[c] = {0.8, 0, 1};
                                }
                            }
                        }
                    }
                }
            }

            // Save link tets from vertices
            std::vector<std::pair<std::vector<OpenVolumeMesh::VertexHandle>, int>> tets_and_star_indices;
            for (auto refined_v : refined.boundary_vertices()) {
                std::vector<std::pair<int, mpq_class>> bary = refined_bary[refined_v];
                if (bary.size() == 1) {
                    OpenVolumeMesh::VertexHandle v = v_c2m[OpenVolumeMesh::VertexHandle(bary[0].first)];
                    for (auto vc : mesh.vertex_cells(v)) {
                        bool edge_connected = false;
                        for (auto vce : mesh.cell_edges(vc)) {
                            for (auto vceh : mesh.edge_halffaces(vce)) {
                                if (h_m2c.count(vceh) == 1) {
                                    edge_connected = true;
                                    break;
                                }
                            }
                            if (edge_connected) {
                                break;
                            }
                        }
                        if (!edge_connected) {
                            OpenVolumeMesh::HalfFaceHandle h = mesh.vertex_opposite_halfface(v, vc);
                            std::vector<OpenVolumeMesh::VertexHandle> tet;
                            for (auto hv : mesh.halfface_vertices(h)) {
                                if (v_m2c.count(hv) == 1) {
                                    hv = v_r2n[v_c2r[v_m2c[hv]]];
                                }
                                tet.push_back(hv);
                            }
                            tet.push_back(v_r2n[refined_v]);
                            tets_and_star_indices.push_back({tet, star_idx[vc]});
                        }
                    }
                }
            }

            // Delete old tets in the link
            for (auto cutout_h : cutout.halffaces()) {
                for (auto hv : mesh.halfface_vertices(h_c2m[cutout_h])) {
                    if (!mesh.is_deleted(hv)) {
                        mesh.delete_vertex(hv);
                    }
                }
            }

            // Add saved tets
            for (auto [tet, tet_star_idx] : tets_and_star_indices) {
                OpenVolumeMesh::CellHandle c = mesh.add_cell(tet, true);
                if (c.is_valid()) {
                    if (tet_star_idx != -1) {
                        star_idx[c] = tet_star_idx;
                        stars[tet_star_idx].push_back(c);
                        color_c[c] = {0.8, 0, 1};
                    }
                }
            }
        }
    }

    mesh.collect_garbage();
    std::cout << "Finished" << std::endl;
    return true;
}
