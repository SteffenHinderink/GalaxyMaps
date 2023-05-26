#include "canvas.h"

nanogui::Matrix4f Canvas::to_nanogui(Eigen::Matrix4d m) {
    nanogui::Matrix4f n;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            // Column major
            n.m[j][i] = m(i, j);
        }
    }
    return n;
}

// Eigen::Vector3d Canvas::hsv_to_rgb(Eigen::Vector3d hsv) {
//     int h = (int) (hsv(0) / (M_PI / 3)) % 6;
//     double f = hsv(0) / (M_PI / 3) - h;
//     double p = hsv(2) * (1 - hsv(1));
//     double q = hsv(2) * (1 - hsv(1) * f);
//     double t = hsv(2) * (1 - hsv(1) * (1 - f));
//     switch (h) {
//     case 0:
//         return {hsv(2), t, p};
//     case 1:
//         return {q, hsv(2), p};
//     case 2:
//         return {p, hsv(2), t};
//     case 3:
//         return {p, q, hsv(2)};
//     case 4:
//         return {t, p, hsv(2)};
//     case 5:
//         return {hsv(2), p, q};
//     }
//     return {0, 0, 0};
// }

Canvas::Canvas(std::string dir, nanogui::Widget* parent, Mesh& mesh, bool parameter_space)
    : nanogui::Canvas(parent), mesh_(mesh), parameter_space_(parameter_space), last_frame_(std::chrono::steady_clock::now()), rot_(false), rot_fps_(false), trans_(false),
      size_(0.5), show_cells_(true), show_faces_(false), show_edges_(false), show_vertices_(false), show_boundary_(false), show_coords_(false), show_colored_(false), white_background_(false), select_(false) {
    {
        std::lock_guard sync(mutex_);

        std::ifstream file(dir + "../shader/vs.vert");
        std::string vs((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
        file = std::ifstream(dir + "../shader/fs.frag");
        std::string fs((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
        shader_ = new nanogui::Shader(render_pass(), "shader", vs, fs);
        this->set_fixed_size({SIZE, SIZE});
        if (!parameter_space_) {
            request_focus();
        }
        selected_c_ = mesh_.request_cell_property<bool>("selected_c");
        selected_f_ = mesh_.request_face_property<bool>("selected_f");
        selected_e_ = mesh_.request_edge_property<bool>("selected_e");
        selected_v_ = mesh_.request_vertex_property<bool>("selected_v");
    }
    init();
}

void Canvas::set_other(Canvas* other) {
    other_ = other;
}

void Canvas::init() {
    {
        std::lock_guard sync(mutex_);

        Eigen::Vector3d min = Eigen::Vector3d::Constant(std::numeric_limits<double>::max());
        Eigen::Vector3d max = Eigen::Vector3d::Constant(std::numeric_limits<double>::lowest());
        OpenVolumeMesh::VertexPropertyT<Eigen::Vector3d> parameter = mesh_.request_vertex_property<Eigen::Vector3d>("parameter");
        for (auto v : mesh_.vertices()) {
            Eigen::Vector3d p = parameter_space_ ? parameter[v] : mesh_.vertex(v);
            for (int i = 0; i < 3; i++) {
                if (p(i) < min(i)) {
                    min(i) = p(i);
                }
                if (p(i) > max(i)) {
                    max(i) = p(i);
                }
            }
        }
        scale_ = 0;
        for (int i = 0; i < 3; i++) {
            if (max(i) - min(i) > scale_) {
                scale_ = max(i) - min(i);
            }
        }
        scale_ *= 1.1;
        double near = scale_ / 100;
        double far = scale_ * 100;

        Eigen::Matrix4d projection_matrix = Eigen::Matrix4d::Identity();
        projection_matrix.block<2, 2>(2, 2) << (-far - near) / (far - near), -2 * far * near / (far - near), -1, 0;
        shader_->set_uniform("uProjectionMatrix", to_nanogui(projection_matrix));
        Eigen::Vector3d c = (min + max) / 2;
        transformation_matrix_ = Eigen::Matrix4d::Identity();
        transformation_matrix_.block<3, 1>(0, 3) = -c;
    }
    update();
}

void Canvas::update() {
    if (!mutex_.try_lock()) {
        return;
    }
    std::lock_guard sync(update_mutex_);

    if (white_background_) {
        this->set_background_color({255, 255, 255, 255});
    } else {
        this->set_background_color({32, 32, 32, 255});
    }

    indices_.clear();
    positions_.clear();
    normals_.clear();
    colors_.clear();
    sources_.clear();

    auto add_triangle = [&](std::vector<Eigen::Vector3d> tri_positions, std::vector<Eigen::Vector3d> tri_colors, std::pair<Element, int> source = {Element::OTHER, -1}, std::vector<Eigen::Vector3d> tri_normals = {}) {
        if (tri_normals.size() == 0) {
            tri_normals = std::vector<Eigen::Vector3d>(3, (tri_positions[1] - tri_positions[0]).cross(tri_positions[2] - tri_positions[0]).normalized());
        }
        for (int i = 0; i < 3; i++) {
            indices_.push_back(indices_.size());
            for (int j = 0; j < 3; j++) {
                positions_.push_back(tri_positions[i](j));
                normals_.push_back(tri_normals[i](j));
                colors_.push_back(tri_colors[i](j));
            }
        }
        sources_.push_back(source);
    };

    auto add_sphere = [&](Eigen::Vector3d c, double r, Eigen::Vector3d color, std::pair<Element, int> source = {Element::OTHER, -1}, int res = 2) {
        double phi = 1.61803398875; // = (1 + sqrt(5)) / 2 (Golden ratio)
        std::vector<Eigen::Vector3d> vertices = {
            {0, 1, phi},
            {0, 1, -phi},
            {0, -1, phi},
            {0, -1, -phi},
            {1, phi, 0},
            {1, -phi, 0},
            {-1, phi, 0},
            {-1, -phi, 0},
            {phi, 0, 1},
            {-phi, 0, 1},
            {phi, 0, -1},
            {-phi, 0, -1}};
        std::vector<std::vector<int>> faces = {
            {0, 4, 6},
            {4, 1, 6},
            {6, 1, 11},
            {1, 3, 11},
            {11, 3, 7},
            {3, 5, 7},
            {7, 5, 2},
            {5, 8, 2},
            {2, 8, 0},
            {8, 4, 0},
            {9, 0, 6},
            {9, 6, 11},
            {9, 11, 7},
            {9, 7, 2},
            {9, 2, 0},
            {10, 4, 8},
            {10, 8, 5},
            {10, 5, 3},
            {10, 3, 1},
            {10, 1, 4}};
        for (int i = 0; i < faces.size(); i++) {
            std::vector<Eigen::Vector3d> triangle(3);
            for (int j = 0; j < 3; j++) {
                triangle[j] = vertices[faces[i][j]];
            }
            Eigen::Vector3d ab = (triangle[1] - triangle[0]) / res;
            Eigen::Vector3d ac = (triangle[2] - triangle[0]) / res;
            std::vector<std::vector<Eigen::Vector3d>> p(res + 1);
            std::vector<std::vector<Eigen::Vector3d>> n(res + 1);
            for (int j = 0; j <= res; j++) {
                p[j] = std::vector<Eigen::Vector3d>(res - j + 1);
                n[j] = std::vector<Eigen::Vector3d>(res - j + 1);
                for (int k = 0; k <= res - j; k++) {
                    n[j][k] = (triangle[0] + j * ab + k * ac).normalized();
                    p[j][k] = c + r * n[j][k];
                }
            }
            for (int j = 0; j < res; j++) {
                for (int k = 0; k < res - j; k++) {
                    add_triangle({p[j][k], p[j + 1][k], p[j][k + 1]}, {color, color, color}, source, {n[j][k], n[j + 1][k], n[j][k + 1]});
                    if (k < res - j - 1) {
                        add_triangle({p[j][k + 1], p[j + 1][k], p[j + 1][k + 1]}, {color, color, color}, source, {n[j][k + 1], n[j + 1][k], n[j + 1][k + 1]});
                    }
                }
            }
        }
    };

    auto add_cylinder = [&](Eigen::Vector3d a, Eigen::Vector3d b, double r, Eigen::Vector3d color, std::pair<Element, int> source = {Element::OTHER, -1}, int res = 8) {
        Eigen::Vector3d x = (b - a).normalized();
        Eigen::Vector3d y = x.cross(std::abs(x(0)) > 0.999 ? Eigen::Vector3d(0, 1, 0) : Eigen::Vector3d(1, 0, 0)).normalized();
        Eigen::Vector3d z = x.cross(y).normalized();
        Eigen::Vector3d n0 = z;
        Eigen::Vector3d a0 = a + r * n0;
        Eigen::Vector3d b0 = b + r * n0;
        for (int i = 0; i < res; i++) {
            double angle = 2 * M_PI * (i + 1) / res;
            Eigen::Vector3d n1 = std::sin(angle) * y + std::cos(angle) * z;
            Eigen::Vector3d a1 = a + r * n1;
            Eigen::Vector3d b1 = b + r * n1;
            add_triangle({a, a0, a1}, {color, color, color}, source);
            add_triangle({a0, b0, b1}, {color, color, color}, source, {n0, n0, n1});
            add_triangle({a0, b1, a1}, {color, color, color}, source, {n0, n1, n1});
            add_triangle({b, b1, b0}, {color, color, color}, source);
            n0 = n1;
            a0 = a1;
            b0 = b1;
        }
    };

    auto add_cone = [&](Eigen::Vector3d a, Eigen::Vector3d b, double r, Eigen::Vector3d color, std::pair<Element, int> source = {Element::OTHER, -1}, int res = 16) {
        Eigen::Vector3d ab = b - a;
        double l = ab.norm();
        Eigen::Vector3d x = ab / l;
        Eigen::Vector3d y = x.cross(std::abs(x(0)) > 0.999 ? Eigen::Vector3d(0, 1, 0) : Eigen::Vector3d(1, 0, 0)).normalized();
        Eigen::Vector3d z = x.cross(y).normalized();
        Eigen::Vector3d n0 = (l * z + r * x).normalized();
        Eigen::Vector3d a0 = a + r * z;
        Eigen::Vector3d b0 = b + r * z;
        for (int i = 0; i < res; i++) {
            double half_angle = 2 * M_PI * (i + 0.5) / res;
            Eigen::Vector3d n = (l * (std::sin(half_angle) * y + std::cos(half_angle) * z) + r * x).normalized();
            double angle = 2 * M_PI * (i + 1) / res;
            Eigen::Vector3d d = std::sin(angle) * y + std::cos(angle) * z;
            Eigen::Vector3d n1 = (l * d + r * x).normalized();
            Eigen::Vector3d a1 = a + r * d;
            Eigen::Vector3d b1 = b + r * d;
            add_triangle({a, a0, a1}, {color, color, color}, source);
            add_triangle({a0, b, a1}, {color, color, color}, source, {n0, n, n1});
            n0 = n1;
            a0 = a1;
            b0 = b1;
        }
    };

    Eigen::Vector3d selection_color = {0, 1, 0.8};
    OpenVolumeMesh::VertexPropertyT<Eigen::Vector3d> parameter = mesh_.request_vertex_property<Eigen::Vector3d>("parameter");
    if (show_cells_) {
        OpenVolumeMesh::CellPropertyT<Eigen::Vector3d> color_c = mesh_.request_cell_property<Eigen::Vector3d>("color_c", Eigen::Vector3d{1, 1, 1});
        for (auto c : mesh_.cells()) {
            if ((!show_boundary_ || mesh_.is_boundary(c)) && (!show_colored_ || color_c[c] != Eigen::Vector3d{1, 1, 1}) && color_c[c] != Eigen::Vector3d{-1, -1, -1}) {
                Eigen::Vector3d center = Eigen::Vector3d::Zero();
                for (auto cv : mesh_.cell_vertices(c)) {
                    center += parameter_space_ ? parameter[cv] : mesh_.vertex(cv);
                }
                center /= 4;
                for (auto cf : mesh_.cell_faces(c)) {
                    std::vector<Eigen::Vector3d> tri_positions;
                    std::vector<Eigen::Vector3d> tri_colors;
                    for (auto cfv : mesh_.face_vertices(cf)) {
                        Eigen::Vector3d p = parameter_space_ ? parameter[cfv] : mesh_.vertex(cfv);
                        tri_positions.push_back(size_ * p + (1 - size_) * center);
                        Eigen::Vector3d color;
                        if (selected_c_[c]) {
                            color = selection_color;
                        } else {
                            color = color_c[c];
                        }
                        tri_colors.push_back(color);
                    }
                    add_triangle(tri_positions, tri_colors, {Element::CELL, c.idx()});
                }
            }
        }
    }
    if (show_faces_) {
        OpenVolumeMesh::FacePropertyT<Eigen::Vector3d> color_f = mesh_.request_face_property<Eigen::Vector3d>("color_f", Eigen::Vector3d{1, 1, 1});
        for (auto f : mesh_.faces()) {
            if ((!show_boundary_ || mesh_.is_boundary(f)) && (!show_colored_ || color_f[f] != Eigen::Vector3d{1, 1, 1}) && color_f[f] != Eigen::Vector3d{-1, -1, -1}) {
                Eigen::Vector3d center = Eigen::Vector3d::Zero();
                for (auto fv : mesh_.face_vertices(f)) {
                    center += parameter_space_ ? parameter[fv] : mesh_.vertex(fv);
                }
                center /= 3;
                std::vector<Eigen::Vector3d> tri_positions;
                std::vector<Eigen::Vector3d> tri_colors;
                for (auto fv : mesh_.face_vertices(f)) {
                    Eigen::Vector3d p = parameter_space_ ? parameter[fv] : mesh_.vertex(fv);
                    tri_positions.push_back(size_ * p + (1 - size_) * center);
                    Eigen::Vector3d color;
                    if (selected_f_[f]) {
                        color = selection_color;
                    } else {
                        color = color_f[f];
                    }
                    tri_colors.push_back(color);
                }
                add_triangle(tri_positions, tri_colors, {Element::FACE, f.idx()});
            }
        }
    }
    if (show_edges_) {
        OpenVolumeMesh::EdgePropertyT<Eigen::Vector3d> color_e = mesh_.request_edge_property<Eigen::Vector3d>("color_e", Eigen::Vector3d{1, 1, 1});
        for (auto e : mesh_.edges()) {
            if ((!show_boundary_ || mesh_.is_boundary(e)) && (!show_colored_ || color_e[e] != Eigen::Vector3d{1, 1, 1}) && color_e[e] != Eigen::Vector3d{-1, -1, -1}) {
                OpenVolumeMesh::VertexHandle v0 = mesh_.edge_vertices(e)[0];
                OpenVolumeMesh::VertexHandle v1 = mesh_.edge_vertices(e)[1];
                Eigen::Vector3d p0 = parameter_space_ ? parameter[v0] : mesh_.vertex(v0);
                Eigen::Vector3d p1 = parameter_space_ ? parameter[v1] : mesh_.vertex(v1);
                add_cylinder(p0, p1, scale_ / 4096, selected_e_[e] ? selection_color : color_e[e], {Element::EDGE, e.idx()});
            }
        }
    }
    if (show_vertices_) {
        OpenVolumeMesh::VertexPropertyT<Eigen::Vector3d> color_v = mesh_.request_vertex_property<Eigen::Vector3d>("color_v", Eigen::Vector3d{1, 1, 1});
        for (auto v : mesh_.vertices()) {
            if ((!show_boundary_ || mesh_.is_boundary(v)) && (!show_colored_ || color_v[v] != Eigen::Vector3d{1, 1, 1}) && color_v[v] != Eigen::Vector3d{-1, -1, -1}) {
                Eigen::Vector3d p = parameter_space_ ? parameter[v] : mesh_.vertex(v);
                add_sphere(p, scale_ / 1024, selected_v_[v] ? selection_color : color_v[v], {Element::VERTEX, v.idx()});
            }
        }
    }
    if (show_coords_) {
        for (int i = 0; i < 3; i++) {
            Eigen::Vector3d unit = Eigen::Vector3d::Zero();
            unit(i) = 1;
            add_cylinder({0, 0, 0}, unit * 0.75, 0.025, unit, {Element::OTHER, i});
            add_cone(unit * 0.75, unit, 0.075, unit, {Element::OTHER, i});
        }
    }
    // if (!parameter_space_ && mesh_.cell_property_exists<std::vector<int>>("d")) {
    //     OpenVolumeMesh::CellPropertyT<std::vector<int>> d = mesh_.request_cell_property<std::vector<int>>("d");
    //     for (auto c : mesh_.cells()) {
    //         std::vector<int> d_c = d[c];
    //         if (!(d_c[0] == 0 && d_c[1] == 0 && d_c[2] == 0 && d_c[3] == 0)) {
    //             Eigen::Vector3d a = Eigen::Vector3d::Zero();
    //             Eigen::Vector3d b = Eigen::Vector3d::Zero();
    //             int i = 0;
    //             for (auto cv : mesh_.cell_vertices(c)) {
    //                 a += 0.25 * mesh_.vertex(cv);
    //                 b += (0.25 + d_c[i]) * mesh_.vertex(cv);
    //                 i++;
    //             }
    //             double length = scale_ / 32;
    //             Eigen::Vector3d dir = (b - a).normalized() * length / 2;
    //             add_cylinder(a - dir, a + dir / 2, length * 0.025, {1, 1, 1});
    //             add_cone(a + dir / 2, a + dir, length * 0.075, {1, 1, 1});
    //         }
    //     }
    // }

    mutex_.unlock();
}

void Canvas::draw_contents() {
    std::lock_guard sync(update_mutex_);

    std::chrono::steady_clock::time_point curr_frame = std::chrono::steady_clock::now();

    if (focused()) {
        Eigen::Matrix4d translation_matrix = Eigen::Matrix4d::Identity();
        double distance = scale_ * std::chrono::duration_cast<std::chrono::microseconds>(curr_frame - last_frame_).count() / 1000000;
        if (distance > scale_ / 100) {
            distance = scale_ / 100;
        }
        if (glfwGetKey(screen()->glfw_window(), GLFW_KEY_W) == GLFW_PRESS) {
            translation_matrix(2, 3) += distance;
        }
        if (glfwGetKey(screen()->glfw_window(), GLFW_KEY_A) == GLFW_PRESS) {
            translation_matrix(0, 3) += distance;
        }
        if (glfwGetKey(screen()->glfw_window(), GLFW_KEY_S) == GLFW_PRESS) {
            translation_matrix(2, 3) -= distance;
        }
        if (glfwGetKey(screen()->glfw_window(), GLFW_KEY_D) == GLFW_PRESS) {
            translation_matrix(0, 3) -= distance;
        }
        if (glfwGetKey(screen()->glfw_window(), GLFW_KEY_SPACE) == GLFW_PRESS) {
            translation_matrix(1, 3) -= distance;
        }
        if (glfwGetKey(screen()->glfw_window(), GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) { // GLFW_KEY_GRAVE_ACCENT / GLFW_KEY_LEFT_SHIFT
            translation_matrix(1, 3) += distance;
        }
        transformation_matrix_ = translation_matrix * transformation_matrix_;
    }

    size_t n_indices = indices_.size();
    shader_->set_buffer("indices", nanogui::VariableType::UInt32, {n_indices}, &indices_[0]);
    shader_->set_buffer("aPosition", nanogui::VariableType::Float32, {n_indices, 3}, &positions_[0]);
    shader_->set_buffer("aNormal", nanogui::VariableType::Float32, {n_indices, 3}, &normals_[0]);
    shader_->set_buffer("aColor", nanogui::VariableType::Float32, {n_indices, 3}, &colors_[0]);

    glDisable(GL_CULL_FACE);
    Eigen::Matrix4d look_at_matrix = Eigen::Matrix4d::Identity();
    look_at_matrix(2, 3) = -scale_;
    Eigen::Matrix4d model_view_matrix = look_at_matrix * transformation_matrix_;
    Eigen::Matrix4d normal_matrix = model_view_matrix.inverse().transpose();
    shader_->set_uniform("uModelViewMatrix", to_nanogui(model_view_matrix));
    shader_->set_uniform("uNormalMatrix", to_nanogui(normal_matrix));

    shader_->begin();
    shader_->draw_array(nanogui::Shader::PrimitiveType::Triangle, 0, n_indices, true);
    shader_->end();

    last_frame_ = curr_frame;
}

bool Canvas::mouse_button_event(const nanogui::Vector2i& p, int button, bool down, int modifiers) {
    if (down) {
        request_focus();
    }
    if (select_ && button == GLFW_MOUSE_BUTTON_1 && down) {
        // Ray
        Eigen::Matrix4d look_at_matrix = Eigen::Matrix4d::Identity();
        look_at_matrix(2, 3) = -scale_;
        Eigen::Matrix4d model_view_matrix_inverse = (look_at_matrix * transformation_matrix_).inverse();
        Eigen::Vector3d a = (model_view_matrix_inverse * Eigen::Vector4d{0, 0, 0, 1}).hnormalized();
        int x = p[0] - SIZE / 2;
        if (parameter_space_) {
            x -= SIZE;
        }
        int y = SIZE / 2 - p[1];
        Eigen::Vector3d b = (model_view_matrix_inverse * Eigen::Vector4d{x, y, -SIZE / 2, 1}).hnormalized();
        Eigen::Vector3d r = (b - a).normalized();

        // Cast
        double min = std::numeric_limits<double>::max();
        int idx = -1;
        for (int i = 0; i < indices_.size() / 3; i++) {
            std::vector<Eigen::Vector3d> triangle(3);
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    triangle[j](k) = positions_[3 * indices_[3 * i + j] + k];
                }
            }
            Eigen::Vector3d n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]).normalized();
            double nTr = n.dot(r);
            if (nTr == 0) {
                continue;
            }
            double lambda = (n.dot(triangle[0]) - n.dot(a)) / nTr;
            if (lambda < 0) {
                continue;
            }
            if (lambda < min) {
                Eigen::Vector3d s = a + lambda * r;
                for (int j = 0; j < 3; j++) {
                    triangle[j] -= s;
                }
                std::vector<Eigen::Vector3d> normals(3);
                for (int j = 0; j < 3; j++) {
                    normals[j] = triangle[j].cross(triangle[(j + 1) % 3]);
                }
                if (normals[0].dot(normals[1]) >= 0 && normals[0].dot(normals[2]) >= 0) {
                    min = lambda;
                    idx = i;
                }
            }
        }

        // Select
        if (idx != -1) {
            switch (sources_[idx].first) {
            case Element::CELL: {
                OpenVolumeMesh::CellHandle c(sources_[idx].second);
                selected_c_[c] = !selected_c_[c];
                if (selected_c_[c]) {
                    std::cout << "Cell " << c.idx() << std::endl;
                }
                break;
            }
            case Element::FACE: {
                OpenVolumeMesh::FaceHandle f(sources_[idx].second);
                selected_f_[f] = !selected_f_[f];
                if (selected_f_[f]) {
                    std::cout << "Face " << f.idx() << std::endl;
                }
                break;
            }
            case Element::EDGE: {
                OpenVolumeMesh::EdgeHandle e(sources_[idx].second);
                selected_e_[e] = !selected_e_[e];
                if (selected_e_[e]) {
                    std::cout << "Edge " << e.idx() << std::endl;
                }
                break;
            }
            case Element::VERTEX: {
                OpenVolumeMesh::VertexHandle v(sources_[idx].second);
                selected_v_[v] = !selected_v_[v];
                if (selected_v_[v]) {
                    std::cout << "Vertex " << v.idx() << std::endl;
                }
                break;
            }
            default:
                std::cout << "Other " << sources_[idx].second << std::endl;
                break;
            }
        }
        update();
        other_->update();
        return true;
    }
    if (button == GLFW_MOUSE_BUTTON_1 && !rot_fps_ && !trans_) {
        if (down) {
            rot_ = true;
            x0_ = p[0] - SIZE / 2;
            if (parameter_space_) {
                x0_ -= SIZE;
            }
            y0_ = SIZE / 2 - p[1];
        } else {
            rot_ = false;
        }
        return true;
    }
    if (button == GLFW_MOUSE_BUTTON_2 && !rot_ && !trans_) {
        if (down) {
            rot_fps_ = true;
            glfwSetInputMode(screen()->glfw_window(), GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        } else {
            rot_fps_ = false;
            glfwSetInputMode(screen()->glfw_window(), GLFW_CURSOR, GLFW_CURSOR_NORMAL);
        }
        return true;
    }
    if (button == GLFW_MOUSE_BUTTON_3 && !rot_ && !rot_fps_) {
        if (down) {
            trans_ = true;
            x0_ = p[0] - SIZE / 2;
            if (parameter_space_) {
                x0_ -= SIZE;
            }
            y0_ = SIZE / 2 - p[1];
        } else {
            trans_ = false;
        }
        return true;
    }
    return false;
}

bool Canvas::mouse_drag_event(const nanogui::Vector2i& p, const nanogui::Vector2i& rel, int button, int modifiers) {
    if (rot_) {
        int x1 = p[0] - SIZE / 2;
        if (parameter_space_) {
            x1 -= SIZE;
        }
        int y1 = SIZE / 2 - p[1];
        if (x1 != x0_ || y1 != y0_) {
            Eigen::Vector3d v0 = {x0_, y0_, SIZE / 2};
            v0.normalize();
            Eigen::Vector3d v1 = {x1, y1, SIZE / 2};
            v1.normalize();
            Eigen::Vector3d x = v0.cross(v1);
            double sin = x.norm();
            double cos = v0.dot(v1);
            Eigen::Vector3d n = x / sin;
            Eigen::Matrix3d xn;
            xn << 0, -n(2), n(1),
                n(2), 0, -n(0),
                -n(1), n(0), 0;
            Eigen::Matrix4d rotation_matrix = Eigen::Matrix4d::Identity();
            rotation_matrix.block<3, 3>(0, 0) = cos * Eigen::Matrix3d::Identity() + (1 - cos) * n * n.transpose() + sin * xn;
            transformation_matrix_ = rotation_matrix * transformation_matrix_;
            x0_ = x1;
            y0_ = y1;
        }
        return true;
    }
    if (rot_fps_) {
        if (rel[0] != 0 || rel[1] != 0) {
            Eigen::Vector3d v0 = {0, 0, 1};
            Eigen::Vector3d v1 = {rel[0], -rel[1], SIZE / 2};
            v1.normalize();
            Eigen::Vector3d x = v0.cross(v1);
            double sin = x.norm();
            double cos = v0.dot(v1);
            Eigen::Vector3d n = x / sin;
            Eigen::Matrix3d xn;
            xn << 0, -n(2), n(1),
                n(2), 0, -n(0),
                -n(1), n(0), 0;
            Eigen::Matrix4d rotation_matrix = Eigen::Matrix4d::Identity();
            rotation_matrix.block<3, 3>(0, 0) = cos * Eigen::Matrix3d::Identity() + (1 - cos) * n * n.transpose() + sin * xn;
            Eigen::Matrix4d look_at_matrix = Eigen::Matrix4d::Identity();
            look_at_matrix(2, 3) = -scale_;
            transformation_matrix_ = look_at_matrix.inverse() * rotation_matrix * look_at_matrix * transformation_matrix_;
        }
        return true;
    }
    return false;
}

bool Canvas::mouse_motion_event(const nanogui::Vector2i& p, const nanogui::Vector2i& rel, int button, int modifiers) {
    if (trans_) {
        int x1 = p[0] - SIZE / 2;
        if (parameter_space_) {
            x1 -= SIZE;
        }
        int y1 = SIZE / 2 - p[1];
        if (x1 != x0_ || y1 != y0_) {
            Eigen::Vector3d t = {(x1 - x0_) * 2 * scale_ / SIZE, (y1 - y0_) * 2 * scale_ / SIZE, 0};
            Eigen::Matrix4d translation_matrix = Eigen::Matrix4d::Identity();
            translation_matrix.block<3, 1>(0, 3) = t;
            transformation_matrix_ = translation_matrix * transformation_matrix_;
            x0_ = x1;
            y0_ = y1;
        }
        return true;
    }
    return false;
}

bool Canvas::scroll_event(const nanogui::Vector2i& p, const nanogui::Vector2f& rel) {
    if (glfwGetMouseButton(screen()->glfw_window(), GLFW_MOUSE_BUTTON_1) == GLFW_RELEASE &&
        glfwGetMouseButton(screen()->glfw_window(), GLFW_MOUSE_BUTTON_2) == GLFW_RELEASE &&
        glfwGetMouseButton(screen()->glfw_window(), GLFW_MOUSE_BUTTON_3) == GLFW_RELEASE) {
        request_focus();
        if (glfwGetKey(screen()->glfw_window(), GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS && focused()) {
            if (rel[0] > 0 || rel[1] > 0) {
                size_ *= 1.148698355; // = 2^(1/5) -> 5 scrolls from 1 to 0.5
                if (size_ > 1) {
                    size_ = 1;
                }
            } else {
                size_ /= 1.148698355;
            }
            update();
        } else {
            Eigen::Matrix4d scaling_matrix = Eigen::Matrix4d::Identity();
            if (rel[0] > 0 || rel[1] > 0) {
                scaling_matrix.block<3, 3>(0, 0) *= 1.1;
            } else {
                scaling_matrix.block<3, 3>(0, 0) /= 1.1;
            }
            transformation_matrix_ = scaling_matrix * transformation_matrix_;
        }
        return true;
    }
    return false;
}

bool Canvas::keyboard_event(int key, int scancode, int action, int modifiers) {
    if (key == GLFW_KEY_LEFT_ALT && action == GLFW_PRESS) {
        select_ = true;
        return true;
    }
    if (key == GLFW_KEY_LEFT_ALT && action == GLFW_RELEASE) {
        select_ = false;
        return true;
    }
    if (focused()) {
        if (key == GLFW_KEY_C && action == GLFW_PRESS) {
            show_cells_ ^= true;
            update();
            return true;
        }
        if (key == GLFW_KEY_F && action == GLFW_PRESS) {
            show_faces_ ^= true;
            update();
            return true;
        }
        if (key == GLFW_KEY_E && action == GLFW_PRESS) {
            show_edges_ ^= true;
            update();
            return true;
        }
        if (key == GLFW_KEY_V && action == GLFW_PRESS) {
            show_vertices_ ^= true;
            update();
            return true;
        }
        if (key == GLFW_KEY_B && action == GLFW_PRESS) {
            show_boundary_ ^= true;
            update();
            return true;
        }
        if (key == GLFW_KEY_X && action == GLFW_PRESS) {
            show_coords_ ^= true;
            update();
            return true;
        }
        if (key == GLFW_KEY_M && action == GLFW_PRESS) {
            show_colored_ ^= true;
            update();
            return true;
        }
        if (key == GLFW_KEY_K && action == GLFW_PRESS) {
            white_background_ ^= true;
            update();
            return true;
        }
        if (key == GLFW_KEY_R && action == GLFW_PRESS) {
            init();
            return true;
        }
        if (key == GLFW_KEY_U && action == GLFW_PRESS) {
            for (auto c : mesh_.cells()) {
                selected_c_[c] = false;
            }
            for (auto f : mesh_.faces()) {
                selected_f_[f] = false;
            }
            for (auto e : mesh_.edges()) {
                selected_e_[e] = false;
            }
            for (auto v : mesh_.vertices()) {
                selected_v_[v] = false;
            }
            update();
            other_->update();
            return true;
        }
    }
    return false;
}

std::mutex& Canvas::get_mutex() {
    return mutex_;
}
