#include "util.h"

Mesh open(std::string filename, int& orientation) {
    Mesh mesh;
    std::function<void(std::vector<OpenVolumeMesh::VertexHandle>)> set_orientation = [&](std::vector<OpenVolumeMesh::VertexHandle> tet) {
        Matrix3q m;
        for (int i = 0; i < 3; i++) {
            m.col(i) = mesh.vertex(tet[i + 1]).cast<mpq_class>() - mesh.vertex(tet[0]).cast<mpq_class>();
        }
        orientation = m.determinant() > 0 ? 1 : -1;
    };
    std::vector<std::vector<std::string>> input;
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> split_line;
        std::istringstream ss(line);
        std::string s;
        while (std::getline(ss, s, ' ')) {
            split_line.push_back(s);
        }
        input.push_back(split_line);
    }
    int vertices_idx = -1;
    int n_vertices = -1;
    int cells_idx = -1;
    int n_cells = -1;
    for (int i = 0; i < input.size(); i++) {
        if (input[i].size() > 0 && input[i][0] == "POINTS") {
            vertices_idx = i + 1;
            n_vertices = std::stoi(input[i][1]);
        }
        if (input[i].size() > 0 && input[i][0] == "CELLS") {
            cells_idx = i + 1;
            n_cells = std::stoi(input[i][1]);
        }
    }
    std::vector<OpenVolumeMesh::VertexHandle> vertices;
    if (input[vertices_idx].size() == 3) {
        for (int i = vertices_idx; i < vertices_idx + n_vertices; i++) {
            vertices.push_back(mesh.add_vertex({std::stod(input[i][0]),
                                                std::stod(input[i][1]),
                                                std::stod(input[i][2])}));
        }
    } else {
        for (int i = 0; i < n_vertices; i++) {
            vertices.push_back(mesh.add_vertex({std::stod(input[vertices_idx][3 * i]),
                                                std::stod(input[vertices_idx][3 * i + 1]),
                                                std::stod(input[vertices_idx][3 * i + 2])}));
        }
    }
    std::vector<std::vector<OpenVolumeMesh::VertexHandle>> tets;
    if (input[cells_idx][0] == "OFFSETS") {
        n_cells--;
        int offset = cells_idx + n_cells + 3;
        for (int i = 0; i < n_cells; i++) {
            tets.push_back({vertices[std::stoi(input[offset + 4 * i][0])],
                            vertices[std::stoi(input[offset + 4 * i + 1][0])],
                            vertices[std::stoi(input[offset + 4 * i + 2][0])],
                            vertices[std::stoi(input[offset + 4 * i + 3][0])]});
        }
    } else {
        if (input[cells_idx].size() == 5) {
            for (int i = cells_idx; i < cells_idx + n_cells; i++) {
                tets.push_back({vertices[std::stoi(input[i][1])],
                                vertices[std::stoi(input[i][2])],
                                vertices[std::stoi(input[i][3])],
                                vertices[std::stoi(input[i][4])]});
            }
        } else {
            for (int i = 0; i < n_cells; i++) {
                tets.push_back({vertices[std::stoi(input[cells_idx + 5 * i + 1][0])],
                                vertices[std::stoi(input[cells_idx + 5 * i + 2][0])],
                                vertices[std::stoi(input[cells_idx + 5 * i + 3][0])],
                                vertices[std::stoi(input[cells_idx + 5 * i + 4][0])]});
            }
        }
    }
    if (orientation == 0) {
        set_orientation(tets[0]);
    }
    for (int i = 0; i < tets.size(); i++) {
        mesh.add_cell(tets[i][0],
                      tets[i][1],
                      tets[i][orientation > 0 ? 2 : 3],
                      tets[i][orientation > 0 ? 3 : 2]);
    }
    return mesh;
}

bool tutte3d(Mesh& mesh) {
    if (!mesh.vertex_property_exists<Eigen::Vector3d>("parameter")) {
        return false;
    }
    OpenVolumeMesh::VertexPropertyT<Eigen::Vector3d> parameter = mesh.request_vertex_property<Eigen::Vector3d>("parameter");
    OpenVolumeMesh::VertexPropertyT<int> idx = mesh.request_vertex_property<int>("idx");
    int k = 0;
    for (auto v : mesh.vertices()) {
        if (!mesh.is_boundary(v)) {
            idx[v] = k;
            k++;
        }
    }
    if (k > 0) {
        std::vector<Eigen::Triplet<double>> trips;
        Eigen::MatrixXd b = Eigen::MatrixXd::Zero(k, 3);
        for (auto v : mesh.vertices()) {
            if (!mesh.is_boundary(v)) {
                int n_neighbors = 0;
                for (auto vv : mesh.vertex_vertices(v)) {
                    if (mesh.is_boundary(vv)) {
                        b.row(idx[v]) += parameter[vv];
                    } else {
                        trips.push_back(Eigen::Triplet<double>(idx[v], idx[vv], -1));
                    }
                    n_neighbors++;
                }
                trips.push_back(Eigen::Triplet<double>(idx[v], idx[v], n_neighbors));
            }
        }
        Eigen::SparseMatrix<double> a(k, k);
        a.setFromTriplets(trips.begin(), trips.end());
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(a);
        if (solver.info() != Eigen::Success) {
            return false;
        }
        Eigen::MatrixXd x = solver.solve(b);
        for (auto v : mesh.vertices()) {
            if (!mesh.is_boundary(v)) {
                parameter[v] = x.row(idx[v]);
            }
        }
    }
    return true;
}

bool parameter_from_mesh(Mesh& mesh, Mesh& parameter_mesh) {
    if (mesh.n_vertices() != parameter_mesh.n_vertices()) {
        return false;
    }
    OpenVolumeMesh::VertexPropertyT<Eigen::Vector3d> parameter = mesh.request_vertex_property<Eigen::Vector3d>("parameter");
    OpenVolumeMesh::VertexPropertyT<Eigen::Vector3d> surface_parameter = mesh.request_vertex_property<Eigen::Vector3d>("surface_parameter");
    auto v_it = parameter_mesh.vertices_begin();
    for (auto v : mesh.vertices()) {
        parameter[v] = parameter_mesh.vertex(*v_it);
        if (mesh.is_boundary(v)) {
            surface_parameter[v] = parameter[v];
        }
        v_it++;
    }
    mesh.set_persistent(parameter);
    mesh.set_persistent(surface_parameter);
    return true;
}

void save(Mesh& mesh, std::string object_name, std::string parameter_name, std::string handles_name) {
    OpenVolumeMesh::VertexPropertyT<Eigen::Vector3d> parameter = mesh.request_vertex_property<Eigen::Vector3d>("parameter");
    std::vector<std::string> filenames;
    if (object_name != "") {
        filenames.push_back(object_name);
    }
    if (parameter_name != "") {
        filenames.push_back(parameter_name);
    }
    for (std::string filename : filenames) {
        std::ofstream file;
        file.precision(std::numeric_limits<double>::max_digits10);
        file.open(filename);
        file << "# vtk DataFile Version 2.0" << std::endl
             << filename << std::endl
             << "ASCII" << std::endl
             << "DATASET UNSTRUCTURED_GRID" << std::endl
             << "POINTS " << mesh.n_vertices() << " double" << std::endl;
        for (auto v : mesh.vertices()) {
            Eigen::Vector3d p;
            if (filename == object_name) {
                p = mesh.vertex(v);
            } else {
                p = parameter[v];
            }
            for (int i = 0; i < 3; i++) {
                if (std::abs(p(i)) < 1e-100) {
                    p(i) = 0;
                }
            }
            file << p(0) << " " << p(1) << " " << p(2) << std::endl;
        }
        file << "CELLS " << mesh.n_cells() << " " << mesh.n_cells() * 5 << std::endl;
        for (auto c : mesh.cells()) {
            std::vector<int> tet;
            for (auto cv : mesh.tet_vertices(c)) {
                tet.push_back(cv.idx());
            }
            file << "4 " << tet[0] << " " << tet[1] << " " << tet[3] << " " << tet[2] << std::endl;
        }
        file << "CELL_TYPES " << mesh.n_cells() << std::endl;
        for (int i = 0; i < mesh.n_cells(); i++) {
            file << "10" << std::endl;
        }
        file.close();
    }
    if (handles_name != "") {
        std::ofstream handles_file;
        handles_file.open(handles_name);
        for (auto v : mesh.vertices()) {
            if (mesh.is_boundary(v)) {
                handles_file << v.idx() << std::endl;
            }
        }
        handles_file.close();
    }
}
