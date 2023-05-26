#pragma once

#include "foliation_para.h"
#include "util.h"

#include <CGAL/QP_functions.h>
#include <CGAL/QP_models.h>

#include <random>

bool inverted(Mesh& mesh, OpenVolumeMesh::CellHandle c, bool parameter_space = true);

Eigen::Vector3d kernel_chebyshev_center(Mesh& mesh, std::vector<OpenVolumeMesh::HalfFaceHandle> polyhedron);

void decimate_refined(Mesh& mesh);

bool galaxy_map(Mesh& mesh);
