#pragma once

#include "mesh.h"
#include "vectorq.h"

#include <Eigen/Sparse>

#include <fstream>

Mesh open(std::string filename, int& orientation);

bool tutte3d(Mesh& mesh);

bool parameter_from_mesh(Mesh& mesh, Mesh& parameter_mesh);

void save(Mesh& mesh, std::string object_name, std::string parameter_name, std::string handles_name = "");
