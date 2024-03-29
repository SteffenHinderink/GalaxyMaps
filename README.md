# Galaxy Maps

<img src=galaxy-maps-image.jpg width=375 align=right>

This is an implementation of the volumetric mapping method presented in the SIGGRAPH 2023 paper
"[**Galaxy Maps: Localized Foliations for Bijective Volumetric Mapping**](http://graphics.cs.uos.de/papers/galaxy-maps-preprint.pdf)".

Given a tetrahedral mesh (ball-topology) and prescribed 3D mapping coordinates for its boundary vertices (convex or, more generally, star-shaped),
it addresses the problem of computing mapping coordinates also for the interior vertices,
such that together they define a **bijective** piecewise linear map.

Previous methods, while sometimes showing very high success rates, do not guarantee a bijective result.
In case of a non-bijective result, some fraction of the tetrahedra are inverted or degenerated by the map,
or even mapped outside of the prescribed boundary.

The Galaxy Maps approach takes such imperfect (non-bijective) maps as input and effectively turns them into bijective maps.
Defective spots are identified, enclosed in suitably shaped *stars*, and the contained part of the map exchanged for a bijective replacement.
These bijective local maps are reliably computed using a generalization of a reliable foliation-based mapping idea [Foliations].

If no such input map is readily available, the code includes the option to compute a simple 3D Tutte embedding [Tutte3D] as starting point.

*Note*: The output mesh is a refined version of the input mesh (exclusively splits; no general remeshing).
Note that some refinement can actually be inevitable to permit any piecewise linear map at all.
However, the code sometimes performs much more refinement than would strictly be necessary.

## Citation

If you find this code useful, please consider citing our paper:

```bibtex
@article{GalaxyMaps,
  author  = {Steffen Hinderink and Marcel Campen},
  title   = {Galaxy Maps: Localized Foliations for Bijective Volumetric Mapping},
  year    = {2023},
  journal = {ACM Trans. Graph.},
  volume  = {42},
  number  = {4},
  pages   = {129:1--129:16}
}
```

## Prerequisites

The code can be used as either a **library** or an **executable**.

In either case, the following libraries need to be installed:

- [CGAL](https://www.cgal.org)
- [Eigen](https://eigen.tuxfamily.org)
- [GMP](https://gmplib.org)

The following libraries are included in the [ext](ext) directory as submodules:

- [OpenVolumeMesh](https://gitlab.vci.rwth-aachen.de:9000/OpenVolumeMesh/OpenVolumeMesh.git)
- [NanoGUI](https://github.com/mitsuba-renderer/nanogui.git)

Remember to clone this repository recursively:\
```git clone --recursive https://github.com/SteffenHinderink/GalaxyMaps.git```

## Library

### CMake

The library can be included via CMake like this:

```cmake
add_subdirectory(path/to/GalaxyMaps)
target_link_libraries(yourTarget GalaxyMaps)
```

### Usage

The function ```galaxy_map``` implements the Galaxy Map method.
It parametrizes and refines a given mesh object.
[OpenVolumeMesh](https://www.graphics.rwth-aachen.de/software/openvolumemesh/)
is used to represent the tetrahedral meshes.
The parametrization is stored as the vertex property ```parameter``` of type ```Eigen::Vector3d``` in floating point
or ```qarameter``` of type ```Vector3q``` in rational numbers.
Both can be used as input for the ```galaxy_map``` function.
The 3D Tutte embedding is implemented in the function ```tutte3d``` using floating point numbers.

```cpp
#include <galaxy.h>

Mesh mesh = ...;
auto parameter = mesh.request_vertex_property<Eigen::Vector3d>("parameter");
for (auto v : mesh.vertices()) { // set boundary map (and possibly initial imperfect interior map)
    parameter[v] = ...;
}
mesh.set_persistent(parameter);

tutte3d(mesh); // optional, in case no initial mapping for the interior was set

galaxy_map(mesh);

auto rational_parameter = mesh.request_vertex_property<Vector3q>("qarameter");
for (auto v : mesh.vertices()) {
    std::cout << "(" << rational_parameter[v].transpose() << ")";
    std::cout << " ≈ (" << parameter[v].transpose() << ")" << std::endl;
}
```

The foliation-based mapping method can also be used independently (globally, without the galaxy approach):

```cpp
Mesh star = ...;
auto rational_parameter = star.request_vertex_property<Vector3q>("qarameter");
for (auto v : star.vertices()) {
    if (star.is_boundary(v)) { // set boundary map
        rational_parameter[v] = ...; // must be star-shaped with respect to origin (0, 0, 0)
    }
}
star.set_persistent(rational_parameter);

FoliationParameterization foliation_para(star);
foliation_para.set_surface_parameterization();
foliation_para.generate_directions();

// direct evaluation of the foliation-based map (without conversion to an exlicit PL map)
MeshPoint p = ...; // barycentric coordinates in a tetrahedron
Vector3q x = foliation_para.psi(p);

std::cout << "Psi(p) = (" << x.transpose() << ")" << std::endl;

MeshPoint q = foliation_para.psi_inverse(x); // = p

// convert to a PL map on a refined mesh
Mesh refined_star = foliation_para.refine();
decimate_refined(refined_star);

rational_parameter = refined_star.request_vertex_property<Vector3q>("qarameter");
for (auto v : refined_star.vertices()) {
    std::cout << "(" << rational_parameter[v].transpose() << ")" << std::endl;
}
```

## Executable

### Building

- ```mkdir build```
- ```cd build```
- ```cmake [-DGUI=OFF] ..``` (The option ```GUI``` controls if the program is built with a GUI)
- ```make -j```

### Usage

```./galaxy_map <mesh> <parameter_mesh> [-t] [-o <output>]```

- ```<mesh>```:
Tetrahedral mesh file in ```.vtk``` format to be parametrized.
- ```<parameter_mesh>```:
Tetrahedral mesh file in ```.vtk``` format, containing the initial (non-bijective) parametrization.
The vertices of the parameter mesh are assumed to correspond to the vertices of the input mesh with the same indices.
- ```-t```:
Option to use 3D Tutte embedding instead of the parameter mesh for the initial parametrization.
The parameter mesh is used for the surface parametrization.
- ```-o <output>```:
Optional argument to output the files ```<output>_object.vtk``` and ```<output>_parameter.vtk``` containing the refined mesh and its parametrization.
*Warning*: The truncation to floating point numbers required for this may introduce tiny inversions or degeneracies.

e.g. ```./galaxy_map ../meshes/example_object.vtk ../meshes/example_parameter.vtk```

## Examples

The [meshes](meshes) directory includes exemplary instances from different datasets:
- ```example```:
Small 3D Tutte counterexample based on [Tutte3D].
- ```tlc_polycube_bird```:
Mesh of a bird from the dataset from [TLC].
- ```tgs_100349```:
Triangle mesh of a cat from [Thingi10K] converted into a tetrahedral mesh with [TetGen] with the surface parametrized to a sphere that neither [TLC] nor [FFM] succeed to parametrize.
- ```cub_rockerarm291```:
Block of a motorcycle complex [MC3D] that neither [TLC] nor [FFM] succeed to parametrize.

## Dataset

The dataset of volumetric mapping problem instances that were created for benchmarking in the Galaxy Maps paper can be downloaded here:
[**Dataset**](http://graphics.cs.uos.de/papersdata/galaxy-maps-data.zip).

## References

- [FFM] Garanzha, V., Kaporin, I., Kudryavtseva, L., Protais, F., Ray, N., & Sokolov, D. (2021). Foldover-free maps in 50 lines of code.
- [Foliations] Campen, M., Silva, C. T., & Zorin, D. (2016). Bijective maps from simplicial foliations.
- [MC3D] Brückler, H., Gupta, O., Mandad, M., & Campen, M. (2022). The 3d motorcycle complex for structured volume decomposition.
- [TetGen] Hang, S. (2015). TetGen, a Delaunay-based quality tetrahedral mesh generator.
- [Thingi10K] Zhou, Q., & Jacobson, A. (2016). Thingi10K: A dataset of 10,000 3d-printing models.
- [TLC] Du, X., Aigerman, N., Zhou, Q., Kovalsky, S. Z., Yan, Y., Kaufman, D. M., & Ju, T. (2020). Lifting simplices to find injectivity.
- [Tutte3D] Floater, M. S., & Pham-Trong, V. (2006). Convex combination maps over triangulations, tilings, and tetrahedral meshes.
