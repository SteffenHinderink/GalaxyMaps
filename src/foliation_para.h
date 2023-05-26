#pragma once

#include "mesh.h"
#include "vectorq.h"

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/gmpxx.h>

#include <queue>

/**
 * Exception for the parameterization.
 * This is thrown when the parameterization cannot be evaluated
 * because there is no surface parameterization or foliation of the tetrahedral mesh.
 */
struct FoliationParameterizationException : public std::exception {
    /**
     * Returns the exception message.
     * @return exception message
     */
    const char* what() const noexcept override {
        return "Surface parameterization or direction field missing";
    }
};

/**
 * Struct for a point in a tetrahedral mesh in barycentric coordinates in a tetrahedron.
 */
struct MeshPoint {
    /// Tetrahedron in which the point lies and for which the barycentric coordinates are defined.
    OpenVolumeMesh::CellHandle c;
    /// Barycentric coordinates in the order of the vertices of the tetrahedron.
    std::vector<mpq_class> bary;
};

/**
 * Class for the parameterization of tetrahedral meshes through foliations.
 */
class FoliationParameterization {
public:
    /**
     * Constructor of the Para class.
     * @param mesh tetrahedral mesh that is parameterized.
     */
    FoliationParameterization(Mesh& mesh);

    /**
     * Sets the surface parameterization of the mesh through a vertex property.
     * The surface parameterization must map to a star-shaped polyhedron
     * so that it can be used for the tetrahedral mesh parameterization.
     * @param property name of the vertex property containing the surface parameterization
     * @return true if the surface parameterization maps to a star-shaped polyhedron, false otherwise
     */
    bool set_surface_parameterization(std::string property = "qarameter");

    /**
     * Generates the foliation directions.
     * @return true if the foliation direction generation was successful, false otherwise
     */
    bool generate_directions();

    /**
     * Evaluates the parameterization.
     * This requires a surface parameterization and the foliation.
     * @param p point for which the parameterization is evaluated
     * @return point to which p is mapped
     * @throws ParaException if the parameterization cannot be evaluated
     */
    Vector3q psi(MeshPoint p);

    /**
     * Evaluates the inverse of the parameterization.
     * This requires a surface parameterization and the foliation.
     * @param p point for which the inverse is evaluated
     * @return point that maps to q
     * @throws ParaException if the inverse cannot be evaluated
     */
    MeshPoint psi_inverse(Vector3q q);

    /**
     * Refines the mesh so that the parameterization can be evaluated through linear interpolation.
     * The parameterization is saved in the vertex property "psi" as vectors of rational numbers in the refined mesh.
     * The parameterization evaluation requires a surface parameterization and the foliation.
     * @param verbose flag that controls if the refinement process is logged
     * @return refined tetrahedral mesh
     * @throws ParaException if the parameterization cannot be evaluated
     */
    Mesh refine(bool verbose = false);

private:
    bool is_center(OpenVolumeMesh::CellHandle c);

    int& direction_coordinate(OpenVolumeMesh::CellHandle c, OpenVolumeMesh::VertexHandle v);

    mpq_class direction_length(OpenVolumeMesh::CellHandle c);

    void align_directions();

    /**
     * Projects a point along the foliation directions onto the boundary of the tetrahedron in which the point lies.
     * @param p point that is projected
     * @param in direction of the projection
     * @return projected point and length of the projection
     */
    std::pair<std::vector<mpq_class>, mpq_class> project(MeshPoint p, bool in = true);

    /**
     * Converts a point in a tetrahedron into the next tetrahedron along the foliation.
     * @param p point that is converted
     * @param in direction of the conversion
     * @return p in the next tetrahedron
     */
    MeshPoint next(MeshPoint p, bool in = true);

    /// Tetrahedral mesh that is parameterized.
    Mesh& mesh_;

    /// Flag that indicates if the surface parameterization has been set.
    bool surface_parameterized_;

    /// Flag that indicates if the foliation directions have been generated.
    bool shelled_;

    /// Tetrahedron property containing the foliation directions.
    OpenVolumeMesh::CellPropertyT<std::vector<int>> d_;

    OpenVolumeMesh::CellPropertyT<mpq_class> h_;

    /// Vertex property containing the surface parameterization.
    OpenVolumeMesh::VertexPropertyT<Vector3q> surface_parameter_;
};
