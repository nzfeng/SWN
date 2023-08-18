#pragma once

#include "geometrycentral/surface/common_subdivision.h"
#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/intrinsic_triangulation.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/mutation_manager.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/curve_network.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

inline double SPECIAL_VAL = std::numeric_limits<double>::quiet_NaN();

// ===================== NUMERICAL

template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

int mod(int a, int b); // Pythonic modulus

int roundToNearestInteger(double x);

double shift(const CornerData<double>& func, const std::vector<Halfedge>& curve);

FaceData<double> round(const CornerData<double>& f, const std::vector<Halfedge>& curve);

/* Get min and max of a vector, ignoring NaNs. */
std::tuple<double, double> minMax(const Vector<double>& vec);

// ===================== OPERATORS

template <typename T>
SparseMatrix<T> b1(SurfaceMesh& mesh);

template <typename T>
SparseMatrix<T> b2(SurfaceMesh& mesh);

CornerData<Vector2> toProjectiveCoordinates(const CornerData<double>& u);

CornerData<double> fromProjectiveCoordinates(const CornerData<Vector2>& u);

// ===================== I/O

/* Get home directory of a mesh file. */
std::string getHomeDirectory(const std::string& filepath);

bool isStringTrue(const std::string& input);

/* Read in curve, encoded either as OBJ line objects, or as dual edges. */
void readLines(SurfaceMesh& mesh, const std::string& filepath, std::vector<SurfacePoint>& curveNodes,
               std::vector<std::array<size_t, 2>>& curveEdges, std::vector<std::array<Face, 2>>& dualChain,
               int offset = -1);

/* Read in curve data (from a file NOT containing mesh data), encoded either as segments between barycentric points, or
 * as dual edges. */
void readCurves(SurfaceMesh& mesh, const std::string& filepath, std::vector<SurfacePoint>& curveNodes,
                std::vector<std::array<size_t, 2>>& curveEdges, std::vector<std::array<Face, 2>>& dualChain,
                int offset = -1);

/* Input curves are already assumed to have been organized into connected components. */
void exportCurves(const VertexData<Vector3>& vertexPositions, const std::vector<SurfacePoint>& curveNodes,
                  const std::vector<std::vector<std::array<size_t, 2>>>& curveEdges, const std::string& filename);

void exportCurves(const VertexData<Vector3>& vertexPositions, const std::vector<std::vector<Halfedge>>& curveHalfedges,
                  const std::string& filename);

void exportFunction(EmbeddedGeometryInterface& geom, const CornerData<double>& u, const std::string& filename);

/* Export on common subdivision. */
void exportFunction(IntegerCoordinatesIntrinsicTriangulation& intTri, VertexPositionGeometry& manifoldGeom,
                    const CornerData<double>& u, const std::string& filename);

// ===================== CURVE MANIPULATION

Halfedge determineHalfedgeFromVertices(const Vertex& vA, const Vertex& vB);

SurfacePoint reinterpretTo(const SurfacePoint& p, SurfaceMesh& otherMesh);

std::vector<Halfedge> convertToHalfedges(const std::vector<SurfacePoint>& curveNodes,
                                         const std::vector<std::array<size_t, 2>>& curveEdges);

Vector<int> convertToChain(const SurfaceMesh& mesh, const std::vector<Halfedge>& curve);

std::vector<std::vector<Halfedge>> getCurveComponents(IntrinsicGeometryInterface& geom,
                                                      const std::vector<Halfedge>& curveHalfedges,
                                                      bool useEndpoints = false);

std::vector<std::vector<std::array<size_t, 2>>>
getCurveComponents(SurfaceMesh& mesh, const std::vector<SurfacePoint>& curveNodes,
                   const std::vector<std::array<size_t, 2>>& curveEdges);

std::tuple<std::vector<Halfedge>, std::vector<Halfedge>>
getCurveDecomposition(const std::vector<Halfedge>& curveHalfedges, const CornerData<double>& vFunc,
                      double epsilon = 1e-1);

std::vector<Halfedge> getJumpLocus(const std::vector<Halfedge>& curveHalfedges, const FaceData<double>& func,
                                   double epsilon = 1e-1);

std::vector<Halfedge> getJumpLocus(const std::vector<Halfedge>& curveHalfedges, const CornerData<double>& func,
                                   double epsilon = 1e-1);

std::vector<Halfedge> getCompletedBoundingLoops(const std::vector<Halfedge>& curveHalfedges,
                                                const CornerData<double>& func, double epsilon = 1e-1);


// ===================== MESH MUTATION

void resetIntrinsicTriangulationAndMarkEdges(std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation>& intTri,
                                             ManifoldSurfaceMesh& manifoldMesh, VertexPositionGeometry& manifoldGeom,
                                             const std::vector<Halfedge>& curveHalfedgesOnManifold);

std::vector<Halfedge> determineHalfedgesInIntrinsicTriangulation(IntegerCoordinatesIntrinsicTriangulation& intTri,
                                                                 const std::vector<Halfedge>& curveHalfedgesOnManifold);

// ===================== VISUALIZATION

void displayCurves(const VertexPositionGeometry& geometry, const std::vector<Halfedge>& curveHalfedges,
                   const std::vector<SurfacePoint>& curveNodes, const std::vector<std::array<size_t, 2>>& curveEdges,
                   const std::vector<std::array<Face, 2>>& dualChain);

/*
 * Linearly interpolate CornerData defined on the intrinsic triangulation to the common subdivision.
 */
template <typename T>
CornerData<T> interpolateAcrossB(CommonSubdivision& cs, const CornerData<T>& dataB) {

    CornerData<T> interp(*cs.mesh);
    for (Vertex v : cs.mesh->vertices()) {
        SurfacePoint posB = cs.sourcePoints[v]->posB;
        if (posB.type == SurfacePointType::Vertex || posB.type == SurfacePointType::Edge) {
            for (Corner c : v.adjacentCorners()) {
                Face f = c.face();
                Face fB = cs.sourceFaceB[f];
                SurfacePoint pB_face = posB.inFace(fB);
                T val = 0.;
                size_t idx = 0;
                for (Halfedge he : fB.adjacentHalfedges()) {
                    Corner cB = he.corner();
                    val += dataB[cB] * pB_face.faceCoords[idx];
                    idx++;
                }
                interp[c] = val;
            }
        } else {
            // this vertex in the common subdivison lies within a face in the intrinsic mesh
            SurfacePoint pB_face = posB.inSomeFace();
            Face fB = pB_face.face;
            T val = 0.;
            size_t idx = 0;
            for (Halfedge he : fB.adjacentHalfedges()) {
                Corner cB = he.corner();
                val += dataB[cB] * pB_face.faceCoords[idx];
                idx++;
            }
            for (Corner c : v.adjacentCorners()) {
                interp[c] = val;
            }
        }
    }
    return interp;
}

CornerData<Vector2> interpolateVector2AcrossB(CommonSubdivision& cs, const CornerData<Vector2>& dataB);