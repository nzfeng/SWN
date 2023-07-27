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

// ===================== NUMERICAL

int mod(int a, int b); // Pythonic modulus

int roundToNearestInteger(double x);

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

// write output scalar functions; also functions for outputting curves (curve completions + completion of nonbounding
// loops)

/* Input curves are already assumed to have been organized into connected components. */
void exportCurvesAsOBJ(const VertexData<Vector3>& vertexPositions, const std::vector<SurfacePoint>& curveNodes,
                       const std::vector<std::vector<std::array<size_t, 2>>>& curveEdges, const std::string& filename);

// ===================== CURVE MANIPULATION

Halfedge determineHalfedgeFromVertices(const Vertex& vA, const Vertex& vB);

SurfacePoint reinterpretTo(const SurfacePoint& p, SurfaceMesh& otherMesh);

std::vector<Halfedge> convertToHalfedges(const std::vector<SurfacePoint>& curveNodes,
                                         const std::vector<std::array<size_t, 2>>& curveEdges);

Vector<int> convertToChain(IntrinsicGeometryInterface& geom, const std::vector<Halfedge>& curve);

std::vector<std::vector<Halfedge>> getCurveComponents(IntrinsicGeometryInterface& geom,
                                                      const std::vector<Halfedge>& curveHalfedges);

std::vector<std::vector<std::array<size_t, 2>>>
getCurveComponents(SurfaceMesh& mesh, const std::vector<SurfacePoint>& curveNodes,
                   const std::vector<std::array<size_t, 2>>& curveEdges);


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