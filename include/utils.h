#pragma once

#include "geometrycentral/surface/common_subdivision.h"
#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/intrinsic_triangulation.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/mutation_manager.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

// ===================== NUMERICAL

int mod(int a, int b); // Pythonic modulus

int roundToNearestInteger(double x);

// ===================== I/O

Halfedge determineHalfedgeFromVertices(const Vertex& vA, const Vertex& vB);

// write output scalar functions; also functions for outputting curves (curve completions + completion of nonbounding
// loops)

// ===================== MESH MUTATION

// SurfacePoint reinterpretTo(const SurfacePoint& p, SurfaceMesh& otherMesh);

// std::vector<Halfedge> remeshIfNecessary(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry,
//                                         const std::vector<Halfedge>& curveHalfedges,
//                                         const std::vector<SurfacePoint>& curveNodes,
//                                         const std::vector<std::array<size_t, 2>>& curveEdges);

// std::vector<Halfedge> remeshToConforming(ManifoldSurfaceMesh& manifoldMesh, VertexPositionGeometry& manifoldGeom,
//                                          const std::vector<SurfacePoint>& curveNodes,
//                                          const std::vector<std::array<size_t, 2>>& curveEdges);

// std::vector<Halfedge> splitEdgesIfNecessary(ManifoldSurfaceMesh& manifoldMesh, VertexPositionGeometry& manifoldGeom,
//                                             const std::vector<Halfedge>& curveHalfedges,
//                                             std::vector<size_t>& markedEdges, int maxEdgeSplits = 3);

// void resetIntrinsicTriangulationAndMarkEdges(std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation>& intTri,
//                                              ManifoldSurfaceMesh& manifoldMesh, VertexPositionGeometry& manifoldGeom,
//                                              const std::vector<Halfedge>& curveHalfedgesOnManifold);

// std::vector<Halfedge> determineHalfedgesInIntrinsicTriangulation(IntegerCoordinatesIntrinsicTriangulation& intTri,
//                                                                  const std::vector<Halfedge>&
//                                                                  curveHalfedgesOnManifold, std::vector<size_t>&
//                                                                  markedEdges);

// std::vector<Halfedge> determineHalfedgesInIntrinsicTriangulation(IntegerCoordinatesIntrinsicTriangulation& intTri,
//                                                                  const std::vector<Halfedge>&
//                                                                  curveHalfedgesOnManifold);

// ===================== VISUALIZATION