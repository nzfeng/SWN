#include "utils.h"

// ===================== NUMERICAL

int mod(int a, int b) {
    return (b + (a % b)) % b;
}

/* Round to nearest int. Add +- 0.5 because compiler always truncates. */
int roundToNearestInteger(double x) {

    double y = x + 0.5 - (x < 0);
    return (int)y;
}


// ===================== I/O

Halfedge determineHalfedgeFromVertices(const Vertex& vA, const Vertex& vB) {

    for (Halfedge he : vA.outgoingHalfedges()) {
        if (he.tipVertex() == vB) return he;
    }
    return Halfedge();
}

/* Read in curve, encoded as OBJ line objects. */
void readLines(const SurfaceMesh& mesh, const std::string& filepath, std::vector<SurfacePoint>& curveNodes,
               std::vector<std::array<size_t, 2>>& curveEdges, int offset) {

    std::ifstream curr_file(filepath.c_str());
    std::string line;
    std::string X; // re-written variable for holding the leading (indicator) char
    size_t idx;    // re-written variable for holding element index

    if (curr_file.is_open()) {

        while (!curr_file.eof()) {

            getline(curr_file, line);

            // Ignore any newlines
            if (line == "") {
                continue;
            }

            std::istringstream iss(line);
            iss >> X; // get the first char in the line

            if (X == "l") {
                std::vector<SurfacePoint> nodes;
                while (true) {
                    if (iss.eof()) break;
                    iss >> idx;
                    idx += offset; // OBJ elements are 1-indexed, whereas geometry-central uses 0-indexing
                    nodes.emplace_back(mesh.vertex(idx));
                }
                size_t N = nodes.size();
                size_t M = curveNodes.size();
                for (size_t i = 0; i < N - 1; i++) {
                    curveEdges.push_back({M + i, M + (i + 1)});
                }
                curveNodes.insert(curveNodes.back(), nodes.begin(), nodes.end());
            }
        }
        curr_file.close();
    } else {
        std::cerr << "Could not open file <" << filepath << ">" << std::endl;
    }
}

/* Read in curve data from a separate file, encoded as segments between barycentric points. */
void readCurves(const SurfaceMesh& mesh, const std::string& filepath, std::vector<SurfacePoint>& curveNodes,
                std::vector<std::array<size_t, 2>>& curveEdges, int offset) {

    std::ifstream curr_file(filepath.c_str());
    std::string line;
    std::string X; // re-written variable for holding the leading (indicator) char
    size_t idx;    // re-written variable for holding element index
    double a, b, c;

    if (curr_file.is_open()) {

        while (!curr_file.eof()) {

            getline(curr_file, line);

            // Ignore any newlines
            if (line == "") {
                continue;
            }

            std::istringstream iss(line);
            iss >> X; // get the first char in the line

            if (X == "v") {
                iss >> idx;
                idx += offset;
                curveNodes.emplace_back(mesh.vertex(idx));
            } else if (X == "e") {
                // WARNING: The edge index is assumed to agree with geometry-central's edge-indexing scheme. Depending
                // on which software you used to export curve data on the mesh, your edge indices may disagree with the
                // indices geometry-central computes from the vertices.
                iss >> idx >> a;
                idx += offset;
                curveNodes.emplace_back(mesh.edge(idx), a);
            } else if (X == "f") {
                iss >> idx >> a >> b >> c;
                idx += offset;
                curveNodes.emplace_back(mesh.face(idx), {a, b, c});
            } else if (X == "l") {
                std::vector<size_t> indices;
                while (true) {
                    if (iss.eof()) break;
                    iss >> idx;
                    idx += offset;
                    indices.push_back(idx);
                }
                size_t N = indices.size();
                for (size_t i = 0; i < N - 1; i++) {
                    curveEdges.push_back({indices[i], indices[i + 1]});
                }
            }
        }
        curr_file.close();
    } else {
        std::cerr << "Could not open file <" << filepath << ">" << std::endl;
    }
}


// ===================== MESH MUTATION


// ===================== VISUALIZATION

void displayCurves(const VertexPositionGeometry& geometry, const std::vector<SurfacePoint>& curveNodes,
                   const std::vector<std::array<size_t, 2>>& curveEdges) {

    std::vector<Vector3> nodes;
    for (const SurfacePoint& p : curveNodes) nodes.push_back(p.interpolate(geometry.vertexPositions));
    polyscope::registerCurveNetwork("input curve", nodes, curveEdges)->setColor({0, 0, 0})->setEnabled(true);
}