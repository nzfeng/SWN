#include "utils.h"

// ===================== NUMERICAL

int mod(int a, int b) {
    return (b + (a % b)) % b;
}

/*
 * Round to nearest int. Add +- 0.5 because compiler always truncates.
 */
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


// ===================== MESH MUTATION


// ===================== VISUALIZATION