///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for triangle-rectangle distance
 *	\file		CTC_TriangleRectangleDistance.h
 *	\author		Pierre Terdiman
 *	\date		January, 13, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef CTCTRIANGLERECTANGLEDISTANCE_H
#define CTCTRIANGLERECTANGLEDISTANCE_H

// Triangle-rectangle squared distance
CONTACT_API float TriangleRectangleSqrDist(const Point& p0,
                                           const Point& p1,
                                           const Point& p2,
                                           const Rectangle3& rectangle,
                                           float* s = null,
                                           float* t = null,
                                           float* u = null,
                                           float* v = null);

inline_ float TriangleRectangleSqrDist(const IndexedTriangle& triangle,
                                       const Point* verts,
                                       const Rectangle3& rectangle,
                                       float* s = null,
                                       float* t = null,
                                       float* u = null,
                                       float* v = null) {
    return Ctc::TriangleRectangleSqrDist(verts[triangle.mRef[0]], verts[triangle.mRef[1]], verts[triangle.mRef[2]],
                                         rectangle, s, t, u, v);
}

#endif  // CTCTRIANGLERECTANGLEDISTANCE_H
