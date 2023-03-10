///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for segment-triangle distance
 *	\file		CTC_SegmentTriangleDistance.h
 *	\author		Pierre Terdiman
 *	\date		January, 13, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef CTCSEGMENTTRIANGLEDISTANCE_H
#define CTCSEGMENTTRIANGLEDISTANCE_H

CONTACT_API float SegmentTriangleSqrDist(const Segment& segment,
                                         const Point& p0,
                                         const Point& p1,
                                         const Point& p2,
                                         float* t = null,
                                         float* u = null,
                                         float* v = null);

inline_ float SegmentTriangleSqrDist(const Segment& segment,
                                     const IndexedTriangle& triangle,
                                     const Point* verts,
                                     float* t = null,
                                     float* u = null,
                                     float* v = null) {
    return Ctc::SegmentTriangleSqrDist(segment, verts[triangle.mRef[0]], verts[triangle.mRef[1]],
                                       verts[triangle.mRef[2]], t, u, v);
}

#endif  // CTCSEGMENTTRIANGLEDISTANCE_H
