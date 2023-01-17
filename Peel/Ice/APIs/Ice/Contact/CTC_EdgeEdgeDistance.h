///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for edge-edge distance
 *	\file		CTC_EdgeEdgeDistance.h
 *	\author		Pierre Terdiman (original code from Eric Larsen in PQP)
 *	\date		January, 13, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef CTCEDGEEDGEDISTANCE_H
#define CTCEDGEEDGEDISTANCE_H

CONTACT_API void EdgeEdgeDist(
        Point& x, Point& y, const Point& p0, const Point& dir0, const Point& p1, const Point& dir1);

#endif  // CTCEDGEEDGEDISTANCE_H
