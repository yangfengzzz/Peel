///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for Catmull-Rom splines.
 *	\file		IceCatmullRom.h
 *	\author		Pierre Terdiman
 *	\date		October, 3, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICECATMULLROM_H
#define ICECATMULLROM_H

ICEMATHS_API void ComputeCatmullRom0(
        Point& q, const Point& p0, const Point& p1, const Point& p2, const Point& p3, float t);
ICEMATHS_API void ComputeCatmullRom1(
        Point& q, const Point& p0, const Point& p1, const Point& p2, const Point& p3, float t);

#endif  // ICECATMULLROM_H
