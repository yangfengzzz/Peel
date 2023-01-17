///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for swept ellipsoid intersection
 *	\file		CTC_SweptEllipsoid.h
 *	\author		Pierre Terdiman
 *	\date		January, 13, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef CTCSWEPTELLIPSOID_H
#define CTCSWEPTELLIPSOID_H

CONTACT_API bool LSSTriangleOverlap(const Triangle& triangle,
                                    const Point& p0,
                                    const Point& dir,
                                    float radius,
                                    Point* hit_point = null,
                                    float* t = null);

#ifdef OLDIES
CONTACT_API bool LSETriangleOverlap(const Triangle& triangle,
                                    const Point& p0,
                                    const Point& dir,
                                    const Point& axis0,
                                    const Point& axis1,
                                    const Point& axis2,
                                    Point* hit_point = null,
                                    float* t = null);
#endif

#endif  // CTCSWEPTELLIPSOID_H
