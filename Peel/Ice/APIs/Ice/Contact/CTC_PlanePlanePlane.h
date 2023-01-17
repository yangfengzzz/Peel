///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for intersection of 3 planes
 *	\file		CTC_PlanePlanePlane.h
 *	\author		Pierre Terdiman
 *	\date		January, 13, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef CTCPLANEPLANEPLANE_H
#define CTCPLANEPLANEPLANE_H

CONTACT_API bool ComputeIntersection(const Plane& p1, const Plane& p2, const Plane& p3, Point& v);

#endif  // CTCPLANEPLANEPLANE_H
