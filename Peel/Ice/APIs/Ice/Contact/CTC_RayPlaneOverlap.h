///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for ray-plane intersection
 *	\file		CTC_RayPlaneOverlap.h
 *	\author		Pierre Terdiman
 *	\date		January, 13, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef CTCRAYPLANEOVERLAP_H
#define CTCRAYPLANEOVERLAP_H

// Ray-Plane intersection
CONTACT_API Point RayPlane(const Point& v1, const Point& v2, const Plane& plane, float& d);
CONTACT_API bool RayPlane(const Ray& line, const Plane& plane, float& distance_along_line, Point& point_on_plane);

#endif  // CTCRAYPLANEOVERLAP_H
