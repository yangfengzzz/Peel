///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for ray-sphere intersection
 *	\file		CTC_RaySphereOverlap.h
 *	\author		Pierre Terdiman
 *	\date		January, 13, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef CTCRAYSPHEREOVERLAP_H
#define CTCRAYSPHEREOVERLAP_H

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Computes a segment-sphere intersection.
 *	Based on GD Mag code, but now works correctly when origin is inside the sphere.
 *	\param		origin		[in] ray origin
 *	\param		dir			[in] ray direction
 *	\param		length		[in] ray length
 *	\param		center		[in] sphere center
 *	\param		radius		[in] sphere radius
 *	\param		dist		[out] distance of intersection point
 *	\param		hit_pos		[out] intersection point
 *	\return		true if ray intersects sphere
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CONTACT_API bool SegmentSphere(const Point& origin,
                               const Point& dir,
                               float length,
                               const Point& center,
                               float radius,
                               float& dist,
                               Point& hit_pos);

inline_ bool SegmentSphere(
        const Point& origin, const Point& dir, float length, const Sphere& sphere, float& dist, Point& hit_pos) {
    return Ctc::SegmentSphere(origin, dir, length, sphere.mCenter, sphere.mRadius, dist, hit_pos);
}

CONTACT_API bool SegmentSphereOverlap(
        const Point& origin, const Point& dir, float length, const Point& center, float radius);
CONTACT_API bool RobustSegmentSphere(const Point& origin,
                                     const Point& dir,
                                     float length,
                                     const Point& center,
                                     float radius,
                                     float& dist,
                                     Point& hit_pos);

#endif  // CTCRAYSPHEREOVERLAP_H
