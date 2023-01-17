///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICEHULLFROMPLANES_H
#define ICEHULLFROMPLANES_H

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Creates hull vertices from a set of planes.
 *	\relates	ConvexHull
 *	\fn			CreateHullFromPlanes(Vertices& vertices, const Planes& planes)
 *	\param		vertices		[out] hull vertices
 *	\param		planes			[in] input planes
 *	\param		inside_point	[in] a point inside the volume defined by the set of planes
 *	\return		true if success
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FUNCTION MESHMERIZER_API bool CreateHullFromPlanes(Vertices& vertices,
                                                   const Planes& planes,
                                                   const Point* inside_point = null);

// ### test!
FUNCTION MESHMERIZER_API bool CreateHullFromPlanes2(Vertices& vertices,
                                                    const ConvexHull& hull,
                                                    const Point* inside_point = null);

#endif  // ICEHULLFROMPLANES_H
