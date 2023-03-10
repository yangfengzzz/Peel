///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for hill-climbing (a.k.a. Local Search)
 *	\file		IceHillClimbing.h
 *	\author		Pierre Terdiman
 *	\date		January, 29, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICEHILLCLIMBING_H
#define ICEHILLCLIMBING_H

// Forward declarations
class IndexedSurface;
class Valencies;

//! This structure holds the cached information used in the separating vector algorithm.
struct MESHMERIZER_API SVCache : Pair {
    Point SV;           //!< Last separating vector
    udword pid;         //!< Last supporting vertex for polytope P
    udword qid;         //!< Last supporting vertex for polytope Q
    IndexedSurface* P;  //!< Polytope P
    IndexedSurface* Q;  //!< Polytope Q
                        //		BOOL				PState;	//!< Collision state
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Computes supporting vertex in a given direction. This is the brute-force O(n) version.
 *	\relates	CollisionHull
 *	\fn			ComputeSupportingVertex(udword nb_verts, const Point* verts, const Point& dir, const
 *Matrix4x4* mat) \param		nb_verts	[in] number of vertices \param		verts
 *[in] array of vertices \param		dir			[in] the direction \param		mat
 *[in] a possible matrix for vertices \return		udword	index of supporting vertex
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FUNCTION MESHMERIZER_API udword ComputeSupportingVertex(udword nb_verts,
                                                        const Point* verts,
                                                        const Point& dir,
                                                        const Matrix4x4* mat = null);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Performs local search within a list of vertices.
 *	\relates	CollisionHull
 *	\fn			LocalSearch(udword& id, const Point& dir, const Point* verts, const Valencies* val)
 *	\param		id			[in/out] initial vertex index / final supporting vertex index
 *	\param		dir			[in] separating vector
 *	\param		verts		[in] list of vertices
 *	\param		val			[in] valencies for the list of vertices
 *	\return		true if success
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FUNCTION MESHMERIZER_API bool LocalSearch(udword& id, const Point& dir, const Point* verts, const Valencies* val);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Performs local search within a list of vertices, using extra timestamps.
 *	\relates	CollisionHull
 *	\fn			LocalSearchTimestamps(udword& id, const Point& dir, const Point* verts, const Valencies*
 *val, udword time, udword* stamps)
 *	\param		id			[in/out] initial vertex index / final supporting vertex index
 *	\param		dir			[in] separating vector
 *	\param		verts		[in] list of vertices
 *	\param		val			[in] valencies for the list of vertices
 *	\param		time		[in] current timestamp
 *	\param		stamps		[in] array of timestamps (one/vertex)
 *	\return		true if success
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FUNCTION MESHMERIZER_API bool LocalSearchTimestamps(
        udword& id, const Point& dir, const Point* verts, const Valencies* val, udword time, udword* stamps);

#endif  // ICEHILLCLIMBING_H
