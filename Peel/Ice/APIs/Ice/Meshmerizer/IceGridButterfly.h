///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains a hardcoded Butterfly subdivision for a grid. (a.k.a. on-the-fly Butterfly)
 *	\file		IceGridButterfly.h
 *	\author		Pierre Terdiman
 *	\date		January, 29, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICEGRIDBUTTERFLY_H
#define ICEGRIDBUTTERFLY_H

//! Only works with this FVF
struct MESHMERIZER_API GridVertex {
    Point p;
    Point n;
    float u, v;
};

FUNCTION MESHMERIZER_API bool GridButterfly(udword nbu,
                                            udword nbv,
                                            const Point* p,
                                            const Point* n,
                                            GridVertex* vertices,
                                            uword* faces,
                                            bool& init_done,
                                            udword& nb_verts,
                                            udword& nb_indices);

#endif  // ICEGRIDBUTTERFLY_H
