///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code to optimize tri-lists for cache coherence.
 *	\file		IceListOptimizer.h
 *	\author		Pierre Terdiman, Tom Forsyth
 *	\date		March, 8, 2002
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICELISTOPTIMIZER_H
#define ICELISTOPTIMIZER_H

MESHMERIZER_API bool OptimiseTriList(udword nb_faces, uword* list);
MESHMERIZER_API void OptimiseVertexCoherencyTriList(WORD* list, int nb_tris);
MESHMERIZER_API void OptimizeFacesLRU(int faceCount, short* list);

#endif  // ICELISTOPTIMIZER_H
