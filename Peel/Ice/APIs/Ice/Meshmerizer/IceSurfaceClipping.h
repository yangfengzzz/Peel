///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for surface clipping
 *	\file		IceSurfaceClipping.h
 *	\author		Pierre Terdiman
 *	\date		February, 20, 2007
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICESURFACECLIPPING_H
#define ICESURFACECLIPPING_H

typedef void (*ClippedMeshCallback)(const FatTriSurface& mesh, const udword* mapping, void* user_data);

FUNCTION MESHMERIZER_API bool CreatedClippedMesh(const IndexedSurface* surface,
                                                 const IndexedSurface* uv_surface,
                                                 const IndexedSurface* rgb_surface,
                                                 const Matrix4x4* world,
                                                 const Plane* clip_planes,
                                                 const AABB* clip_box,
                                                 ClippedMeshCallback callback,
                                                 void* user_data);

#endif  // ICESURFACECLIPPING_H
