///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for tangent-space maths.
 *	\file		IceTextureBasis.h
 *	\author		Pierre Terdiman
 *	\date		September, 20, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICETEXTUREBASIS_H
#define ICETEXTUREBASIS_H

ICERENDERER_API bool GenerateBases(udword nb_verts,
                                   const ubyte* positions,
                                   const ubyte* uvs,
                                   udword stride,
                                   const uword* indices,
                                   udword nb_faces,
                                   Matrix3x3* basis_matrices);
ICERENDERER_API bool GenerateBases(udword nb_verts,
                                   const ubyte* positions,
                                   const ubyte* uvs,
                                   const ubyte* normals,
                                   udword stride,
                                   const uword* indices,
                                   udword nb_faces,
                                   Matrix3x3* basis_matrices);

#endif  // ICETEXTUREBASIS_H
