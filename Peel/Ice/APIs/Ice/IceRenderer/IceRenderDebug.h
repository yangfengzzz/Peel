///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains rendering debug-helpers
 *	\file		IceRenderDebug.h
 *	\author		Pierre Terdiman
 *	\date		April, 11, 2003
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICERENDERDEBUG_H
#define ICERENDERDEBUG_H

// Forward references
class Renderer;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Checks the vertex format is compatible with current renderer settings. A message box is displayed otherwise.
 *	\fn			CheckFormat(VertexFormat fvf, Renderer* rnd)
 *	\param		fvf		[in] flexible vertex format
 *	\param		rnd		[in] renderer
 *	\return		true (Just so that I can write "ASSERT(CheckFormat(...)))
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ICERENDERER_API bool CheckFormat(VertexFormat fvf, Renderer* rnd);

#endif  // ICERENDERDEBUG_H
