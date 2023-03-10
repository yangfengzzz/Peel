///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	CAME chunk for ZB2 format.
 *	\file		ChunkCAME.h
 *	\author		Pierre Terdiman
 *	\date		August, 29, 2001
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef CHUNKCAME_H
#define CHUNKCAME_H

#define CAME_VERSION 1

enum ZCB2_CAME_Flag {
    ZCB2_CAME_ORTHO = (1 << 0),
    ZCB2_CAME_INFINITE_CLIP = (1 << 1),
};

class ZCB2_API CAMEChunk : public PRSChunk {
    DECLARE_CHUNK(CAMEChunk, mCAMECore)

    DECLARE_STD_MEMBER(FOV, float)
    DECLARE_STD_MEMBER(NearClip, float)
    DECLARE_STD_MEMBER(FarClip, float)

    DECLARE_STD_FLAG(OrthoCam, ZCB2_CAME_ORTHO)
    DECLARE_STD_FLAG(InfiniteClip, ZCB2_CAME_INFINITE_CLIP)
};

#endif  // CHUNKCAME_H
