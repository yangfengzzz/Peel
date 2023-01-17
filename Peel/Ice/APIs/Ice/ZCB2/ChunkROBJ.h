///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	ROBJ chunk for ZB2 format.
 *	\file		ChunkROBJ.h
 *	\author		Pierre Terdiman
 *	\date		September, 13, 2001
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef CHUNKROBJ_H
#define CHUNKROBJ_H

#define ROBJ_VERSION 1

class ZCB2_API ROBJChunk : public BaseChunk {
    DECLARE_CHUNK(ROBJChunk, mROBJCore)

    DECLARE_SUBCHUNK(Controllers, XIDSChunk)
};

#endif  // CHUNKROBJ_H
