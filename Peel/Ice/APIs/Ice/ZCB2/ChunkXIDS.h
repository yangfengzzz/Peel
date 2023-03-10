///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	XIDS chunk for ZB2 format.
 *	\file		ChunkXIDS.h
 *	\author		Pierre Terdiman
 *	\date		September, 11, 2001
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef CHUNKXIDS_H
#define CHUNKXIDS_H

#define XIDS_VERSION 1

class ZCB2_API XIDSChunk : public BaseChunk {
    DECLARE_CHUNK(XIDSChunk, mXIDSCore)

    DECLARE_STD_ARRAY(DIDs, udword)
    DECLARE_STD_ARRAY(WIDs, uword)
};

#endif  // CHUNKXIDS_H
