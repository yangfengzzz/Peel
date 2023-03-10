///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	CURV chunk for ZB2 format.
 *	\file		ChunkCURV.h
 *	\author		Pierre Terdiman
 *	\date		September, 11, 2001
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef CHUNKCURV_H
#define CHUNKCURV_H

#define CURV_VERSION 1

enum ZCB2_CURV_Flag {
    ZCB2_CURV_CLOSED = (1 << 0),  //!< true : the curve is closed
};

class ZCB2_API CURVChunk : public PNTSChunk {
    DECLARE_CHUNK(CURVChunk, mCURVCore)

    DECLARE_STD_FLAG(Closed, ZCB2_CURV_CLOSED)
};

#endif  // CHUNKCURV_H
