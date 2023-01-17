
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef CHUNKDLIT_H
#define CHUNKDLIT_H

#define DLIT_VERSION 1

class ZCB2_API DLITChunk : public LITEChunk {
    DECLARE_CHUNK(DLITChunk, mDLITCore)

    DECLARE_STD_MEMBER(TDist, float)
};

#endif  // CHUNKDLIT_H
