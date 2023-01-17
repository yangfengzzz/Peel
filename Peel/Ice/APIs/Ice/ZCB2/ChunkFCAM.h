///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	FCAM chunk for ZB2 format.
 *	\file		ChunkFCAM.h
 *	\author		Pierre Terdiman
 *	\date		August, 30, 2001
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef CHUNKFCAM_H
#define CHUNKFCAM_H

	#define FCAM_VERSION	1

	class ZCB2_API FCAMChunk : public CAMEChunk
	{
		DECLARE_CHUNK(FCAMChunk, mFCAMCore)
	};

#endif // CHUNKFCAM_H