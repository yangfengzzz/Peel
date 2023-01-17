///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	TCAM chunk for ZB2 format.
 *	\file		ChunkTCAM.h
 *	\author		Pierre Terdiman
 *	\date		August, 30, 2001
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef CHUNKTCAM_H
#define CHUNKTCAM_H

	#define TCAM_VERSION	1

	class ZCB2_API TCAMChunk : public CAMEChunk
	{
		DECLARE_CHUNK(TCAMChunk, mTCAMCore)
	};

#endif // CHUNKTCAM_H