///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains path related code.
 *	\file		IcePath.h
 *	\author		Pierre Terdiman
 *	\date		April, 4, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICEPATH_H
#define ICEPATH_H

struct VirtualFileHandle;

FUNCTION ICECORE_API bool RegisterPath(const char* path);
FUNCTION ICECORE_API bool ReleasePaths();
FUNCTION ICECORE_API bool FindFile(const char* filename, VirtualFileHandle* handle = null, bool search_archives = true);

FUNCTION ICECORE_API bool RemovePath(String& filename);
//	FUNCTION ICECORE_API	bool		GetPath(const String& filename, String& path);
FUNCTION ICECORE_API bool GetPath(const char* filename, String& path);
FUNCTION ICECORE_API bool GetPath2(const char* filename, String& path);

#endif  // ICEPATH_H
