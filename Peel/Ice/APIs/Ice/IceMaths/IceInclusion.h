///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for various inclusion functions.
 *	\file		IceInclusion.h
 *	\author		Pierre Terdiman
 *	\date		September, 25, 2003
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICEINCLUSION_H
#define ICEINCLUSION_H

ICEMATHS_API IceInterval Fabs(const IceInterval& interval);
ICEMATHS_API void CosinusInclusionSnyder92(IceInterval& dst, float a, float b);
ICEMATHS_API void CosinusInclusionRevisited(IceInterval& dst, float a, float b);

#endif  // ICEINCLUSION_H
