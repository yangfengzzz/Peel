///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for lightnings.
 *	\file		IceLightning.h
 *	\author		Pierre Terdiman
 *	\date		March, 5, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICELIGHTNING_H
#define ICELIGHTNING_H

	ICETERRAIN_API udword	ComputeLightningNbPts(udword level);
	ICETERRAIN_API bool		GenerateLightning(udword seed, udword level, float dispersion, const Point& p0, const Point& p1, Vertices& lightning);

#endif // ICELIGHTNING_H