///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Compilation flags for IceCharacter.
 *	\file		IceCharacterSettings.h
 *	\author		Pierre Terdiman
 *	\date		October, 19, 2002
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICECHARACTERSETTINGS_H
#define ICECHARACTERSETTINGS_H

#define MAX_BASIC_EVENTS 2  // 8
#define SUPPORT_COMPLEX_EVENTS
//	#define SUPPORT_TIME_BLEND_MAPPING
#define SUPPORT_HEIGHT_TRACK
// #define SUPPORT_MOTION_INTERVAL
#define SUPPORT_MOTION_MODIFIERS 16
//	#define SUPPORT_WEAPON_MODIFIERS	2
#define SUPPORT_WEAPON_MODIFIERS2 3
//	#define SUPPORT_MOTION_LINKS		32

// New design too painful with target override => remove this feature
//	#define SUPPORT_CELL_ACTIVE_TRANS
//	#define SUPPORT_TARGET_OVERRIDE
#define SUPPORT_NEW_DESIGN
#define SUPPORT_NEW_DESIGN2

#define SUPPORT_NEW_MOTIONDATA_DESIGN

#endif  // ICECHARACTERSETTINGS_H
