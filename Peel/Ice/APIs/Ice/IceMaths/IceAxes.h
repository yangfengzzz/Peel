///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains axes definition.
 *	\file		IceAxes.h
 *	\author		Pierre Terdiman
 *	\date		January, 29, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICEAXES_H
#define ICEAXES_H

enum AxisIndex {
    X_ = 0,
    Y_ = 1,
    Z_ = 2,
    W_ = 3,

    AXIS_FORCE_DWORD = 0x7fffffff
};

enum AxisOrder {
    AXES_XYZ = (X_) | (Y_ << 2) | (Z_ << 4),
    AXES_XZY = (X_) | (Z_ << 2) | (Y_ << 4),
    AXES_YXZ = (Y_) | (X_ << 2) | (Z_ << 4),
    AXES_YZX = (Y_) | (Z_ << 2) | (X_ << 4),
    AXES_ZXY = (Z_) | (X_ << 2) | (Y_ << 4),
    AXES_ZYX = (Z_) | (Y_ << 2) | (X_ << 4),

    AXES_FORCE_DWORD = 0x7fffffff
};

class ICEMATHS_API Axes : public Allocateable {
public:
    inline_ Axes(AxisOrder order) {
        mAxis0 = (order)&3;
        mAxis1 = (order >> 2) & 3;
        mAxis2 = (order >> 4) & 3;
    }
    inline_ ~Axes() {}

    udword mAxis0;
    udword mAxis1;
    udword mAxis2;
};

#endif  // ICEAXES_H
