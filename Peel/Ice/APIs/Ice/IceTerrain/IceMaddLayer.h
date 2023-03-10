///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for multiply-add layers. This replaces the previous "constant layer" and "scale layer".
 *	\file		IceMaddLayer.h
 *	\author		Pierre Terdiman
 *	\date		December, 25, 2020
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICEMADDLAYER_H
#define ICEMADDLAYER_H

struct ICETERRAIN_API MADDLAYERCREATE : HEIGHTLAYERCREATE {
    MADDLAYERCREATE() : mAltitude(0.0f) {}

    float mAltitude;
};

class ICETERRAIN_API MaddLayer : public HeightLayer {
public:
    MaddLayer();
    virtual ~MaddLayer();

    virtual bool Init(const MADDLAYERCREATE& create);
    virtual bool Update(float* field, udword width, udword height) const;

    inline_ float GetAltitude() const { return mAltitude; }
    inline_ void SetAltitude(float h) { mAltitude = h; }

protected:
    float mAltitude;
};

#endif  // ICEMADDLAYER_H
