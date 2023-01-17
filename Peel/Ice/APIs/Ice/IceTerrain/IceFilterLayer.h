///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for filter layers.
 *	\file		IceFilterLayer.h
 *	\author		Pierre Terdiman
 *	\date		March, 5, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICEFILTERLAYER_H
#define ICEFILTERLAYER_H

struct ICETERRAIN_API FILTERLAYERCREATE : HEIGHTLAYERCREATE {
    FILTERLAYERCREATE() : mFilter(Idt) {}

    Matrix3x3 mFilter;
};

class ICETERRAIN_API FilterLayer : public HeightLayer {
public:
    FilterLayer();
    virtual ~FilterLayer();

    virtual bool Init(const FILTERLAYERCREATE& create);
    virtual bool Update(float* field, udword width, udword height) const;

    inline_ const Matrix3x3& GetFilter() const { return mFilter; }
    inline_ void SetFilter(const Matrix3x3& m) { mFilter = m; }

protected:
    Matrix3x3 mFilter;
};

#endif  // ICEFILTERLAYER_H
