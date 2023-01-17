///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for bitmap layers.
 *	\file		IceBitmapLayer.h
 *	\author		Pierre Terdiman
 *	\date		March, 5, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICEBITMAPLAYER_H
#define ICEBITMAPLAYER_H

/*
        struct ICETERRAIN_API CONSTANTLAYERCREATE : HEIGHTLAYERCREATE
        {
                float		ConstantHeight;		//!< Constant height
        };

        class ICETERRAIN_API ConstantLayer ; public HeightLayer
        {
                public:
                                                                                ConstantLayer();
                virtual							~ConstantLayer();

                virtual			bool			Init(const CONSTANTLAYERCREATE* create);
                virtual			bool			Update(float* field, udword width, udword height);

                // Data access
                inline_	udword			GetConstantHeight()			{ return mConstantHeight;
}

                protected:
                                                float			mConstantHeight;
        };

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        //
Bitmap Layer Class Definition
        //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        struct BITMAPLAYERCREATE : HEIGHTLAYERCREATE
        {
                                        ubyte*				BitmapHeight;		// Initialize from a
byte array...
//					Picture*			PictureHeight;		// ...or from a picture
        };

        class BitmapLayer : public HeightLayer
        {
                public:
                                                                                BitmapLayer();
                        virtual						~BitmapLayer();

                        // Run-time type information (::Node)
                        DECLARENODE(BitmapLayer, HeightLayer, &guidClass3DDatabase, &guidSubClassGroup,
&guidNodeBitmapLayer)

                        virtual	bool				Init(HEIGHTLAYERCREATE* create);
                        virtual HeightLayer&		Update(float* field, udword width, udword height);

                protected:
                        // Bitmap parameters
                                        ubyte*				mBitmapHeights;
        };
*/

#endif  // ICEBITMAPLAYER_H
