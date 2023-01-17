///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains base code for field controllers.
 *	\file		IceFieldController.h
 *	\author		Pierre Terdiman
 *	\date		April, 4, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICEFIELDCONTROLLER_H
#define ICEFIELDCONTROLLER_H

enum CtrlMode {
    CTRL_SAMPLES = 1,     //!< Samples
    CTRL_KEYFRAMES = 2,   //!< Keyframes
    CTRL_PROCEDURAL = 3,  //!< Procedural

    CTRL_FORCE_DWORD = 0x7fffffff
};

//! A field controller
class ICECORE_API FieldController {
public:
    inline_ FieldController() {}
    virtual ~FieldController() {}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     *	Gets controlled type.
     *	\return		controlled type
     */
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    virtual FieldType GetControlledType() const = 0;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     *	Gets number of controlled items.
     *	\return		number of controlled items
     */
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    virtual udword GetNbControlledItems() const { return 1; /* Default is 1 item */ }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     *	Checks the controller is compatible with a given field descriptor.
     *	\param		field		[in] field descriptor
     *	\return		true if both are compatible
     */
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    inline_ bool IsCompatible(const FieldDescriptor& field) const {
        if (GetControlledType() != field.Type) return false;
        //							if(GetNbControlledItems()!=field.GetCount())	return
        //false;
        return true;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     *	Gets a value for a given frame.
     *	\param		time		[in] time-related info
     *	\param		dest		[out] a possible user destination buffer
     *	\param		owner		[in] a possible owner, as procedural controller sometimes depends on the caller
     *	\return		address of the correct value (or null)
     */
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    virtual const void* GetValue(const TimeInfo& time, void* dest = null, Cell* owner = null) = 0;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     *	Executes a controller in 2 possible ways:
     *	- using an accessor if one is available (useful for lazy-evaluation, etc)
     *	- by directly modifying the field within the class
     *
     *	\param		owner		[in] owner cell to modify
     *	\param		field		[in] descriptor for the field to be modified
     *	\param		time		[in] time-related info
     *	\return		true if success
     */
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    inline_ bool Apply(Cell* owner, const FieldDescriptor& field, const TimeInfo& time) {
        if (!owner) return false;

        // Use the accessor if available
        if (field.Accessor) {
            return (field.Accessor)(owner, this, time);
        }

        // Else use the generic slower path

        // Get new data
        const void* Data = GetValue(time);
        if (Data) field.Replace(owner, Data);
        return true;
    }
};

/*
        //! A field controller
        class ICECORE_API FieldController
        {
                public:
                inline_						FieldController(FieldType type=FIELD_FORCE_DWORD)
   { mControlledType = type;	} virtual						~FieldController()
   {}

                inline_			void		SetControlledType(FieldType type)		{
   mControlledType = type; return *this; } inline_			FieldType	GetControlledType()
   { return mControlledType;				}

                                bool				IsCompatible(FieldDescriptor& field)	{ return
   mControlledType==field.Type;	}

                virtual	bool				Update(Cell* cell, FieldDescriptor* field)	= 0;

                private:
                                FieldType			mControlledType;
        };
*/

/*
        class ICECORE_API SinusFloatController : public FieldController
        {
                public:
                                                                        SinusFloatController()	{
   SetControlledType(FIELD_FLOAT); mTime = 0.0f; mAmplitude = 10.0f; }
                virtual						~SinusFloatController()	{}

                virtual	bool				Update(Cell* cell, FieldDescriptor* field)
                                {
                                        mTime += 0.01f;
                                        float Value = mAmplitude * sinf(mTime);

                                        field->Replace(cell, &Value);

                                        return true;
                                }

                private:
                                float				mTime;
                                float				mAmplitude;
        };
*/

#endif  // ICEFIELDCONTROLLER_H
