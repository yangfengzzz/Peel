///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code related to sorting along an axis. Previously scattered everywhere in the codebase.
 *	\file		IceAxisSortHelper.h
 *	\author		Pierre Terdiman
 *	\date		July, 2, 2009
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICEAXISSORTHELPER_H
#define ICEAXISSORTHELPER_H

class ICEMATHS_API AxisSortHelper : public Allocateable {
public:
    AxisSortHelper();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     *	Checks if the axis has changed significantly since last call. If yes, sorting should be done. Otherwise not.
     *	\param		dir		[in] current direction vector (e.g. view vector)
     *	\return		true if sorting is necessary
     */
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    bool NeedsSorting(const Point& dir);

private:
    Point mCachedDir;            //!< Cached direction vector
    static float mSortingLimit;  //!< Limit |cos(angle)| beyond which sorting is performed
};

#endif  // ICEAXISSORTHELPER_H
