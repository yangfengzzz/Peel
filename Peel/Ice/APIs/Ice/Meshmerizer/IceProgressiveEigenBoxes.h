///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for incremental computation of best sorting axis
 *	\file		IceProgressiveEigenBoxes.h
 *	\author		Pierre Terdiman
 *	\date		October, 3, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICEPROGRESSIVEEIGENBOXES_H
#define ICEPROGRESSIVEEIGENBOXES_H

class MESHMERIZER_API ProgressiveEigenBoxes : public ProgressiveEigen {
public:
    ProgressiveEigenBoxes();
    ~ProgressiveEigenBoxes();

    bool ComputeBestSortingAxis(udword nb, const AABB* list);
};

#endif  // ICEPROGRESSIVEEIGENBOXES_H
