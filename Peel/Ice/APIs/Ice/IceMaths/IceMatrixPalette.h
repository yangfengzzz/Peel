///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for a matrix container.
 *	\file		IceMatrixPalette.h
 *	\author		Pierre Terdiman
 *	\date		April, 4, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICEMATRIXPALETTE_H
#define ICEMATRIXPALETTE_H

class ICEMATHS_API MatrixPalette : public Container {
public:
    MatrixPalette() {}
    ~MatrixPalette() {}

    inline_ udword GetNbMatrices() const { return GetNbEntries() / (sizeof(Matrix4x4) / sizeof(float)); }
    inline_ const Matrix4x4* GetMatrices() const { return (const Matrix4x4*)GetEntries(); }

    void AddMatrix(const Matrix4x4& mat) { Add((float*)mat.m, 16); }
};

CHECK_CONTAINER_ITEM(Matrix4x4)

#endif  // ICEMATRIXPALETTE_H
