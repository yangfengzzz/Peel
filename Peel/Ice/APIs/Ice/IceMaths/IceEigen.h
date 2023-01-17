///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	This file contains code for eigen vectors & eigen values. Original code from Numerical Recipes in C.
 *	\file		IceEigen.h
 *	\author		Pierre Terdiman
 *	\date		January, 29, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICEEIGEN_H
#define ICEEIGEN_H

	FUNCTION ICEMATHS_API int Eigen(Matrix3x3& vout, Point& dout, const Matrix3x3& mat);
	FUNCTION ICEMATHS_API bool EigenVectors(const Matrix3x3& M, Matrix3x3& V, Point& d);

#endif	// ICEEIGEN_H

