///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code to quantize vertices.
 *	\file		IceQuantize.h
 *	\author		Pierre Terdiman
 *	\date		January, 29, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICEQUANTIZE_H
#define ICEQUANTIZE_H

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Computes quantization coeffs for a vertex cloud.
 *	\relates	VertexCloud
 *	\fn			ComputeQuantizationCoeffs(udword nb_bits, udword nb_pts, Point* pts, Point& quantizer)
 *	\param		nb_bits			[in] number of bits used to quantize a floating-point value (max is 15,
 *leaving a sign bit) \param		nb_pts			[in] number of vertices to quantize
 *	\param		pts				[in] list of vertices
 *	\param		quantizer		[out] quantization coeffs
 *	\return		true if success
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FUNCTION MESHMERIZER_API bool ComputeQuantizationCoeffs(udword nb_bits,
                                                        udword nb_pts,
                                                        const Point* pts,
                                                        Point& quantizer);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Quantizes a list of vertices.
 *	\relates	VertexCloud
 *	\fn			Quantize(udword nb_bits, udword nb_pts, Point* pts, sword* quantized, Point&
 *dequant_coeff) \param		nb_bits			[in] number of bits used to quantize a floating-point value (max
 *is 15, leaving a sign bit) \param		nb_pts			[in] number of vertices to quantize \param
 *pts				[in/out] list of vertices
 *	\param		quantized		[out] a place to store nb_pts quantized vertices, or null. If null, the
 *original vertex cloud is updated \param		dequant_coeff	[out] dequantization coeffs
 *	\return		true if success
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FUNCTION MESHMERIZER_API bool Quantize(
        udword nb_bits, udword nb_pts, Point* pts, sword* quantized, Point& dequant_coeff);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Dequantizes a list of vertices.
 *	\relates	VertexCloud
 *	\fn			Dequantize(udword nb_pts, const sword* quantized, Point* pts, const Point&
 *dequant_coeff) \param		nb_pts			[in] number of vertices to dequantize
 *	\param		quantized		[in] list of quantized vertices
 *	\param		pts				[out] a place to store nb_pts dequantized vertices
 *	\param		dequant_coeff	[out] dequantization coeffs
 *	\return		true if success
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FUNCTION MESHMERIZER_API bool Dequantize(udword nb_pts, const sword* quantized, Point* pts, const Point& dequant_coeff);

#endif  // ICEQUANTIZE_H
