#ifndef SPECTRAL_H_
#define SPECTRAL_H_

#include "kritapigment_export.h"

KRITAPIGMENT_EXPORT void spectralMix(float srcR, float srcG, float srcB, float factor, float* dstR, float* dstG, float* dstB);
KRITAPIGMENT_EXPORT void linearToReflectance(float r, float g, float b, float* R);
KRITAPIGMENT_EXPORT float reflectanceToLuminance(float* R);
KRITAPIGMENT_EXPORT float luminanceToConcentration(float l1, float l2, float t);
KRITAPIGMENT_EXPORT void reflectanceToXYZ(float *R, float* X, float* Y, float* Z);
KRITAPIGMENT_EXPORT void XYZToLinear(float X, float Y, float Z, float* r, float* g, float* b);

#endif