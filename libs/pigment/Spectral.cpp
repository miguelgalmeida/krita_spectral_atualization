/*
    This is almost a direct port of Spectral.js by Ronald van Wijnen.
    Spectral.js can be found at: https://github.com/rvanwijnen/spectral.js

    10-11-2023 ronald.van.wijnen@gmail.com
*/

#include "Spectral.h"
#include <cmath>

static const int SIZE = 36;
static const float EPSILON = 0.00000001;

//Generated red, green and blue spectral data
//For reference but not used
//http://scottburns.us/fast-rgb-to-spectrum-conversion-for-reflectances/
static const float SPD_R[SIZE] = {0.02159246, 0.02029311, 0.02180791, 0.0238033, 0.02520813, 0.02541496, 0.02462128, 0.02097371, 0.0157528, 0.01116804, 0.00857828, 0.00658188, 0.00517172, 0.00454521, 0.00414512, 0.00434311, 0.00523816, 0.00725194, 0.01254366, 0.02806713, 0.09134228, 0.48408109, 0.87037832, 0.93951313, 0.96092699, 0.96862376, 0.97126388, 0.97228582, 0.97189874, 0.97269186, 0.97173481, 0.97234454, 0.97150339, 0.970858, 0.97055387, 0.9696714};
static const float SPD_G[SIZE] = {0.01054241, 0.01087898, 0.01106351, 0.01073657, 0.01168181, 0.01243472, 0.01498691, 0.02010039, 0.03035626, 0.06338896, 0.17342384, 0.56832114, 0.827792, 0.91656047, 0.95200284, 0.96409645, 0.97059086, 0.97250254, 0.9691482, 0.95534465, 0.89263723, 0.5003641, 0.11623672, 0.04795139, 0.02787353, 0.02005796, 0.01738217, 0.01542911, 0.01543808, 0.01454683, 0.01519777, 0.0142859, 0.01506912, 0.01550626, 0.0155458, 0.01630284};
static const float SPD_B[SIZE] = {0.96786514, 0.96882791, 0.96712858, 0.96546014, 0.96311006, 0.96215032, 0.96039181, 0.9589259, 0.95389094, 0.925443, 0.81799789, 0.42509696, 0.16703627, 0.07889433, 0.04385204, 0.03156044, 0.02417098, 0.02024552, 0.01830814, 0.01658822, 0.01602049, 0.01555481, 0.01338496, 0.01253549, 0.01119948, 0.01131827, 0.01135395, 0.01228507, 0.01266319, 0.01276133, 0.01306743, 0.01336957, 0.01342749, 0.01363574, 0.0138936, 0.01402576};

//CIE (1964-10 degree) color matching functions multiplied by the CIE D65 standard illuminant
//normalized by dividing the values by the sum of CIE_CMF_Y * D65
//these are called X-, Y- and Z-bar.
static const float X_BAR[SIZE] = {0.00006469, 0.00021942, 0.0011206, 0.00376671, 0.01188085, 0.02328702, 0.03456028, 0.03722472, 0.03241918, 0.02123373, 0.01049125, 0.00329592, 0.00050705, 0.0009487, 0.00627387, 0.01686504, 0.02869036, 0.04267588, 0.05625615, 0.06947213, 0.08305522, 0.08612824, 0.09046839, 0.08500598, 0.07090844, 0.05063015, 0.03547485, 0.02146875, 0.01251677, 0.00680475, 0.00346465, 0.00149765, 0.00076972, 0.00040738, 0.00016901, 0.00009523};
static const float Y_BAR[SIZE] = {0.00000184, 0.00000621, 0.00003101, 0.00010475, 0.00035365, 0.0009515, 0.00228232, 0.00420743, 0.00668897, 0.00988864, 0.01524983, 0.02141884, 0.03342376, 0.05131129, 0.07040384, 0.0878409, 0.0942514, 0.09795911, 0.09415453, 0.08678319, 0.0788585, 0.06352829, 0.05374276, 0.04264713, 0.03161814, 0.02088573, 0.01386046, 0.00810284, 0.00463022, 0.00249144, 0.00125933, 0.00054166, 0.00027796, 0.00014711, 0.00006103, 0.00003439};
static const float Z_BAR[SIZE] = {0.00030502, 0.00103683, 0.00531327, 0.01795484, 0.057079, 0.11365445, 0.17336305, 0.19621147, 0.18608701, 0.13995396, 0.08917675, 0.04789741, 0.02814633, 0.01613806, 0.0077593, 0.00429626, 0.00200556, 0.00086149, 0.00036905, 0.00019143, 0.00014956, 0.00009231, 0.00006814, 0.00002883, 0.00001577, 0.00000394, 0.00000158, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

//XYZ to RGB matrix, rounded on 8 decimals
//https://github.com/w3c/csswg-drafts/issues/5922
static const float XYZ_RGB[3][3] = {
    {3.24096994, -1.53738318, -0.49861076},
    {-0.96924364, 1.8759675, 0.04155506},
    {0.05563008, -0.20397696, 1.05697151}
};

void spectralMix(float srcR, float srcG, float srcB, float factor, float* dstR, float* dstG, float* dstB) {
    float R1[36], R2[36], R[36];
    
    linearToReflectance(srcR, srcG, srcB, R1);
    linearToReflectance(*dstR, *dstG, *dstB, R2);

    float l1 = reflectanceToLuminance(R1);
    float l2 = reflectanceToLuminance(R2);

    float c = luminanceToConcentration(l1, l2, factor);
    
    for (int i = 0; i < SIZE; i++) {
        float KS = 0.0;

        KS += (pow(1.0 - R1[i], 2) / (2 * R1[i])) * (1.0 - c);
        KS += (pow(1.0 - R2[i], 2) / (2 * R2[i])) * c;

        float KM = 1.0 + KS - sqrt(pow(KS, 2) + 2 * KS);

        R[i] = KM;
    }

    float X = 0.0, Y = 0.0, Z = 0.0;

    reflectanceToXYZ(R, &X, &Y, &Z);

    XYZToLinear(roundf(X * 1000) / 1000, roundf(Y * 1000) / 1000, roundf(Z * 1000) / 1000, dstR, dstG, dstB);
}

//function for creating a SPD out of the red, green and blue SPD's
void linearToReflectance(float r, float g, float b, float* R) {
    float weightW = fmin(r, fmin(g, b));
    
    float weightR = r - weightW;
    float weightG = g - weightW;
    float weightB = b - weightW;

    for (int i = 0; i < SIZE; i++) {
      R[i] = fmax(EPSILON, weightW + weightR * SPD_R[i] + weightG * SPD_G[i] + weightB * SPD_B[i]);
    }
}

//luminance is the sum of multiplying the reflectance with the Y-bar
float reflectanceToLuminance(float* R) {
    float l = 0.0;

    for (int i = 0; i < SIZE; i++) {
      l += R[i] * Y_BAR[i];
    }

    return l;
}

//formula for calculating the mix factor based on the luminance
//this is needed for generating a perceptually even spread between 0 and 1
//where 0 is color1 and 1 is color2 
float luminanceToConcentration(float l1, float l2, float t) {
    float t1 = l1 * pow(1.0 - t, 2);
    float t2 = l2 * pow(t, 2);

    return t2 / (t1 + t2);
}

void reflectanceToXYZ(float *R, float* X, float* Y, float* Z) {
    for (int i = 0; i < SIZE; i++) {
      *X += R[i] * X_BAR[i];
      *Y += R[i] * Y_BAR[i];
      *Z += R[i] * Z_BAR[i];
    }
}

void XYZToLinear(float X, float Y, float Z, float* r, float* g, float* b) {
    *r = X * XYZ_RGB[0][0] + Y * XYZ_RGB[0][1] + Z * XYZ_RGB[0][2];
    *g = X * XYZ_RGB[1][0] + Y * XYZ_RGB[1][1] + Z * XYZ_RGB[1][2];
    *b = X * XYZ_RGB[2][0] + Y * XYZ_RGB[2][1] + Z * XYZ_RGB[2][2];
}
