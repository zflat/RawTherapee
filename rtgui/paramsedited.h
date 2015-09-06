/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _PARAMEDITED_H_
#define _PARAMEDITED_H_

#include <glibmm.h>
#include <vector>
#include "../rtengine/procparams.h"
#include "../rtengine/rtengine.h"

class GeneralParamsEdited
{

public:
    bool rank;
    bool colorlabel;
    bool intrash;

    operator bool (void) const
    {
        printf("bool: %d\n", rank|colorlabel|intrash);
        return rank|colorlabel|intrash;
    }
    void set(bool v);
};

class ToneCurveParamsEdited
{

public:
    bool curve;
    bool curve2;
    bool curveMode;
    bool curveMode2;
    bool brightness;
    bool black;
    bool contrast;
    bool saturation;
    bool shcompr;
    bool hlcompr;
    bool hlcomprthresh;
    bool autoexp;
    bool clip;
    bool expcomp;
    bool hrenabled;
    bool method;

    operator bool (void) const
    {
        printf("bool: %d\n", curve|curve2|curveMode|curveMode2|brightness|black|contrast|saturation|shcompr|hlcompr|hlcomprthresh|autoexp|clip|expcomp|hrenabled|method);
        return curve|curve2|curveMode|curveMode2|brightness|black|contrast|saturation|shcompr|hlcompr|hlcomprthresh|autoexp|clip|expcomp|hrenabled|method;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class LCurveParamsEdited
{
public:
    bool brightness;
    bool contrast;
    bool chromaticity;
    bool avoidcolorshift;
    bool rstprotection;
    bool lcurve;
    bool acurve;
    bool bcurve;
    bool lcredsk;
    bool cccurve;
    bool chcurve;
    bool lhcurve;
    bool hhcurve;
    bool lccurve;
    bool clcurve;

    bool enabled;
    bool method;

    operator bool (void) const
    {
        printf("bool: %d\n", brightness|contrast|chromaticity|avoidcolorshift|rstprotection|lcurve|acurve|bcurve|lcredsk|cccurve|chcurve|lhcurve|hhcurve|lccurve|clcurve|enabled|method);
        return brightness|contrast|chromaticity|avoidcolorshift|rstprotection|lcurve|acurve|bcurve|lcredsk|cccurve|chcurve|lhcurve|hhcurve|lccurve|clcurve|enabled|method;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class RGBCurvesParamsEdited
{

public:
    bool lumamode;
    bool rcurve;
    bool gcurve;
    bool bcurve;

    operator bool (void) const
    {
        printf("bool: %d\n", lumamode|rcurve|gcurve|bcurve);
        return lumamode|rcurve|gcurve|bcurve;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class ColorToningEdited
{

public:
    bool enabled;
    bool opacityCurve;
    bool colorCurve;
    bool clcurve;
    bool method;
    bool autosat;
    bool satprotectionthreshold;
    bool saturatedopacity;
    bool strength;
    bool shadowsColSat;
    bool hlColSat;
    bool balance;
    bool twocolor;
    bool cl2curve;
    bool redlow;
    bool greenlow;
    bool bluelow;
    bool redmed;
    bool greenmed;
    bool bluemed;
    bool redhigh;
    bool greenhigh;
    bool bluehigh;
    bool satlow;
    bool sathigh;
    bool lumamode;

    operator bool (void) const
    {
        printf("bool: %d\n", enabled|opacityCurve|colorCurve|clcurve|method|autosat|satprotectionthreshold|saturatedopacity|strength|shadowsColSat|hlColSat|balance|twocolor|cl2curve|redlow|greenlow|bluelow|redmed|greenmed|bluemed|redhigh|greenhigh|bluehigh|satlow|sathigh|lumamode);
        return enabled|opacityCurve|colorCurve|clcurve|method|autosat|satprotectionthreshold|saturatedopacity|strength|shadowsColSat|hlColSat|balance|twocolor|cl2curve|redlow|greenlow|bluelow|redmed|greenmed|bluemed|redhigh|greenhigh|bluehigh|satlow|sathigh|lumamode;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class SharpenEdgeParamsEdited
{

public :
    bool enabled;
    bool passes;
    bool amount;
    bool threechannels;

    operator bool (void) const
    {
        printf("bool: %d\n", enabled|passes|amount|threechannels);
        return enabled|passes|amount|threechannels;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class SharpenMicroParamsEdited
{
public :
    bool enabled;
    bool matrix;
    bool amount;
    bool uniformity;

    operator bool (void) const
    {
        printf("bool: %d\n", enabled|matrix|amount|uniformity);
        return enabled|matrix|amount|uniformity;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class SharpeningParamsEdited
{

public:
    bool enabled;
    bool radius;
    bool amount;
    bool threshold;
    bool edgesonly;
    bool edges_radius;
    bool edges_tolerance;
    bool halocontrol;
    bool halocontrol_amount;

    bool method;
    bool deconvamount;
    bool deconvradius;
    bool deconviter;
    bool deconvdamping;

    operator bool (void) const
    {
        printf("bool: %d\n", enabled|radius|amount|threshold|edgesonly|edges_radius|edges_tolerance|halocontrol|halocontrol_amount|method|deconvamount|deconvradius|deconviter|deconvdamping);
        return enabled|radius|amount|threshold|edgesonly|edges_radius|edges_tolerance|halocontrol|halocontrol_amount|method|deconvamount|deconvradius|deconviter|deconvdamping;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class VibranceParamsEdited
{

public:
    bool enabled;
    bool pastels;
    bool saturated;
    bool psthreshold;
    bool protectskins;
    bool avoidcolorshift;
    bool pastsattog;
    bool skintonescurve;

    operator bool (void) const
    {
        printf("bool: %d\n", enabled|pastels|saturated|psthreshold|protectskins|avoidcolorshift|pastsattog|skintonescurve);
        return enabled|pastels|saturated|psthreshold|protectskins|avoidcolorshift|pastsattog|skintonescurve;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class WBParamsEdited
{

public:
    bool method;
    bool temperature;
    bool green;
    bool equal;

    operator bool (void) const
    {
        printf("bool: %d\n", method|temperature|green|equal);
        return method|temperature|green|equal;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class DefringeParamsEdited
{

public:
    bool enabled;
    bool radius;
    bool threshold;
    bool huecurve;

    operator bool (void) const
    {
        printf("bool: %d\n", enabled|radius|threshold|huecurve);
        return enabled|radius|threshold|huecurve;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class ImpulseDenoiseParamsEdited
{

public:
    bool enabled;
    bool thresh;

    operator bool (void) const
    {
        printf("bool: %d\n", enabled|thresh);
        return enabled|thresh;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class ColorAppearanceParamsEdited
{

public:
    bool curve;
    bool curve2;
    bool curve3;
    bool curveMode;
    bool curveMode2;
    bool curveMode3;
    bool enabled;
    bool degree;
    bool autodegree;
    bool autoadapscen;
    bool surround;
    bool adapscen;
    bool adaplum;
    bool badpixsl;
    bool wbmodel;
    bool algo;
    bool jlight;
    bool qbright;
    bool chroma;
    bool schroma;
    bool mchroma;
    bool contrast;
    bool qcontrast;
    bool colorh;
    bool rstprotection;
    bool surrsource;
    bool gamut;
//  bool badpix;
    bool datacie;
    bool tonecie;
//  bool sharpcie;

    operator bool (void) const
    {
        printf("bool: %d\n", curve|curve2|curve3|curveMode|curveMode2|curveMode3|enabled|degree|autodegree|autoadapscen|surround|adapscen|adaplum|badpixsl|wbmodel|algo|jlight|qbright|chroma|schroma|mchroma|contrast|qcontrast|colorh|rstprotection|surrsource|gamut/*|badpix*/|datacie|tonecie/*|sharpcie*/);
        return curve|curve2|curve3|curveMode|curveMode2|curveMode3|enabled|degree|autodegree|autoadapscen|surround|adapscen|adaplum|badpixsl|wbmodel|algo|jlight|qbright|chroma|schroma|mchroma|contrast|qcontrast|colorh|rstprotection|surrsource|gamut/*|badpix*/|datacie|tonecie/*|sharpcie*/;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class DirPyrDenoiseParamsEdited
{

public:
    bool enabled;
    bool enhance;
    bool median;
    bool autochroma;
    bool Ldetail;
    bool luma;
    bool chroma;
    bool redchro;
    bool bluechro;
    bool gamma;
    bool lcurve;
    bool cccurve;

//    bool perform;
    bool dmethod;
    bool Lmethod;
    bool Cmethod;
    bool C2method;
    bool smethod;
    bool medmethod;
    bool methodmed;
    bool rgbmethod;
    bool passes;

    operator bool (void) const
    {
        printf("bool: %d\n", enabled|enhance|median|autochroma|Ldetail|luma|chroma|redchro|bluechro|gamma|lcurve|cccurve/*|perform*/|dmethod|Lmethod|Cmethod|C2method|smethod|medmethod|methodmed|rgbmethod|passes);
        return enabled|enhance|median|autochroma|Ldetail|luma|chroma|redchro|bluechro|gamma|lcurve|cccurve/*|perform*/|dmethod|Lmethod|Cmethod|C2method|smethod|medmethod|methodmed|rgbmethod|passes;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class EPDParamsEdited
{
public:
    bool enabled;
    bool strength;
    bool gamma;
    bool edgeStopping;
    bool scale;
    bool reweightingIterates;

    operator bool (void) const
    {
        printf("bool: %d\n", enabled|strength|gamma|edgeStopping|scale|reweightingIterates);
        return enabled|strength|gamma|edgeStopping|scale|reweightingIterates;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};


class SHParamsEdited
{

public:
    bool enabled;
    bool hq;
    bool highlights;
    bool htonalwidth;
    bool shadows;
    bool stonalwidth;
    bool localcontrast;
    bool radius;

    operator bool (void) const
    {
        printf("bool: %d\n", enabled|hq|highlights|htonalwidth|shadows|stonalwidth|localcontrast|radius);
        return enabled|hq|highlights|htonalwidth|shadows|stonalwidth|localcontrast|radius;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class CropParamsEdited
{

public:
    bool enabled;
    bool x;
    bool y;
    bool w;
    bool h;
    bool fixratio;
    bool ratio;
    bool orientation;
    bool guide;

    operator bool (void) const
    {
        printf("bool: %d\n", enabled|x|y|w|h|fixratio|ratio|orientation|guide);
        return enabled|x|y|w|h|fixratio|ratio|orientation|guide;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class CoarseTransformParamsEdited
{

public:
    bool rotate;
    bool hflip;
    bool vflip;

    operator bool (void) const
    {
        printf("bool: %d\n", rotate|hflip|vflip);
        return rotate|hflip|vflip;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class CommonTransformParamsEdited
{

public:
    bool autofill;

    operator bool (void) const
    {
        printf("bool: %d\n", autofill);
        return autofill;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class RotateParamsEdited
{

public:
    bool degree;

    operator bool (void) const
    {
        printf("bool: %d\n", degree);
        return degree;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class DistortionParamsEdited
{

public:
    bool amount;

    operator bool (void) const
    {
        printf("bool: %d\n", amount);
        return amount;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class LensProfParamsEdited
{
public:
    bool lcpFile, useDist, useVign, useCA;

    operator bool (void) const
    {
        printf("bool: %d\n", lcpFile|useDist|useVign|useCA);
        return lcpFile|useDist|useVign|useCA;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class PerspectiveParamsEdited
{

public:
    bool horizontal;
    bool vertical;

    operator bool (void) const
    {
        printf("bool: %d\n", horizontal|vertical);
        return horizontal|vertical;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class GradientParamsEdited
{

public:
    bool enabled;
    bool degree;
    bool feather;
    bool strength;
    bool centerX;
    bool centerY;

    operator bool (void) const
    {
        printf("bool: %d\n", enabled|degree|feather|strength|centerX|centerY);
        return enabled|degree|feather|strength|centerX|centerY;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class PCVignetteParamsEdited
{

public:
    bool enabled;
    bool strength;
    bool feather;
    bool roundness;

    operator bool (void) const
    {
        printf("bool: %d\n", enabled|strength|feather|roundness);
        return enabled|strength|feather|roundness;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class VignettingParamsEdited
{

public:
    bool amount;
    bool radius;
    bool strength;
    bool centerX;
    bool centerY;

    operator bool (void) const
    {
        printf("bool: %d\n", amount|radius|strength|centerX|centerY);
        return amount|radius|strength|centerX|centerY;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class ChannelMixerParamsEdited
{

public:
    bool red[3];
    bool green[3];
    bool blue[3];

    operator bool (void) const
    {
        bool retVal = false;
        for (size_t i=0; i<3; ++i) {
            retVal |= red[i]|green[i]|blue[i];
        }
        printf("bool: %d\n", retVal);
        return retVal;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class BlackWhiteParamsEdited
{

public:
    bool enabledcc;
    bool enabled;
    bool method;
    bool filter;
    bool setting;
    bool mixerRed;
    bool mixerOrange;
    bool mixerYellow;
    bool mixerGreen;
    bool mixerCyan;
    bool mixerBlue;
    bool mixerMagenta;
    bool mixerPurple;
    bool gammaRed;
    bool gammaGreen;
    bool gammaBlue;
    bool luminanceCurve;
    bool beforeCurve;
    bool beforeCurveMode;
    bool afterCurve;
    bool afterCurveMode;
    bool autoc;
    bool algo;

    operator bool (void) const
    {
        printf("bool: %d\n", enabledcc|enabled|method|filter|setting|mixerRed|mixerOrange|mixerYellow|mixerGreen|mixerCyan|mixerBlue|mixerMagenta|mixerPurple|gammaRed|gammaGreen|gammaBlue|luminanceCurve|beforeCurve|beforeCurveMode|afterCurve|afterCurveMode|autoc|algo);
        return enabledcc|enabled|method|filter|setting|mixerRed|mixerOrange|mixerYellow|mixerGreen|mixerCyan|mixerBlue|mixerMagenta|mixerPurple|gammaRed|gammaGreen|gammaBlue|luminanceCurve|beforeCurve|beforeCurveMode|afterCurve|afterCurveMode|autoc|algo;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class CACorrParamsEdited
{

public:
    bool red;
    bool blue;

    operator bool (void) const
    {
        printf("bool: %d\n", red|blue);
        return red|blue;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class ResizeParamsEdited
{

public:
    bool scale;
    bool appliesTo;
    bool method;
    bool dataspec;
    bool width;
    bool height;
    bool enabled;

    operator bool (void) const
    {
        printf("bool: %d\n", scale|appliesTo|method|dataspec|width|height|enabled);
        return scale|appliesTo|method|dataspec|width|height|enabled;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class ColorManagementParamsEdited
{

public:
    bool input;
    bool toneCurve;
    bool applyLookTable;
    bool applyBaselineExposureOffset;
    bool applyHueSatMap;
    bool blendCMSMatrix;
    bool dcpIlluminant;
    bool working;
    bool output;
    bool gamma;
    bool gampos;
    bool slpos;
    bool gamfree;
    bool freegamma;

    operator bool (void) const
    {
        printf("bool: %d\n", input|toneCurve|applyLookTable|applyBaselineExposureOffset|applyHueSatMap|blendCMSMatrix|dcpIlluminant|working|output|gamma|gampos|slpos|gamfree|freegamma);
        return input|toneCurve|applyLookTable|applyBaselineExposureOffset|applyHueSatMap|blendCMSMatrix|dcpIlluminant|working|output|gamma|gampos|slpos|gamfree|freegamma;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class WaveletParamsEdited
{

public:
    bool enabled;
    bool strength;
    bool balance;
    bool iter;
    bool median;
    bool medianlev;
    bool linkedg;
    bool cbenab;
    bool lipst;
    bool Medgreinf;
    bool avoid;
    bool tmr;
    bool c[9];
    bool ch[9];
    bool Lmethod;
    bool CHmethod;
    bool CHSLmethod;
    bool EDmethod;
    bool BAmethod;
    bool NPmethod;
    bool TMmethod;
    bool HSmethod;
    bool CLmethod;
    bool Backmethod;
    bool Tilesmethod;
    bool daubcoeffmethod;
    bool Dirmethod;
    bool rescon;
    bool resconH;
    bool reschro;
    bool tmrs;
    bool gamma;
    bool sup;
    bool sky;
    bool thres;
    bool threshold;
    bool threshold2;
    bool edgedetect;
    bool edgedetectthr;
    bool edgedetectthr2;
    bool edgesensi;
    bool edgeampli;
    bool chro;
    bool chroma;
    bool contrast;
    bool edgrad;
    bool edgval;
    bool edgthresh;
    bool thr;
    bool thrH;
    bool skinprotect;
    bool hueskin;
    bool hueskin2;
    bool hllev;
    bool bllev;
    bool edgcont;
    bool level0noise;
    bool level1noise;
    bool level2noise;
    bool level3noise;
    bool ccwcurve;
    bool opacityCurveBY;
    bool opacityCurveRG;
    bool opacityCurveW;
    bool opacityCurveWL;
    bool hhcurve;
    bool Chcurve;
    bool pastlev;
    bool satlev;
    bool wavclCurve;
    bool greenlow;
    bool bluelow;
    bool greenmed;
    bool bluemed;
    bool greenhigh;
    bool bluehigh;
    bool enacont;
    bool expcontrast;
    bool expchroma;
    bool expedge;
    bool expresid;
    bool expfinal;
    bool exptoning;
    bool expnoise;

    operator bool (void) const
    {
        bool retVal = enabled|strength|balance|iter|median|medianlev|linkedg|cbenab|lipst|Medgreinf|avoid|tmr|Lmethod|CHmethod|CHSLmethod|EDmethod|BAmethod|NPmethod|TMmethod|HSmethod|CLmethod|Backmethod|Tilesmethod|daubcoeffmethod|Dirmethod|rescon|resconH|reschro|tmrs|gamma|sup|sky|thres|threshold|threshold2|edgedetect|edgedetectthr|edgedetectthr2|edgesensi|edgeampli|chro|chroma|contrast|edgrad|edgval|edgthresh|thr|thrH|skinprotect|hueskin|hueskin2|hllev|bllev|edgcont|level0noise|level1noise|level2noise|level3noise|ccwcurve|opacityCurveBY|opacityCurveRG|opacityCurveW|opacityCurveWL|hhcurve|Chcurve|pastlev|satlev|wavclCurve|greenlow|bluelow|greenmed|bluemed|greenhigh|bluehigh|enacont|expcontrast|expchroma|expedge|expresid|expfinal|exptoning|expnoise;
        for (size_t i=0; i<9; ++i) {
            retVal |= c[i]|ch[i];
        }
        printf("bool: %d\n", retVal);
        return retVal;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class DirPyrEqualizerParamsEdited
{

public:
    bool enabled;
    bool gamutlab;
    bool mult[6];

    bool threshold;
    bool skinprotect;
    bool hueskin;
    //bool algo;

    operator bool (void) const
    {
        bool retVal = enabled|gamutlab|threshold|skinprotect|hueskin/*|algo*/;
        for (size_t i=0; i<6; ++i) {
            retVal |= mult[i];
        }
        printf("bool: %d\n", retVal);
        return retVal;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class HSVEqualizerParamsEdited
{

public:
    bool hcurve;
    bool scurve;
    bool vcurve;

    operator bool (void) const
    {
        printf("bool: %d\n", hcurve|scurve|vcurve);
        return hcurve|scurve|vcurve;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class FilmSimulationParamsEdited
{
public:
    bool enabled;
    bool clutFilename;
    bool strength;

    operator bool (void) const
    {
        printf("bool: %d\n", enabled|clutFilename|strength);
        return enabled|clutFilename|strength;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class RAWParamsEdited
{

public:
    class BayerSensor
    {

    public:
        bool method;
        bool ccSteps;
        bool exBlack0;
        bool exBlack1;
        bool exBlack2;
        bool exBlack3;
        bool exTwoGreen;
        bool dcbIterations;
        bool dcbEnhance;
        bool lmmseIterations;
        //bool allEnhance;
        bool greenEq;
        bool linenoise;

        operator bool (void) const
        {
            printf("bool: %d\n", method|ccSteps|exBlack0|exBlack1|exBlack2|exBlack3|exTwoGreen|dcbIterations|dcbEnhance|lmmseIterations/*|allEnhance*/|greenEq|linenoise);
            return method|ccSteps|exBlack0|exBlack1|exBlack2|exBlack3|exTwoGreen|dcbIterations|dcbEnhance|lmmseIterations/*|allEnhance*/|greenEq|linenoise;
        }
    };

    class XTransSensor
    {

    public:
        bool method;
        bool ccSteps;
        bool exBlackRed;
        bool exBlackGreen;
        bool exBlackBlue;

        operator bool (void) const
        {
            printf("bool: %d\n", method|ccSteps|exBlackRed|exBlackGreen|exBlackBlue);
            return method|ccSteps|exBlackRed|exBlackGreen|exBlackBlue;
        }
    };

    BayerSensor bayersensor;
    XTransSensor xtranssensor;

    bool caCorrection;
    bool caRed;
    bool caBlue;
    bool hotPixelFilter;
    bool deadPixelFilter;
    bool hotDeadPixelThresh;
    bool darkFrame;
    bool dfAuto;
    bool ff_file;
    bool ff_AutoSelect;
    bool ff_BlurRadius;
    bool ff_BlurType;
    bool ff_AutoClipControl;
    bool ff_clipControl;
    bool exPos;
    bool exPreser;

    operator bool (void) const
    {
        printf("bool: %d\n", bayersensor|xtranssensor|caCorrection|caRed|caBlue|hotPixelFilter|deadPixelFilter|hotDeadPixelThresh|darkFrame|dfAuto|ff_file|ff_AutoSelect|ff_BlurRadius|ff_BlurType|ff_AutoClipControl|ff_clipControl|exPos|exPreser);
        return bayersensor|xtranssensor|caCorrection|caRed|caBlue|hotPixelFilter|deadPixelFilter|hotDeadPixelThresh|darkFrame|dfAuto|ff_file|ff_AutoSelect|ff_BlurRadius|ff_BlurType|ff_AutoClipControl|ff_clipControl|exPos|exPreser;
    }
    void set(bool v);
    void initFrom (std::vector<const void*> elems);
    void combine (void* paramsToEdit, const void* newValues, bool dontforceSet);
};

class ParamsEdited
{

public:
    GeneralParamsEdited           general;
    ToneCurveParamsEdited         toneCurve;
    LCurveParamsEdited            labCurve;
    RGBCurvesParamsEdited         rgbCurves;
    ColorToningEdited             colorToning;
    SharpeningParamsEdited        sharpening;
    SharpeningParamsEdited        prsharpening;
    SharpenEdgeParamsEdited       sharpenEdge;
    SharpenMicroParamsEdited      sharpenMicro;
    VibranceParamsEdited          vibrance;
    ColorAppearanceParamsEdited   colorappearance;
    WBParamsEdited                wb;
    DefringeParamsEdited          defringe;
    DirPyrDenoiseParamsEdited     dirpyrDenoise;
    EPDParamsEdited               epd;
    ImpulseDenoiseParamsEdited    impulseDenoise;
    SHParamsEdited                sh;
    CropParamsEdited              crop;
    CoarseTransformParamsEdited   coarse;
    CommonTransformParamsEdited   commonTrans;
    RotateParamsEdited            rotate;
    DistortionParamsEdited        distortion;
    LensProfParamsEdited          lensProf;
    PerspectiveParamsEdited       perspective;
    GradientParamsEdited          gradient;
    PCVignetteParamsEdited        pcvignette;
    CACorrParamsEdited            cacorrection;
    VignettingParamsEdited        vignetting;
    ChannelMixerParamsEdited      chmixer;
    BlackWhiteParamsEdited        blackwhite;
    ResizeParamsEdited            resize;
    ColorManagementParamsEdited   icm;
    RAWParamsEdited               raw;
    DirPyrEqualizerParamsEdited   dirpyrequalizer;
    WaveletParamsEdited           wavelet;
    HSVEqualizerParamsEdited      hsvequalizer;
    FilmSimulationParamsEdited    filmSimulation;
    bool                          exif;
    bool                          iptc;

    ParamsEdited (bool value = false);

    void set   (bool v, int subPart = rtengine::procparams::ProcParams::FLAGS|rtengine::procparams::ProcParams::EXIF|rtengine::procparams::ProcParams::IPTC|rtengine::procparams::ProcParams::TOOL);

    void initFrom (const std::vector<rtengine::procparams::ProcParams>& src);
    void combine (rtengine::procparams::ProcParams& toEdit, const rtengine::procparams::ProcParams& mods, bool forceSet);

    bool isTagsSet();
    bool isToolSet();
    bool isExifSet();
    bool isIptcSet();

    bool operator== (const ParamsEdited& other);
    bool operator!= (const ParamsEdited& other);
};
#endif
