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
#include "paramsedited.h"
#include <cstring>
#include "options.h"
#include "addsetids.h"

using namespace rtengine;
using namespace rtengine::procparams;


void GeneralParamsEdited::set (bool v)
{
    rank         = v;
    colorlabel   = v;
    intrash      = v;
}

void ToneCurveParamsEdited::set (bool v)
{
    curve      = v;
    curve2     = v;
    curveMode  = v;
    curveMode2 = v;
    brightness = v;
    black      = v;
    contrast   = v;
    saturation = v;
    shcompr    = v;
    hlcompr    = v;
    hlcomprthresh = v;
    autoexp    = v;
    clip       = v;
    expcomp    = v;
    hrenabled   = v;
    method    = v;
}

void LCurveParamsEdited::set (bool v)
{
    lcurve      = v;
    acurve      = v;
    bcurve      = v;
    cccurve     = v;
    chcurve     = v;
    lhcurve     = v;
    hhcurve     = v;
    lccurve    = v;
    clcurve    = v;
    brightness  = v;
    contrast    = v;
    chromaticity    = v;
    avoidcolorshift = v;
    rstprotection   = v;
    lcredsk         = v;
    enabled         = v;
    method          = v;
}

void RGBCurvesParamsEdited::set (bool v)
{
    lumamode       = v;
    rcurve         = v;
    gcurve         = v;
    bcurve         = v;
}

void ColorToningEdited::set (bool v)
{
    enabled      = v;
    autosat      = v;
    opacityCurve = v;
    colorCurve   = v;
    satprotectionthreshold = v;
    saturatedopacity       = v;
    strength               = v;
    shadowsColSat          = v;
    hlColSat   = v;
    balance    = v;
    clcurve    = v;
    method     = v;
    twocolor   = v;
    cl2curve   = v;
    redlow     = v;
    greenlow   = v;
    bluelow    = v;
    satlow     = v;
    sathigh    = v;
    redmed     = v;
    greenmed   = v;
    bluemed    = v;
    redhigh    = v;
    greenhigh  = v;
    bluehigh   = v;
    lumamode   = v;
}

void SharpeningParamsEdited::set (bool v)
{
    enabled            = v;
    radius             = v;
    amount             = v;
    threshold          = v;
    edgesonly          = v;
    edges_radius       = v;
    edges_tolerance    = v;
    halocontrol        = v;
    halocontrol_amount = v;
    method         = v;
    deconvamount   = v;
    deconvradius   = v;
    deconviter     = v;
    deconvdamping  = v;
}

void SharpenEdgeParamsEdited::set (bool v)
{
    enabled       = v;
    passes        = v;
    amount        = v;
    threechannels = v;
}

void SharpenMicroParamsEdited::set (bool v)
{
    enabled      = v;
    matrix       = v;
    amount       = v;
    uniformity   = v;
}

void VibranceParamsEdited::set (bool v)
{
    enabled          = v;
    pastels          = v;
    saturated        = v;
    psthreshold      = v;
    protectskins     = v;
    avoidcolorshift  = v;
    pastsattog       = v;
    skintonescurve   = v;
}

void ColorAppearanceParamsEdited::set (bool v)
{
    enabled    = v;
    degree     = v;
    autodegree = v;
    surround     = v;
    adapscen    = v;
    autoadapscen = v;
    adaplum    = v;
    badpixsl    = v;
    wbmodel    = v;
    algo    = v;
    jlight     = v;
    qbright     = v;
    chroma     = v;
    schroma     = v;
    mchroma     = v;
    contrast     = v;
    qcontrast     = v;
    colorh     = v;
    rstprotection     = v;
    surrsource = v;
    gamut = v;
    //badpix = v;
    datacie = v;
    tonecie = v;
    //sharpcie = v;
    curve      = v;
    curve2     = v;
    curve3     = v;
    curveMode  = v;
    curveMode2 = v;
    curveMode3 = v;
}

void WBParamsEdited::set (bool v)
{
    method                  = v;
    green                   = v;
    temperature             = v;
    equal                   = v;
}

void DefringeParamsEdited::set (bool v)
{
    enabled           = v;
    radius            = v;
    threshold         = v;
    huecurve          = v;
}

void DirPyrDenoiseParamsEdited::set (bool v)
{
    enabled      = v;
    enhance      = v;
    //perform      = v;
    lcurve      = v;
    cccurve      = v;
    median      = v;
    autochroma      = v;
    luma         = v;
    Ldetail      = v;
    chroma       = v;
    redchro      = v;
    bluechro     = v;
    gamma        = v;
    passes        = v;
    dmethod      = v;
    Lmethod      = v;
    Cmethod      = v;
    C2method      = v;
    smethod      = v;
    medmethod      = v;
    methodmed      = v;
    rgbmethod      = v;
}

void EPDParamsEdited::set (bool v)
{
    enabled                = v;
    strength            = v;
    gamma            = v;
    edgeStopping        = v;
    scale               = v;
    reweightingIterates = v;
}

void ImpulseDenoiseParamsEdited::set (bool v)
{
    enabled     = v;
    thresh      = v;
}

void SHParamsEdited::set (bool v)
{
    enabled       = v;
    hq            = v;
    highlights    = v;
    htonalwidth   = v;
    shadows       = v;
    stonalwidth   = v;
    localcontrast = v;
    radius        = v;
}

void CropParamsEdited::set (bool v)
{
    enabled = v;
    x       = v;
    y       = v;
    w       = v;
    h       = v;
    fixratio = v;
    ratio   = v;
    orientation = v;
    guide   = v;
}

void CoarseTransformParamsEdited::set (bool v)
{
    rotate = v;
    hflip = v;
    vflip = v;
}

void CommonTransformParamsEdited::set (bool v)
{
    autofill = v;
}

void RotateParamsEdited::set (bool v)
{
    degree = v;
}

void DistortionParamsEdited::set (bool v)
{
    amount = v;
}

void LensProfParamsEdited::set (bool v)
{
    lcpFile = v;
    useDist = v;
    useVign = v;
    useCA = v;
}

void PerspectiveParamsEdited::set (bool v)
{
    horizontal = v;
    vertical = v;
}

void GradientParamsEdited::set (bool v)
{
    enabled = v;
    degree = v;
    feather = v;
    strength = v;
    centerX = v;
    centerY = v;
}

void PCVignetteParamsEdited::set (bool v)
{
    enabled = v;
    strength = v;
    feather = v;
    roundness = v;
}

void CACorrParamsEdited::set (bool v)
{
    red = v;
    blue = v;
}

void VignettingParamsEdited::set (bool v)
{
    amount = v;
    radius = v;
    strength = v;
    centerX = v;
    centerY = v;
}

void ChannelMixerParamsEdited::set (bool v)
{
    red[0] = v;
    red[1] = v;
    red[2] = v;
    green[0] = v;
    green[1] = v;
    green[2] = v;
    blue[0] = v;
    blue[1] = v;
    blue[2] = v;
}

void BlackWhiteParamsEdited::set (bool v)
{
    enabled   = v;
    enabledcc   = v;
    mixerRed   = v;
    mixerOrange   = v;
    mixerYellow   = v;
    mixerGreen   = v;
    mixerCyan   = v;
    mixerBlue   = v;
    mixerMagenta   = v;
    mixerPurple   = v;
    gammaRed   = v;
    gammaGreen   = v;
    gammaBlue   = v;
    filter   = v;
    setting   = v;
    method   = v;
    luminanceCurve = v;
    beforeCurve      = v;
    beforeCurveMode  = v;
    afterCurve      = v;
    afterCurveMode  = v;
    autoc    = v;
    algo    = v;
}

void ResizeParamsEdited::set (bool v)
{
    scale     = v;
    appliesTo = v;
    method    = v;
    dataspec  = v;
    width     = v;
    height    = v;
    enabled   = v;
}

void ColorManagementParamsEdited::set (bool v)
{
    input          = v;
    toneCurve      = v;
    applyLookTable = v;
    applyBaselineExposureOffset = v;
    applyHueSatMap = v;
    blendCMSMatrix = v;
    dcpIlluminant  = v;
    working        = v;
    output         = v;
    gamma          = v;
    gamfree        = v;
    freegamma      = v;
    gampos         = v;
    slpos          = v;
}

void RAWParamsEdited::set (bool v)
{
    bayersensor.method = v;
    bayersensor.ccSteps = v;
    bayersensor.exBlack0 = v;
    bayersensor.exBlack1 = v;
    bayersensor.exBlack2 = v;
    bayersensor.exBlack3 = v;
    bayersensor.exTwoGreen = v;
    bayersensor.dcbIterations = v;
    bayersensor.dcbEnhance = v;
    //bayersensor.allEnhance = v;
    bayersensor.lmmseIterations = v;
    bayersensor.greenEq = v;
    bayersensor.linenoise = v;
    xtranssensor.method = v;
    xtranssensor.ccSteps = v;
    xtranssensor.exBlackRed = v;
    xtranssensor.exBlackGreen = v;
    xtranssensor.exBlackBlue = v;
    caCorrection = v;
    caBlue  = v;
    caRed   = v;
    hotPixelFilter = v;
    deadPixelFilter = v;
    hotDeadPixelThresh = v;
    darkFrame = v;
    dfAuto = v;
    ff_file = v;
    ff_AutoSelect = v;
    ff_BlurRadius = v;
    ff_BlurType = v;
    ff_AutoClipControl = v;
    ff_clipControl = v;
    exPos = v;
    exPreser = v;
}

void DirPyrEqualizerParamsEdited::set (bool v)
{
    enabled = v;
    gamutlab = v;

    for(int i = 0; i < 6; i++) {
        mult[i] = v;
    }

    threshold = v;
    skinprotect = v;
    hueskin = v;
    //algo = v;
}

void WaveletParamsEdited::set (bool v)
{
    enabled = v;
    strength = v;
    balance = v;
    iter = v;
    median = v;
    medianlev = v;
    linkedg = v;
    cbenab = v;
    enacont = v;
    greenhigh = v;
    greenmed = v;
    greenlow = v;
    bluehigh = v;
    bluemed = v;
    bluelow = v;
    lipst = v;
    Medgreinf = v;
    avoid = v;
    tmr = v;
    Lmethod = v;
    CLmethod = v;
    Backmethod = v;
    Tilesmethod = v;
    daubcoeffmethod = v;
    CHmethod = v;
    CHSLmethod = v;
    EDmethod = v;
    NPmethod = v;
    BAmethod = v;
    TMmethod = v;
    HSmethod = v;
    Dirmethod = v;
    rescon = v;
    resconH = v;
    reschro = v;
    tmrs = v;
    gamma = v;
    sup = v;
    sky = v;
    thres = v;
    threshold = v;
    threshold2 = v;
    edgedetect = v;
    edgedetectthr = v;
    edgedetectthr2 = v;
    edgesensi = v;
    edgeampli = v;
    chroma = v;
    chro = v;
    contrast = v;
    edgrad = v;
    edgval = v;
    edgthresh = v;
    thr = v;
    thrH = v;
    skinprotect = v;
    hueskin = v;
    hueskin2 = v;
    hllev = v;
    bllev = v;
    edgcont = v;
    level0noise = v;
    level1noise = v;
    level2noise = v;
    level3noise = v;
    ccwcurve = v;
    opacityCurveRG   = v;
    opacityCurveBY   = v;
    opacityCurveW   = v;
    opacityCurveWL   = v;
    hhcurve     = v;
    Chcurve     = v;
    wavclCurve     = v;
    pastlev = v;
    satlev = v;

    //enacont = v;
    //enachrom = v;
    //enaedge = v;
    //enares = v;

    expfinal = v;
    expcontrast = v;
    expchroma = v;
    expedge = v;
    expresid = v;
    exptoning = v;
    expnoise = v;

    for(int i = 0; i < 9; i++) {
        c[i] = v;
    }

    for(int i = 0; i < 9; i++) {
        ch[i] = v;
    }
}

void HSVEqualizerParamsEdited::set (bool v)
{
    hcurve = v;
    scurve = v;
    vcurve = v;
}

void FilmSimulationParamsEdited::set (bool v)
{
    enabled = v;
    clutFilename = v;
    strength = v;
}

ParamsEdited::ParamsEdited (bool value)
{
    set (value);
}

// By default set all subparts to 'v'
void ParamsEdited::set (bool v, int subPart)
{
    if (subPart & ProcParams::FLAGS) {
        general.set(v);
    }
    if (subPart & ProcParams::TOOL) {
        toneCurve.set(v);
        labCurve.set(v);
        rgbCurves.set(v);
        colorToning.set(v);
        sharpening.set(v);
        prsharpening.set(v);
        sharpenEdge.set(v);
        sharpenMicro.set(v);
        vibrance.set(v);
        colorappearance.set(v);
        wb.set(v);
        defringe.set(v);
        impulseDenoise.set(v);
        dirpyrDenoise.set(v);
        epd.set(v);
        sh.set(v);
        crop.set(v);
        coarse.set(v);
        commonTrans.set(v);
        rotate.set(v);
        distortion.set(v);
        lensProf.set(v);
        perspective.set(v);
        gradient.set(v);
        pcvignette.set(v);
        cacorrection.set(v);
        vignetting.set(v);
        chmixer.set(v);
        blackwhite.set(v);
        resize.set(v);
        icm.set(v);
        raw.set(v);
        wavelet.set(v);
        dirpyrequalizer.set(v);
        hsvequalizer.set(v);
        filmSimulation.set(v);
    }
    if (subPart & ProcParams::EXIF) {
        exif = v;
    }
    if (subPart & ProcParams::IPTC) {
        iptc = v;
    }
}

void ToneCurveParamsEdited::initFrom (std::vector<const void*> elems)
{
    const ToneCurveParams *e0 = static_cast< const ToneCurveParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const ToneCurveParams* e = static_cast< const ToneCurveParams* >(elems.at(i));
        curve = curve && e0->curve == e->curve;
        curve2 = curve2 && e0->curve2 == e->curve2;
        curveMode = curveMode && e0->curveMode == e->curveMode;
        curveMode2 = curveMode2 && e0->curveMode2 == e->curveMode2;
        brightness = brightness && e0->brightness == e->brightness;
        black = black && e0->black == e->black;
        contrast = contrast && e0->contrast == e->contrast;
        saturation = saturation && e0->saturation == e->saturation;
        shcompr = shcompr && e0->shcompr == e->shcompr;
        hlcompr = hlcompr && e0->hlcompr == e->hlcompr;
        hlcomprthresh = hlcomprthresh && e0->hlcomprthresh == e->hlcomprthresh;
        autoexp = autoexp && e0->autoexp == e->autoexp;
        clip = clip && e0->clip == e->clip;
        expcomp = expcomp && e0->expcomp == e->expcomp;
        hrenabled = hrenabled && e0->hrenabled == e->hrenabled;
        method = method && e0->method == e->method;
    }
}

void LCurveParamsEdited::initFrom (std::vector<const void*> elems)
{
    const LCurveParams *e0 = static_cast< const LCurveParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const LCurveParams* e = static_cast< const LCurveParams* >(elems.at(i));
        lcurve = lcurve && e0->lcurve == e->lcurve;
        acurve = acurve && e0->acurve == e->acurve;
        bcurve = bcurve && e0->bcurve == e->bcurve;
        cccurve = cccurve && e0->cccurve == e->cccurve;
        chcurve = chcurve && e0->chcurve == e->chcurve;
        lhcurve = lhcurve && e0->lhcurve == e->lhcurve;
        hhcurve = hhcurve && e0->hhcurve == e->hhcurve;
        lccurve = lccurve && e0->lccurve == e->lccurve;
        clcurve = clcurve && e0->clcurve == e->clcurve;
        brightness = brightness && e0->brightness == e->brightness;
        contrast = contrast && e0->contrast == e->contrast;
        chromaticity = chromaticity && e0->chromaticity == e->chromaticity;
        avoidcolorshift = avoidcolorshift && e0->avoidcolorshift == e->avoidcolorshift;
        rstprotection = rstprotection && e0->rstprotection == e->rstprotection;
        lcredsk = lcredsk && e0->lcredsk == e->lcredsk;
    }
}

void RGBCurvesParamsEdited::initFrom (std::vector<const void*> elems)
{
    const RGBCurvesParams *e0 = static_cast< const RGBCurvesParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const RGBCurvesParams* e = static_cast< const RGBCurvesParams* >(elems.at(i));
        lumamode = lumamode && e0->lumamode == e->lumamode;
        rcurve = rcurve && e0->rcurve == e->rcurve;
        gcurve = gcurve && e0->gcurve == e->gcurve;
        bcurve = bcurve && e0->bcurve == e->bcurve;
    }
}

void ColorToningEdited::initFrom (std::vector<const void*> elems)
{
    const ColorToningParams *e0 = static_cast< const ColorToningParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const ColorToningParams *e = static_cast< const ColorToningParams* >(elems.at(i));
        enabled = enabled && e0->enabled == e->enabled;
        twocolor = twocolor && e0->twocolor == e->twocolor;
        opacityCurve = opacityCurve && e0->opacityCurve == e->opacityCurve;
        colorCurve = colorCurve && e0->colorCurve == e->colorCurve;
        autosat = autosat && e0->autosat == e->autosat;
        satprotectionthreshold = satprotectionthreshold && e0->satProtectionThreshold == e->satProtectionThreshold;
        saturatedopacity = saturatedopacity && e0->saturatedOpacity == e->saturatedOpacity;
        strength = strength && e0->strength == e->strength;
        shadowsColSat = shadowsColSat && e0->shadowsColSat == e->shadowsColSat;
        hlColSat = hlColSat && e0->hlColSat == e->hlColSat;
        balance = balance && e0->balance == e->balance;
        clcurve = clcurve && e0->clcurve == e->clcurve;
        cl2curve = cl2curve && e0->cl2curve == e->cl2curve;
        method = method && e0->method == e->method;
        redlow = redlow && e0->redlow == e->redlow;
        greenlow = greenlow && e0->greenlow == e->greenlow;
        bluelow = bluelow && e0->bluelow == e->bluelow;
        satlow = satlow && e0->satlow == e->satlow;
        sathigh = sathigh && e0->sathigh == e->sathigh;
        redmed = redmed && e0->redmed == e->redmed;
        greenmed = greenmed && e0->greenmed == e->greenmed;
        bluemed = bluemed && e0->bluemed == e->bluemed;
        redhigh = redhigh && e0->redhigh == e->redhigh;
        greenhigh = greenhigh && e0->greenhigh == e->greenhigh;
        bluehigh = bluehigh && e0->bluehigh == e->bluehigh;
        lumamode = lumamode && e0->lumamode == e->lumamode;
    }
}

void SharpeningParamsEdited::initFrom (std::vector<const void*> elems)
{
    const SharpeningParams *e0 = static_cast< const SharpeningParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const SharpeningParams *e = static_cast< const SharpeningParams* >(elems.at(i));
        enabled = enabled && e0->enabled == e->enabled;
        radius = radius && e0->radius == e->radius;
        amount = amount && e0->amount == e->amount;
        threshold = threshold && e0->threshold == e->threshold;
        edgesonly = edgesonly && e0->edgesonly == e->edgesonly;
        edges_radius = edges_radius && e0->edges_radius == e->edges_radius;
        edges_tolerance = edges_tolerance && e0->edges_tolerance == e->edges_tolerance;
        halocontrol = halocontrol && e0->halocontrol == e->halocontrol;
        halocontrol_amount = halocontrol_amount && e0->halocontrol_amount == e->halocontrol_amount;
        method = method && e0->method == e->method;
        deconvamount = deconvamount && e0->deconvamount == e->deconvamount;
        deconvradius = deconvradius && e0->deconvradius == e->deconvradius;
        deconviter = deconviter && e0->deconviter == e->deconviter;
        deconvdamping = deconvdamping && e0->deconvdamping == e->deconvdamping;
    }
}

void SharpenEdgeParamsEdited::initFrom (std::vector<const void*> elems)
{
    const SharpenEdgeParams *e0 = static_cast< const SharpenEdgeParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const SharpenEdgeParams *e = static_cast< const SharpenEdgeParams* >(elems.at(i));
        enabled = enabled && e0->enabled == e->enabled;
        passes = passes && e0->passes == e->passes;
        amount = amount && e0->amount == e->amount;
        threechannels = threechannels && e0->threechannels == e->threechannels;
    }
}

void SharpenMicroParamsEdited::initFrom (std::vector<const void*> elems)
{
    const SharpenMicroParams *e0 = static_cast< const SharpenMicroParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const SharpenMicroParams *e = static_cast< const SharpenMicroParams* >(elems.at(i));
        enabled = enabled && e0->enabled == e->enabled;
        matrix = matrix && e0->matrix == e->matrix;
        amount = amount && e0->amount == e->amount;
        uniformity = uniformity && e0->uniformity == e->uniformity;
    }
}

void VibranceParamsEdited::initFrom (std::vector<const void*> elems)
{
    const VibranceParams *e0 = static_cast< const VibranceParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const VibranceParams *e = static_cast< const VibranceParams* >(elems.at(i));
        enabled = enabled && e0->enabled == e->enabled;
        pastels = pastels && e0->pastels == e->pastels;
        saturated = saturated && e0->saturated == e->saturated;
        psthreshold = psthreshold && e0->psthreshold == e->psthreshold;
        protectskins = protectskins && e0->protectskins == e->protectskins;
        avoidcolorshift = avoidcolorshift && e0->avoidcolorshift == e->avoidcolorshift;
        pastsattog = pastsattog && e0->pastsattog == e->pastsattog;
        skintonescurve = skintonescurve && e0->skintonescurve == e->skintonescurve;
    }
}

void ColorAppearanceParamsEdited::initFrom (std::vector<const void*> elems)
{
    const ColorAppearanceParams *e0 = static_cast< const ColorAppearanceParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const ColorAppearanceParams *e = static_cast< const ColorAppearanceParams* >(elems.at(i));
        enabled = enabled && e0->enabled == e->enabled;
        degree = degree && e0->degree == e->degree;
        autodegree = autodegree && e0->autodegree == e->autodegree;
        surround = surround && e0->surround == e->surround;
        adapscen = adapscen && e0->adapscen == e->adapscen;
        autoadapscen = autoadapscen && e0->autoadapscen == e->autoadapscen;
        adaplum = adaplum && e0->adaplum == e->adaplum;
        badpixsl = badpixsl && e0->badpixsl == e->badpixsl;
        wbmodel = wbmodel && e0->wbmodel == e->wbmodel;
        algo = algo && e0->algo == e->algo;
        jlight = jlight && e0->jlight == e->jlight;
        qbright = qbright && e0->qbright == e->qbright;
        chroma = chroma && e0->chroma == e->chroma;
        schroma = schroma && e0->schroma == e->schroma;
        mchroma = mchroma && e0->mchroma == e->mchroma;
        rstprotection = rstprotection && e0->rstprotection == e->rstprotection;
        contrast = contrast && e0->contrast == e->contrast;
        qcontrast = qcontrast && e0->qcontrast == e->qcontrast;
        colorh = colorh && e0->colorh == e->colorh;
        surrsource = surrsource && e0->surrsource == e->surrsource;
        gamut = gamut && e0->gamut == e->gamut;
        //badpix = badpix && e0->badpix == e->badpix;
        datacie = datacie && e0->datacie == e->datacie;
        tonecie = tonecie && e0->tonecie == e->tonecie;
        //sharpcie = sharpcie && e0->sharpcie == e->sharpcie;
        curve = curve && e0->curve == e->curve;
        curve3 = curve3 && e0->curve3 == e->curve3;
        curve2 = curve2 && e0->curve2 == e->curve2;
        curveMode = curveMode && e0->curveMode == e->curveMode;
        curveMode2 = curveMode2 && e0->curveMode2 == e->curveMode2;
        curveMode3 = curveMode3 && e0->curveMode3 == e->curveMode3;
    }
}

void WBParamsEdited::initFrom (std::vector<const void*> elems)
{
    const WBParams *e0 = static_cast< const WBParams* >(elems.at(0));
    size_t size = elems.size();
  for (size_t i = 1; i < size; ++i) {
    const WBParams *e = static_cast< const WBParams* >(elems.at(i));
      method = method && e0->method == e->method;
        green = green && e0->green == e->green;
        equal = equal && e0->equal == e->equal;
        temperature = temperature && e0->temperature == e->temperature;
    }
}

void DefringeParamsEdited::initFrom (std::vector<const void*> elems)
{
    const DefringeParams *e0 = static_cast< const DefringeParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const DefringeParams *e = static_cast< const DefringeParams* >(elems.at(i));
        enabled = enabled && e0->enabled == e->enabled;
        radius = radius && e0->radius == e->radius;
        threshold = threshold && e0->threshold == e->threshold;
        huecurve = huecurve && e0->huecurve == e->huecurve;
    }
}

void DirPyrDenoiseParamsEdited::initFrom (std::vector<const void*> elems)
{
    const DirPyrDenoiseParams *e0 = static_cast< const DirPyrDenoiseParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const DirPyrDenoiseParams *e = static_cast< const DirPyrDenoiseParams* >(elems.at(i));
        enabled = enabled && e0->enabled == e->enabled;
        enhance = enhance && e0->enhance == e->enhance;
        median = median && e0->median == e->median;
        autochroma = autochroma && e0->autochroma == e->autochroma;
        //perform = perform && e0->perform == e->perform;
        luma = luma && e0->luma == e->luma;
        lcurve = lcurve && e0->lcurve == e->lcurve;
        cccurve = cccurve && e0->cccurve == e->cccurve;
        Ldetail = Ldetail && e0->Ldetail == e->Ldetail;
        chroma = chroma && e0->chroma == e->chroma;
        redchro = redchro && e0->redchro == e->redchro;
        bluechro = bluechro && e0->bluechro == e->bluechro;
        gamma = gamma && e0->gamma == e->gamma;
        passes = passes && e0->passes == e->passes;
        dmethod = dmethod && e0->dmethod == e->dmethod;
        Lmethod = Lmethod && e0->Lmethod == e->Lmethod;
        Cmethod = Cmethod && e0->Cmethod == e->Cmethod;
        C2method = C2method && e0->C2method == e->C2method;
        smethod = smethod && e0->smethod == e->smethod;
        medmethod = medmethod && e0->medmethod == e->medmethod;
        methodmed = methodmed && e0->methodmed == e->methodmed;
        rgbmethod = rgbmethod && e0->rgbmethod == e->rgbmethod;
    }
}

void EPDParamsEdited::initFrom (std::vector<const void*> elems)
{
    const EPDParams *e0 = static_cast< const EPDParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const EPDParams *e = static_cast< const EPDParams* >(elems.at(i));
        enabled = enabled && e0->enabled == e->enabled;
        strength = strength && e0->strength == e->strength;
        gamma = gamma && e0->gamma == e->gamma;
        edgeStopping = edgeStopping && e0->edgeStopping == e->edgeStopping;
        scale = scale && e0->scale == e->scale;
        reweightingIterates = reweightingIterates && e0->reweightingIterates == e->reweightingIterates;
    }
}

void ImpulseDenoiseParamsEdited::initFrom (std::vector<const void*> elems)
{
    const ImpulseDenoiseParams *e0 = static_cast< const ImpulseDenoiseParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const ImpulseDenoiseParams *e = static_cast< const ImpulseDenoiseParams* >(elems.at(i));
        enabled = enabled && e0->enabled == e->enabled;
        thresh = thresh && e0->thresh == e->thresh;
    }
}

void SHParamsEdited::initFrom (std::vector<const void*> elems)
{
    const SHParams *e0 = static_cast< const SHParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const SHParams *e = static_cast< const SHParams* >(elems.at(i));
        enabled = enabled && e0->enabled == e->enabled;
        hq = hq && e0->hq == e->hq;
        highlights = highlights && e0->highlights == e->highlights;
        htonalwidth = htonalwidth && e0->htonalwidth == e->htonalwidth;
        shadows = shadows && e0->shadows == e->shadows;
        stonalwidth = stonalwidth && e0->stonalwidth == e->stonalwidth;
        localcontrast = localcontrast && e0->localcontrast == e->localcontrast;
        radius = radius && e0->radius == e->radius;
    }
}

void CropParamsEdited::initFrom (std::vector<const void*> elems)
{
    const CropParams *e0 = static_cast< const CropParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const CropParams *e = static_cast< const CropParams* >(elems.at(i));
        enabled = enabled && e0->enabled == e->enabled;
        x = x && e0->x == e->x;
        y = y && e0->y == e->y;
        w = w && e0->w == e->w;
        h = h && e0->h == e->h;
        fixratio = fixratio && e0->fixratio == e->fixratio;
        ratio = ratio && e0->ratio == e->ratio;
        orientation = orientation && e0->orientation == e->orientation;
        guide = guide && e0->guide == e->guide;
    }
}

void CoarseTransformParamsEdited::initFrom (std::vector<const void*> elems)
{
    const CoarseTransformParams *e0 = static_cast< const CoarseTransformParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const CoarseTransformParams *e = static_cast< const CoarseTransformParams* >(elems.at(i));
        rotate = rotate && e0->rotate == e->rotate;
        hflip = hflip && e0->hflip == e->hflip;
        vflip = vflip && e0->vflip == e->vflip;
    }
}

void CommonTransformParamsEdited::initFrom (std::vector<const void*> elems)
{
    const CommonTransformParams *e0 = static_cast< const CommonTransformParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const CommonTransformParams *e = static_cast< const CommonTransformParams* >(elems.at(i));
        autofill = autofill && e0->autofill == e->autofill;
    }
}

void RotateParamsEdited::initFrom (std::vector<const void*> elems)
{
    const RotateParams *e0 = static_cast< const RotateParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const RotateParams *e = static_cast< const RotateParams* >(elems.at(i));
        degree = degree && e0->degree == e->degree;
    }
}

void DistortionParamsEdited::initFrom (std::vector<const void*> elems)
{
    const DistortionParams *e0 = static_cast< const DistortionParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const DistortionParams *e = static_cast< const DistortionParams* >(elems.at(i));
        amount = amount && e0->amount == e->amount;
    }
}

void LensProfParamsEdited::initFrom (std::vector<const void*> elems)
{
    const LensProfParams *e0 = static_cast< const LensProfParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const LensProfParams *e = static_cast< const LensProfParams* >(elems.at(i));
        lcpFile = lcpFile && e0->lcpFile == e->lcpFile;
        useDist = useDist && e0->useDist == e->useDist;
        useVign = useVign && e0->useVign == e->useVign;
        useCA = useCA && e0->useCA == e->useCA;
    }
}

void PerspectiveParamsEdited::initFrom (std::vector<const void*> elems)
{
    const PerspectiveParams *e0 = static_cast< const PerspectiveParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const PerspectiveParams *e = static_cast< const PerspectiveParams* >(elems.at(i));
        horizontal = horizontal && e0->horizontal == e->horizontal;
        vertical = vertical && e0->vertical == e->vertical;
    }
}

void GradientParamsEdited::initFrom (std::vector<const void*> elems)
{
    const GradientParams *e0 = static_cast< const GradientParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const GradientParams *e = static_cast< const GradientParams* >(elems.at(i));
        enabled = enabled && e0->enabled == e->enabled;
        degree = degree && e0->degree == e->degree;
        feather = feather && e0->feather == e->feather;
        strength = strength && e0->strength == e->strength;
        centerX = centerX && e0->centerX == e->centerX;
        centerY = centerY && e0->centerY == e->centerY;
    }
}

void PCVignetteParamsEdited::initFrom (std::vector<const void*> elems)
{
    const PCVignetteParams *e0 = static_cast< const PCVignetteParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const PCVignetteParams *e = static_cast< const PCVignetteParams* >(elems.at(i));
        enabled = enabled && e0->enabled == e->enabled;
        strength = strength && e0->strength == e->strength;
        feather = feather && e0->feather == e->feather;
        roundness = roundness && e0->roundness == e->roundness;
    }
}

void CACorrParamsEdited::initFrom (std::vector<const void*> elems)
{
    const CACorrParams *e0 = static_cast< const CACorrParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const CACorrParams *e = static_cast< const CACorrParams* >(elems.at(i));
        red = red && e0->red == e->red;
        blue = blue && e0->blue == e->blue;
    }
}

void VignettingParamsEdited::initFrom (std::vector<const void*> elems)
{
    const VignettingParams *e0 = static_cast< const VignettingParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const VignettingParams *e = static_cast< const VignettingParams* >(elems.at(i));
        amount = amount && e0->amount == e->amount;
        radius = radius && e0->radius == e->radius;
        strength = strength && e0->strength == e->strength;
        centerX = centerX && e0->centerX == e->centerX;
        centerY = centerY && e0->centerY == e->centerY;
    }
}

void ChannelMixerParamsEdited::initFrom (std::vector<const void*> elems)
{
    const ChannelMixerParams *e0 = static_cast< const ChannelMixerParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const ChannelMixerParams *e = static_cast< const ChannelMixerParams* >(elems.at(i));
        red[0] = red[0] && e0->red[0] == e->red[0];
        red[1] = red[1] && e0->red[1] == e->red[1];
        red[2] = red[2] && e0->red[2] == e->red[2];
        green[0] = green[0] && e0->green[0] == e->green[0];
        green[1] = green[1] && e0->green[1] == e->green[1];
        green[2] = green[2] && e0->green[2] == e->green[2];
        blue[0] = blue[0] && e0->blue[0] == e->blue[0];
        blue[1] = blue[1] && e0->blue[1] == e->blue[1];
        blue[2] = blue[2] && e0->blue[2] == e->blue[2];
    }
}

void BlackWhiteParamsEdited::initFrom (std::vector<const void*> elems)
{
    const BlackWhiteParams *e0 = static_cast< const BlackWhiteParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const BlackWhiteParams *e = static_cast< const BlackWhiteParams* >(elems.at(i));
        enabledcc = enabledcc && e0->enabledcc == e->enabledcc;
        enabled = enabled && e0->enabled == e->enabled;
        mixerRed = mixerRed && e0->mixerRed == e->mixerRed;
        mixerOrange = mixerOrange && e0->mixerOrange == e->mixerOrange;
        mixerYellow = mixerYellow && e0->mixerYellow == e->mixerYellow;
        mixerGreen = mixerGreen && e0->mixerGreen == e->mixerGreen;
        mixerCyan = mixerCyan && e0->mixerCyan == e->mixerCyan;
        mixerBlue = mixerBlue && e0->mixerBlue == e->mixerBlue;
        mixerMagenta = mixerMagenta && e0->mixerMagenta == e->mixerMagenta;
        mixerPurple = mixerPurple && e0->mixerPurple == e->mixerPurple;
        gammaRed = gammaRed && e0->gammaRed == e->gammaRed;
        gammaGreen = gammaGreen && e0->gammaGreen == e->gammaGreen;
        gammaBlue = gammaBlue && e0->gammaBlue == e->gammaBlue;
        filter = filter && e0->filter == e->filter;
        setting = setting && e0->setting == e->setting;
        luminanceCurve = luminanceCurve && e0->luminanceCurve == e->luminanceCurve;
        method = method && e0->method == e->method;
        beforeCurve = beforeCurve && e0->beforeCurve == e->beforeCurve;
        beforeCurveMode = beforeCurveMode && e0->beforeCurveMode == e->beforeCurveMode;
        afterCurve = afterCurve && e0->afterCurve == e->afterCurve;
        afterCurveMode = afterCurveMode && e0->afterCurveMode == e->afterCurveMode;
        autoc = autoc && e0->autoc == e->autoc;
        algo = algo && e0->algo == e->algo;
    }
}

void ResizeParamsEdited::initFrom (std::vector<const void*> elems)
{
    const ResizeParams *e0 = static_cast< const ResizeParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const ResizeParams *e = static_cast< const ResizeParams* >(elems.at(i));
        scale = scale && e0->scale == e->scale;
        appliesTo = appliesTo && e0->appliesTo == e->appliesTo;
        method = method && e0->method == e->method;
        dataspec = dataspec && e0->dataspec == e->dataspec;
        width = width && e0->width == e->width;
        height = height && e0->height == e->height;
        enabled = enabled && e0->enabled == e->enabled;
    }
}

void ColorManagementParamsEdited::initFrom (std::vector<const void*> elems)
{
    const ColorManagementParams *e0 = static_cast< const ColorManagementParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const ColorManagementParams *e = static_cast< const ColorManagementParams* >(elems.at(i));
        input = input && e0->input == e->input;
        toneCurve = toneCurve && e0->toneCurve == e->toneCurve;
        applyLookTable = applyLookTable && e0->applyLookTable == e->applyLookTable;
        applyBaselineExposureOffset = applyBaselineExposureOffset && e0->applyBaselineExposureOffset == e->applyBaselineExposureOffset;
        applyHueSatMap = applyHueSatMap && e0->applyHueSatMap == e->applyHueSatMap;
        blendCMSMatrix = blendCMSMatrix && e0->blendCMSMatrix == e->blendCMSMatrix;
        dcpIlluminant = dcpIlluminant && e0->dcpIlluminant == e->dcpIlluminant;
        working = working && e0->working == e->working;
        output = output && e0->output == e->output;
        gamma = gamma && e0->gamma == e->gamma;
        freegamma = freegamma && e0->freegamma == e->freegamma;
        gampos = gampos && e0->gampos == e->gampos;
        slpos = slpos && e0->slpos == e->slpos;
    }
}

void RAWParamsEdited::initFrom (std::vector<const void*> elems)
{
    const RAWParams *e0 = static_cast< const RAWParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const RAWParams *e = static_cast< const RAWParams* >(elems.at(i));
        bayersensor.method = bayersensor.method && e0->bayersensor.method == e->bayersensor.method;
        bayersensor.ccSteps = bayersensor.ccSteps && e0->bayersensor.ccSteps == e->bayersensor.ccSteps;
        bayersensor.exBlack0 = bayersensor.exBlack0 && e0->bayersensor.black0 == e->bayersensor.black0;
        bayersensor.exBlack1 = bayersensor.exBlack1 && e0->bayersensor.black1 == e->bayersensor.black1;
        bayersensor.exBlack2 = bayersensor.exBlack2 && e0->bayersensor.black2 == e->bayersensor.black2;
        bayersensor.exBlack3 = bayersensor.exBlack3 && e0->bayersensor.black3 == e->bayersensor.black3;
        bayersensor.exTwoGreen = bayersensor.exTwoGreen && e0->bayersensor.twogreen == e->bayersensor.twogreen;
        bayersensor.dcbIterations = bayersensor.dcbIterations && e0->bayersensor.dcb_iterations == e->bayersensor.dcb_iterations;
        bayersensor.dcbEnhance = bayersensor.dcbEnhance && e0->bayersensor.dcb_enhance == e->bayersensor.dcb_enhance;
        //bayersensor.allEnhance = bayersensor.allEnhance && e0->bayersensor.all_enhance == e->bayersensor.all_enhance;
        bayersensor.lmmseIterations = bayersensor.lmmseIterations && e0->bayersensor.lmmse_iterations == e->bayersensor.lmmse_iterations;
        bayersensor.greenEq = bayersensor.greenEq && e0->bayersensor.greenthresh == e->bayersensor.greenthresh;
        bayersensor.linenoise = bayersensor.linenoise && e0->bayersensor.linenoise == e->bayersensor.linenoise;
        xtranssensor.method = xtranssensor.method && e0->xtranssensor.method == e->xtranssensor.method;
        xtranssensor.ccSteps = xtranssensor.ccSteps && e0->xtranssensor.ccSteps == e->xtranssensor.ccSteps;
        xtranssensor.exBlackRed = xtranssensor.exBlackRed && e0->xtranssensor.blackred == e->xtranssensor.blackred;
        xtranssensor.exBlackGreen = xtranssensor.exBlackGreen && e0->xtranssensor.blackgreen == e->xtranssensor.blackgreen;
        xtranssensor.exBlackBlue = xtranssensor.exBlackBlue && e0->xtranssensor.blackblue == e->xtranssensor.blackblue;
        caCorrection = caCorrection && e0->ca_autocorrect == e->ca_autocorrect;
        caRed = caRed && e0->cared == e->cared;
        caBlue = caBlue && e0->cablue == e->cablue;
        hotPixelFilter = hotPixelFilter && e0->hotPixelFilter == e->hotPixelFilter;
        deadPixelFilter = deadPixelFilter && e0->deadPixelFilter == e->deadPixelFilter;
        hotDeadPixelThresh = hotDeadPixelThresh && e0->hotdeadpix_thresh == e->hotdeadpix_thresh;
        darkFrame = darkFrame && e0->dark_frame == e->dark_frame;
        dfAuto = dfAuto && e0->df_autoselect == e->df_autoselect;
        ff_file = ff_file && e0->ff_file == e->ff_file;
        ff_AutoSelect = ff_AutoSelect && e0->ff_AutoSelect == e->ff_AutoSelect;
        ff_BlurRadius = ff_BlurRadius && e0->ff_BlurRadius == e->ff_BlurRadius;
        ff_BlurType = ff_BlurType && e0->ff_BlurType == e->ff_BlurType;
        ff_AutoClipControl = ff_AutoClipControl && e0->ff_AutoClipControl == e->ff_AutoClipControl;
        ff_clipControl = ff_clipControl && e0->ff_clipControl == e->ff_clipControl;
        exPos = exPos && e0->expos == e->expos;
        exPreser = exPreser && e0->preser == e->preser;
    }
}

void DirPyrEqualizerParamsEdited::initFrom (std::vector<const void*> elems)
{
    const DirPyrEqualizerParams *e0 = static_cast< const DirPyrEqualizerParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const DirPyrEqualizerParams *e = static_cast< const DirPyrEqualizerParams* >(elems.at(i));
        enabled = enabled && e0->enabled == e->enabled;
        gamutlab = gamutlab && e0->gamutlab == e->gamutlab;
        for(int i = 0; i < 6; i++) {
            mult[i] = mult[i] && e0->mult[i] == e->mult[i];
        }
        threshold = threshold && e0->threshold == e->threshold;
        skinprotect = skinprotect && e0->skinprotect == e->skinprotect;
        //algo = algo && e0->algo == e->algo;
        hueskin = hueskin && e0->hueskin == e->hueskin;
    }
}

void WaveletParamsEdited::initFrom (std::vector<const void*> elems)
{
    const WaveletParams *e0 = static_cast< const WaveletParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const WaveletParams *e = static_cast< const WaveletParams* >(elems.at(i));
        enabled = enabled && e0->enabled == e->enabled;
        strength = strength && e0->strength == e->strength;
        balance = balance && e0->balance == e->balance;
        iter = iter && e0->iter == e->iter;
        median = median && e0->median == e->median;
        medianlev = medianlev && e0->medianlev == e->medianlev;
        linkedg = linkedg && e0->linkedg == e->linkedg;
        cbenab = cbenab && e0->cbenab == e->cbenab;
        greenmed = greenmed && e0->greenmed == e->greenmed;
        bluemed = bluemed && e0->bluemed == e->bluemed;
        greenhigh = greenhigh && e0->greenhigh == e->greenhigh;
        bluehigh = bluehigh && e0->bluehigh == e->bluehigh;
        greenlow = greenlow && e0->greenlow == e->greenlow;
        bluelow = bluelow && e0->bluelow == e->bluelow;
        lipst = lipst && e0->lipst == e->lipst;
        Medgreinf = Medgreinf && e0->Medgreinf == e->Medgreinf;
        avoid = avoid && e0->avoid == e->avoid;
        tmr = tmr && e0->tmr == e->tmr;
        Lmethod = Lmethod && e0->Lmethod == e->Lmethod;
        CLmethod = CLmethod && e0->CLmethod == e->CLmethod;
        Backmethod = Backmethod && e0->Backmethod == e->Backmethod;
        Tilesmethod = Tilesmethod && e0->Tilesmethod == e->Tilesmethod;
        daubcoeffmethod = daubcoeffmethod && e0->daubcoeffmethod == e->daubcoeffmethod;
        CHmethod = CHmethod && e0->CHmethod == e->CHmethod;
        CHSLmethod = CHSLmethod && e0->CHSLmethod == e->CHSLmethod;
        EDmethod = EDmethod && e0->EDmethod == e->EDmethod;
        NPmethod = NPmethod && e0->NPmethod == e->NPmethod;
        BAmethod = BAmethod && e0->BAmethod == e->BAmethod;
        TMmethod = TMmethod && e0->TMmethod == e->TMmethod;
        HSmethod = HSmethod && e0->HSmethod == e->HSmethod;
        Dirmethod = Dirmethod && e0->Dirmethod == e->Dirmethod;
        rescon = rescon && e0->rescon == e->rescon;
        resconH = resconH && e0->resconH == e->resconH;
        reschro = reschro && e0->reschro == e->reschro;
        tmrs = tmrs && e0->tmrs == e->tmrs;
        gamma = gamma && e0->gamma == e->gamma;
        sup = sup && e0->sup == e->sup;
        sky = sky && e0->sky == e->sky;
        threshold = threshold && e0->threshold == e->threshold;
        threshold2 = threshold2 && e0->threshold2 == e->threshold2;
        edgedetect = edgedetect && e0->edgedetect == e->edgedetect;
        edgedetectthr = edgedetectthr && e0->edgedetectthr == e->edgedetectthr;
        edgedetectthr2 = edgedetectthr2 && e0->edgedetectthr2 == e->edgedetectthr2;
        edgesensi = edgesensi && e0->edgesensi == e->edgesensi;
        edgeampli = edgeampli && e0->edgeampli == e->edgeampli;
        thres = thres && e0->thres == e->thres;
        chroma = chroma && e0->chroma == e->chroma;
        chro = chro && e0->chro == e->chro;
        contrast = contrast && e0->contrast == e->contrast;
        edgrad = edgrad && e0->edgrad == e->edgrad;
        edgval = edgval && e0->edgval == e->edgval;
        edgthresh = edgthresh && e0->edgthresh == e->edgthresh;
        thr = thr && e0->thr == e->thr;
        thrH = thrH && e0->thrH == e->thrH;
        hueskin = hueskin && e0->hueskin == e->hueskin;
        hueskin2 = hueskin2 && e0->hueskin2 == e->hueskin2;
        hllev = hllev && e0->hllev == e->hllev;
        bllev = bllev && e0->bllev == e->bllev;
        edgcont = edgcont && e0->edgcont == e->edgcont;
        level0noise = level0noise && e0->level0noise == e->level0noise;
        level1noise = level1noise && e0->level1noise == e->level1noise;
        level2noise = level2noise && e0->level2noise == e->level2noise;
        level3noise = level3noise && e0->level3noise == e->level3noise;
        pastlev = pastlev && e0->pastlev == e->pastlev;
        satlev = satlev && e0->satlev == e->satlev;
        ccwcurve = ccwcurve && e0->ccwcurve == e->ccwcurve;
        opacityCurveRG = opacityCurveRG && e0->opacityCurveRG == e->opacityCurveRG;
        opacityCurveBY = opacityCurveBY && e0->opacityCurveBY == e->opacityCurveBY;
        opacityCurveW = opacityCurveW && e0->opacityCurveW == e->opacityCurveW;
        opacityCurveWL = opacityCurveWL && e0->opacityCurveWL == e->opacityCurveWL;
        wavclCurve = wavclCurve && e0->wavclCurve == e->wavclCurve;
        hhcurve = hhcurve && e0->hhcurve == e->hhcurve;
        Chcurve = Chcurve && e0->Chcurve == e->Chcurve;
        skinprotect = skinprotect && e0->skinprotect == e->skinprotect;
        //enacont = enacont && e0->enacont == e->enacont;
        expcontrast = expcontrast && e0->expcontrast == e->expcontrast;
        expchroma = expchroma && e0->expchroma == e->expchroma;
        expedge = expedge && e0->expedge == e->expedge;
        expresid = expresid && e0->expresid == e->expresid;
        expfinal = expfinal && e0->expfinal == e->expfinal;
        exptoning = exptoning && e0->exptoning == e->exptoning;
        expnoise = expnoise && e0->expnoise == e->expnoise;
        for(int i = 0; i < 9; i++) {
            c[i] = c[i] && e0->c[i] == e->c[i];
        }
        for(int i = 0; i < 9; i++) {
            ch[i] = ch[i] && e0->ch[i] == e->ch[i];
        }
    }
}

void HSVEqualizerParamsEdited::initFrom (std::vector<const void*> elems)
{
    const HSVEqualizerParams *e0 = static_cast< const HSVEqualizerParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const HSVEqualizerParams *e = static_cast< const HSVEqualizerParams* >(elems.at(i));
        hcurve = hcurve && e0->hcurve == e->hcurve;
        scurve = scurve && e0->scurve == e->scurve;
        vcurve = vcurve && e0->vcurve == e->vcurve;
    }
}

void FilmSimulationParamsEdited::initFrom (std::vector<const void*> elems)
{
    const FilmSimulationParams *e0 = static_cast< const FilmSimulationParams* >(elems.at(0));
    size_t size = elems.size();
    for (size_t i = 1; i < size; ++i) {
        const FilmSimulationParams *e = static_cast< const FilmSimulationParams* >(elems.at(i));
        enabled = enabled && e0->enabled == e->enabled;
        clutFilename = clutFilename && e0->clutFilename == e->clutFilename;
        strength = strength && e0->strength == e->strength;
    }
}

void ParamsEdited::initFrom (const std::vector<rtengine::procparams::ProcParams>& src)
{

    set (true);

    if (src.empty()) {
        return;
    }

    size_t vectSize = src.size();
    std::vector<const void*> elems(vectSize);

    // *INDENT-OFF*

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).toneCurve);
    toneCurve.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).labCurve);
    labCurve.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).rgbCurves);
    rgbCurves.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).colorToning);
    colorToning.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).sharpenEdge);
    sharpenEdge.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).sharpenMicro);
    sharpenMicro.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).sharpening);
    sharpening.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).prsharpening);
    prsharpening.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).vibrance);
    vibrance.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).colorappearance);
    colorappearance.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).wb);
    wb.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).defringe);
    defringe.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).impulseDenoise);
    impulseDenoise.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).dirpyrDenoise);
    dirpyrDenoise.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).epd);
    epd.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).sh);
    sh.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).crop);
    crop.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).coarse);
    coarse.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).commonTrans);
    commonTrans.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).rotate);
    rotate.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).distortion);
    distortion.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).lensProf);
    lensProf.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).perspective);
    perspective.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).gradient);
    gradient.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).pcvignette);
    pcvignette.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).cacorrection);
    cacorrection.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).vignetting);
    vignetting.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).chmixer);
    chmixer.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).blackwhite);
    blackwhite.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).resize);
    resize.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).icm);
    icm.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).raw);
    raw.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).wavelet);
    wavelet.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).dirpyrequalizer);
    dirpyrequalizer.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).hsvequalizer);
    hsvequalizer.initFrom(elems);

    for (size_t i=0; i<vectSize; ++i) elems.at(i) = &(src.at(i).filmSimulation);
    filmSimulation.initFrom(elems);

    // *INDENT-ON*

    // Handling EXIF & IPTC
    /*
     * HOMBRE: How can we handle that ???
     *
    for (size_t i=0; i<vectSize; ++i) {
        const ProcParams& other = src[i];

        exif = exif && p.exif==other.exif
        iptc = other.iptc;
    }
    */
}

void ToneCurveParamsEdited::combine (ToneCurveParams* toEdit, const ToneCurveParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (curve) toEdit->curve = mods->curve;
    if (curve2) toEdit->curve2 = mods->curve2;
    if (curveMode) toEdit->curveMode = mods->curveMode;
    if (curveMode2) toEdit->curveMode2 = mods->curveMode2;
    if (brightness) toEdit->brightness = dontforceSet && options.baBehav[ADDSET_TC_BRIGHTNESS] ? toEdit->brightness + mods->brightness : mods->brightness;
    if (black) toEdit->black = dontforceSet && options.baBehav[ADDSET_TC_BLACKLEVEL] ? toEdit->black + mods->black : mods->black;
    if (contrast) toEdit->contrast = dontforceSet && options.baBehav[ADDSET_TC_CONTRAST] ? toEdit->contrast + mods->contrast : mods->contrast;
    if (saturation) toEdit->saturation = dontforceSet && options.baBehav[ADDSET_TC_SATURATION] ? toEdit->saturation + mods->saturation : mods->saturation;
    if (shcompr) toEdit->shcompr = dontforceSet && options.baBehav[ADDSET_TC_SHCOMP] ? toEdit->shcompr + mods->shcompr : mods->shcompr;
    if (autoexp) toEdit->autoexp = mods->autoexp;
    if (clip) toEdit->clip = mods->clip;
    if (expcomp) toEdit->expcomp = dontforceSet && options.baBehav[ADDSET_TC_EXPCOMP] ? toEdit->expcomp + mods->expcomp : mods->expcomp;
    if (hlcompr) toEdit->hlcompr = dontforceSet && options.baBehav[ADDSET_TC_HLCOMPAMOUNT] ? toEdit->hlcompr + mods->hlcompr : mods->hlcompr;
    if (hlcomprthresh) toEdit->hlcomprthresh = dontforceSet && options.baBehav[ADDSET_TC_HLCOMPTHRESH] ? toEdit->hlcomprthresh + mods->hlcomprthresh : mods->hlcomprthresh;
    if (hrenabled) toEdit->hrenabled = mods->hrenabled;
    if (method) toEdit->method = mods->method;
    // *INDENT-ON*
}

void LCurveParamsEdited::combine (LCurveParams* toEdit, const LCurveParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (lcurve) toEdit->lcurve = mods->lcurve;
    if (acurve) toEdit->acurve = mods->acurve;
    if (bcurve) toEdit->bcurve = mods->bcurve;
    if (cccurve) toEdit->cccurve = mods->cccurve;
    if (chcurve) toEdit->chcurve = mods->chcurve;
    if (lhcurve) toEdit->lhcurve = mods->lhcurve;
    if (hhcurve) toEdit->hhcurve = mods->hhcurve;
    if (lccurve) toEdit->lccurve = mods->lccurve;
    if (clcurve) toEdit->clcurve = mods->clcurve;
    if (brightness) toEdit->brightness = dontforceSet && options.baBehav[ADDSET_LC_BRIGHTNESS] ? toEdit->brightness + mods->brightness : mods->brightness;
    if (contrast) toEdit->contrast = dontforceSet && options.baBehav[ADDSET_LC_CONTRAST] ? toEdit->contrast + mods->contrast : mods->contrast;
    if (chromaticity) toEdit->chromaticity = dontforceSet && options.baBehav[ADDSET_LC_CHROMATICITY] ? toEdit->chromaticity + mods->chromaticity : mods->chromaticity;
    if (avoidcolorshift) toEdit->avoidcolorshift = mods->avoidcolorshift;
    if (rstprotection) toEdit->rstprotection = mods->rstprotection;
    if (lcredsk) toEdit->lcredsk = mods->lcredsk;
    // *INDENT-ON*
}

void RGBCurvesParamsEdited::combine (RGBCurvesParams* toEdit, const RGBCurvesParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (lumamode) toEdit->lumamode = mods->lumamode;
    if (rcurve) toEdit->rcurve = mods->rcurve;
    if (gcurve) toEdit->gcurve = mods->gcurve;
    if (bcurve) toEdit->bcurve = mods->bcurve;
    // *INDENT-ON*
}

void ColorToningEdited::combine (ColorToningParams* toEdit, const ColorToningParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (enabled) toEdit->enabled = mods->enabled;
    if (twocolor) toEdit->twocolor = mods->twocolor;
    if (opacityCurve) toEdit->opacityCurve = mods->opacityCurve;
    if (colorCurve) toEdit->colorCurve = mods->colorCurve;
    if (enabled) toEdit->enabled = mods->enabled;
    if (opacityCurve) toEdit->opacityCurve = mods->opacityCurve;
    if (satprotectionthreshold) toEdit->satProtectionThreshold = dontforceSet && options.baBehav[ADDSET_COLORTONING_SATTHRESHOLD] ? toEdit->satProtectionThreshold + mods->satProtectionThreshold : mods->satProtectionThreshold;
    if (autosat) toEdit->autosat = mods->autosat;
    if (saturatedopacity) toEdit->saturatedOpacity = dontforceSet && options.baBehav[ADDSET_COLORTONING_SATOPACITY] ? toEdit->saturatedOpacity + mods->saturatedOpacity : mods->saturatedOpacity;
    if (strength) toEdit->strength = dontforceSet && options.baBehav[ADDSET_COLORTONING_STRENGTH] ? toEdit->strength + mods->strength : mods->strength;
    if (shadowsColSat) toEdit->shadowsColSat = mods->shadowsColSat;
    if (hlColSat) toEdit->hlColSat = mods->hlColSat;
    if (balance) toEdit->balance = dontforceSet && options.baBehav[ADDSET_COLORTONING_BALANCE] ? toEdit->balance + mods->balance : mods->balance;
    if (clcurve) toEdit->clcurve = mods->clcurve;
    if (method) toEdit->method = mods->method;
    if (cl2curve) toEdit->cl2curve = mods->cl2curve;
    if (lumamode) toEdit->lumamode = mods->lumamode;
    if (satlow) toEdit->satlow = dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit->satlow + mods->satlow : mods->satlow;
    if (sathigh) toEdit->sathigh = dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit->sathigh + mods->sathigh : mods->sathigh;
    if (redlow) toEdit->redlow = dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit->redlow + mods->redlow : mods->redlow;
    if (greenlow) toEdit->greenlow = dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit->greenlow + mods->greenlow : mods->greenlow;
    if (bluelow) toEdit->bluelow = dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit->bluelow + mods->bluelow : mods->bluelow;
    if (redmed) toEdit->redmed = dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit->redmed + mods->redmed : mods->redmed;
    if (greenmed) toEdit->greenmed = dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit->greenmed + mods->greenmed : mods->greenmed;
    if (bluemed) toEdit->bluemed = dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit->bluemed + mods->bluemed : mods->bluemed;
    if (redhigh) toEdit->redhigh = dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit->redhigh + mods->redhigh : mods->redhigh;
    if (greenhigh) toEdit->greenhigh = dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit->greenhigh + mods->greenhigh : mods->greenhigh;
    if (bluehigh) toEdit->bluehigh = dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit->bluehigh + mods->bluehigh : mods->bluehigh;
    // *INDENT-ON*
}

void SharpeningParamsEdited::combine (SharpeningParams* toEdit, const SharpeningParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (enabled) toEdit->enabled = mods->enabled;
    if (radius) toEdit->radius = mods->radius;
    if (amount) toEdit->amount = dontforceSet && options.baBehav[ADDSET_SHARP_AMOUNT] ? toEdit->amount + mods->amount : mods->amount;
    if (threshold) toEdit->threshold = mods->threshold;
    if (edgesonly) toEdit->edgesonly = mods->edgesonly;
    if (edges_radius) toEdit->edges_radius = mods->edges_radius;
    if (edges_tolerance) toEdit->edges_tolerance = mods->edges_tolerance;
    if (halocontrol) toEdit->halocontrol = mods->halocontrol;
    if (halocontrol_amount) toEdit->halocontrol_amount = mods->halocontrol_amount;
    if (method) toEdit->method = mods->method;
    if (deconvamount) toEdit->deconvamount = dontforceSet && options.baBehav[ADDSET_SHARP_AMOUNT] ? toEdit->deconvamount + mods->deconvamount : mods->deconvamount;
    if (deconvradius) toEdit->deconvradius = mods->deconvradius;
    if (deconviter) toEdit->deconviter = mods->deconviter;
    if (deconvdamping) toEdit->deconvdamping = mods->deconvdamping;
    // *INDENT-ON*
}

void SharpenEdgeParamsEdited::combine (SharpenEdgeParams* toEdit, const SharpenEdgeParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (enabled) toEdit->enabled = mods->enabled;
    if (passes) toEdit->passes = dontforceSet && options.baBehav[ADDSET_SHARPENEDGE_PASS] ? toEdit->passes + mods->passes : mods->passes;
    if (amount) toEdit->amount = dontforceSet && options.baBehav[ADDSET_SHARPENEDGE_AMOUNT] ? toEdit->amount + mods->amount : mods->amount;
    if (threechannels) toEdit->threechannels = mods->threechannels;
    // *INDENT-ON*
}

void SharpenMicroParamsEdited::combine (SharpenMicroParams* toEdit, const SharpenMicroParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (enabled) toEdit->enabled = mods->enabled;
    if (matrix) toEdit->matrix = mods->matrix;
    if (amount) toEdit->amount = dontforceSet && options.baBehav[ADDSET_SHARPENMICRO_AMOUNT] ? toEdit->amount + mods->amount : mods->amount;
    if (uniformity) toEdit->uniformity = dontforceSet && options.baBehav[ADDSET_SHARPENMICRO_UNIFORMITY] ? toEdit->uniformity + mods->uniformity : mods->uniformity;
    // *INDENT-ON*
}

void VibranceParamsEdited::combine (VibranceParams* toEdit, const VibranceParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (enabled) toEdit->enabled = mods->enabled;
    if (pastels) toEdit->pastels = dontforceSet && options.baBehav[ADDSET_VIBRANCE_PASTELS] ? toEdit->pastels + mods->pastels : mods->pastels;
    if (saturated) toEdit->saturated = dontforceSet && options.baBehav[ADDSET_VIBRANCE_SATURATED] ? toEdit->saturated + mods->saturated : mods->saturated;
    if (psthreshold) toEdit->psthreshold = mods->psthreshold;
    if (protectskins) toEdit->protectskins = mods->protectskins;
    if (avoidcolorshift) toEdit->avoidcolorshift = mods->avoidcolorshift;
    if (pastsattog) toEdit->pastsattog = mods->pastsattog;
    if (skintonescurve) toEdit->skintonescurve = mods->skintonescurve;
    // *INDENT-ON*
}

void ColorAppearanceParamsEdited::combine (ColorAppearanceParams* toEdit, const ColorAppearanceParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (curve) toEdit->curve = mods->curve;
    if (curve2) toEdit->curve2 = mods->curve2;
    if (curve3) toEdit->curve3 = mods->curve3;
    if (curveMode) toEdit->curveMode = mods->curveMode;
    if (curveMode2) toEdit->curveMode2 = mods->curveMode2;
    if (curveMode3) toEdit->curveMode3 = mods->curveMode3;
    if (enabled) toEdit->enabled = mods->enabled;
    if (degree) toEdit->degree = dontforceSet && options.baBehav[ADDSET_CAT_DEGREE] ? toEdit->degree + mods->degree : mods->degree;
    if (autodegree) toEdit->autodegree = mods->autodegree;
    if (surround) toEdit->surround = mods->surround;
    if (autoadapscen) toEdit->autoadapscen = mods->autoadapscen;
    if (adapscen) toEdit->adapscen = mods->adapscen;
    if (adaplum) toEdit->adaplum = dontforceSet && options.baBehav[ADDSET_CAT_ADAPTVIEWING] ? toEdit->adaplum + mods->adaplum : mods->adaplum;
    if (badpixsl) toEdit->badpixsl = dontforceSet && options.baBehav[ADDSET_CAT_BADPIX] ? toEdit->badpixsl + mods->badpixsl : mods->badpixsl;
    if (wbmodel) toEdit->wbmodel = mods->wbmodel;
    if (algo) toEdit->algo = mods->algo;
    if (jlight) toEdit->jlight = dontforceSet && options.baBehav[ADDSET_CAT_LIGHT] ? toEdit->jlight + mods->jlight : mods->jlight;
    if (qbright) toEdit->qbright = dontforceSet && options.baBehav[ADDSET_CAT_BRIGHT] ? toEdit->qbright + mods->qbright : mods->qbright;
    if (chroma) toEdit->chroma = dontforceSet && options.baBehav[ADDSET_CAT_CHROMA] ? toEdit->chroma + mods->chroma : mods->chroma;
    if (schroma) toEdit->schroma = dontforceSet && options.baBehav[ADDSET_CAT_CHROMA_S] ? toEdit->schroma + mods->schroma : mods->schroma;
    if (mchroma) toEdit->mchroma = dontforceSet && options.baBehav[ADDSET_CAT_CHROMA_M] ? toEdit->mchroma + mods->mchroma : mods->mchroma;
    if (contrast) toEdit->contrast = dontforceSet && options.baBehav[ADDSET_CAT_CONTRAST] ? toEdit->contrast + mods->contrast : mods->contrast;
    if (qcontrast) toEdit->qcontrast = dontforceSet && options.baBehav[ADDSET_CAT_CONTRAST_Q] ? toEdit->qcontrast + mods->qcontrast : mods->qcontrast;
    if (colorh) toEdit->colorh = dontforceSet && options.baBehav[ADDSET_CAT_HUE] ? toEdit->colorh + mods->colorh : mods->colorh;
    if (rstprotection) toEdit->rstprotection = dontforceSet && options.baBehav[ADDSET_CAT_RSTPRO] ? toEdit->rstprotection + mods->rstprotection : mods->rstprotection;
    if (surrsource) toEdit->surrsource = mods->surrsource;
    if (gamut) toEdit->gamut = mods->gamut;
    // if (badpix) toEdit->badpix = mods->badpix;
    if (datacie) toEdit->datacie = mods->datacie;
    if (tonecie) toEdit->tonecie = mods->tonecie;
    // if (sharpcie) toEdit->sharpcie = mods->sharpcie;
    // *INDENT-ON*
}

void WBParamsEdited::combine (WBParams* toEdit, const WBParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (method) toEdit->method = mods->method;
    if (equal) toEdit->equal = dontforceSet && options.baBehav[ADDSET_WB_EQUAL] ? toEdit->equal + mods->equal : mods->equal;
    if (green) toEdit->green = dontforceSet && options.baBehav[ADDSET_WB_GREEN] ? toEdit->green + mods->green : mods->green;
    if (temperature) toEdit->temperature = dontforceSet && options.baBehav[ADDSET_WB_TEMPERATURE] ? toEdit->temperature + mods->temperature : mods->temperature;
    // *INDENT-ON*
}

/*
void ColorShiftParamsEdited::combine (ColorShiftParams* toEdit, const ColorShiftParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (a) toEdit->a = dontforceSet && options.baBehav[ADDSET_CS_BLUEYELLOW] ? toEdit->a + mods->a : mods->a;
    if (b) toEdit->b = dontforceSet && options.baBehav[ADDSET_CS_GREENMAGENTA] ? toEdit->b + mods->b : mods->b;
    // *INDENT-ON*
}

void LumaDenoiseParamsEdited::combine (LumaDenoiseParams* toEdit, const LumaDenoiseParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (enabled) toEdit->enabled = mods->enabled;
    if (radius) toEdit->radius = mods->radius;
    if (edgetolerance) toEdit->edgetolerance = dontforceSet && options.baBehav[ADDSET_LD_EDGETOLERANCE] ? toEdit->edgetolerance + mods->edgetolerance : mods->edgetolerance;
    // *INDENT-ON*
}

void ColorDenoiseParamsEdited::combine (ColorDenoiseParams* toEdit, const ColorDenoiseParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (enabled) toEdit->enabled = mods->enabled;
    if (amount) toEdit->amount = mods->amount;
    // *INDENT-ON*
}
*/

void DefringeParamsEdited::combine (DefringeParams* toEdit, const DefringeParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (enabled) toEdit->enabled = mods->enabled;
    if (radius) toEdit->radius = mods->radius;
    if (threshold) toEdit->threshold = mods->threshold;
    if (huecurve) toEdit->huecurve = mods->huecurve;
    // *INDENT-ON*
}

void DirPyrDenoiseParamsEdited::combine (DirPyrDenoiseParams* toEdit, const DirPyrDenoiseParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (enabled) toEdit->enabled = mods->enabled;
    if (enhance) toEdit->enhance = mods->enhance;
    if (median) toEdit->median = mods->median;
    if (autochroma) toEdit->autochroma = mods->autochroma;
    if (luma) toEdit->luma = dontforceSet && options.baBehav[ADDSET_DIRPYRDN_LUMA] ? toEdit->luma + mods->luma : mods->luma;
    if (lcurve) toEdit->lcurve = mods->lcurve;
    if (cccurve) toEdit->cccurve = mods->cccurve;
    if (Ldetail) toEdit->Ldetail = dontforceSet && options.baBehav[ADDSET_DIRPYRDN_LUMDET] ? toEdit->Ldetail + mods->Ldetail : mods->Ldetail;
    if (chroma) toEdit->chroma = dontforceSet && options.baBehav[ADDSET_DIRPYRDN_CHROMA] ? toEdit->chroma + mods->chroma : mods->chroma;
    if (redchro) toEdit->redchro = dontforceSet && options.baBehav[ADDSET_DIRPYRDN_CHROMARED] ? toEdit->redchro + mods->redchro : mods->redchro;
    if (bluechro) toEdit->bluechro = dontforceSet && options.baBehav[ADDSET_DIRPYRDN_CHROMABLUE] ? toEdit->bluechro + mods->bluechro : mods->bluechro;
    if (gamma) toEdit->gamma = dontforceSet && options.baBehav[ADDSET_DIRPYRDN_GAMMA] ? toEdit->gamma + mods->gamma : mods->gamma;
    if (passes) toEdit->passes = dontforceSet && options.baBehav[ADDSET_DIRPYRDN_PASSES] ? toEdit->passes + mods->passes : mods->passes;
    // if (perform) toEdit->perform = mods->perform;
    if (dmethod) toEdit->dmethod = mods->dmethod;
    if (Lmethod) toEdit->Lmethod = mods->Lmethod;
    if (Cmethod) toEdit->Cmethod = mods->Cmethod;
    if (C2method) toEdit->C2method = mods->C2method;
    if (smethod) toEdit->smethod = mods->smethod;
    if (medmethod) toEdit->medmethod = mods->medmethod;
    if (methodmed) toEdit->methodmed = mods->methodmed;
    if (rgbmethod) toEdit->rgbmethod = mods->rgbmethod;
    // *INDENT-ON*
}

void EPDParamsEdited::combine (EPDParams* toEdit, const EPDParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (enabled) toEdit->enabled = mods->enabled;
    if (strength) toEdit->strength = mods->strength;
    if (gamma) toEdit->gamma = mods->gamma;
    if (edgeStopping) toEdit->edgeStopping = mods->edgeStopping;
    if (scale) toEdit->scale = mods->scale;
    if (reweightingIterates) toEdit->reweightingIterates = mods->reweightingIterates;
    // *INDENT-ON*
}

void ImpulseDenoiseParamsEdited::combine (ImpulseDenoiseParams* toEdit, const ImpulseDenoiseParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (enabled) toEdit->enabled = mods->enabled;
    if (thresh) toEdit->thresh = mods->thresh;
    // *INDENT-ON*
}

void SHParamsEdited::combine (SHParams* toEdit, const SHParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (enabled) toEdit->enabled = mods->enabled;
    if (hq) toEdit->hq = mods->hq;
    if (highlights) toEdit->highlights = dontforceSet && options.baBehav[ADDSET_SH_HIGHLIGHTS] ? toEdit->highlights + mods->highlights : mods->highlights;
    if (htonalwidth) toEdit->htonalwidth = mods->htonalwidth;
    if (shadows) toEdit->shadows = dontforceSet && options.baBehav[ADDSET_SH_SHADOWS] ? toEdit->shadows + mods->shadows : mods->shadows;
    if (stonalwidth) toEdit->stonalwidth = mods->stonalwidth;
    if (localcontrast) toEdit->localcontrast = dontforceSet && options.baBehav[ADDSET_SH_LOCALCONTRAST] ? toEdit->localcontrast + mods->localcontrast : mods->localcontrast;
    if (radius) toEdit->radius = mods->radius;
    // *INDENT-ON*
}

void CropParamsEdited::combine (CropParams* toEdit, const CropParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (enabled) toEdit->enabled = mods->enabled;
    if (x) toEdit->x = mods->x;
    if (y) toEdit->y = mods->y;
    if (w) toEdit->w = mods->w;
    if (h) toEdit->h = mods->h;
    if (fixratio) toEdit->fixratio = mods->fixratio;
    if (ratio) toEdit->ratio = mods->ratio;
    if (orientation) toEdit->orientation = mods->orientation;
    if (guide) toEdit->guide = mods->guide;
    // *INDENT-ON*
}

void CoarseTransformParamsEdited::combine (CoarseTransformParams* toEdit, const CoarseTransformParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (rotate) toEdit->rotate = mods->rotate;
    if (hflip) toEdit->hflip = mods->hflip;
    if (vflip) toEdit->vflip = mods->vflip;
    // *INDENT-ON*
}

void CommonTransformParamsEdited::combine (CommonTransformParams* toEdit, const CommonTransformParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (autofill) toEdit->autofill = mods->autofill;
    // *INDENT-ON*
}

void RotateParamsEdited::combine (RotateParams* toEdit, const RotateParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (degree) toEdit->degree = dontforceSet && options.baBehav[ADDSET_ROTATE_DEGREE] ? toEdit->degree + mods->degree : mods->degree;
    // *INDENT-ON*
}

void DistortionParamsEdited::combine (DistortionParams* toEdit, const DistortionParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (amount) toEdit->amount = dontforceSet && options.baBehav[ADDSET_DIST_AMOUNT] ? toEdit->amount + mods->amount : mods->amount;
    // *INDENT-ON*
}

void LensProfParamsEdited::combine (LensProfParams* toEdit, const LensProfParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (lcpFile) toEdit->lcpFile = mods->lcpFile;
    if (useDist) toEdit->useDist = mods->useDist;
    if (useVign) toEdit->useVign = mods->useVign;
    if (useCA) toEdit->useCA = mods->useCA;
    // *INDENT-ON*
}

void PerspectiveParamsEdited::combine (PerspectiveParams* toEdit, const PerspectiveParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (horizontal) toEdit->horizontal = dontforceSet && options.baBehav[ADDSET_PERSPECTIVE] ? toEdit->horizontal + mods->horizontal : mods->horizontal;
    if (vertical) toEdit->vertical = dontforceSet && options.baBehav[ADDSET_PERSPECTIVE] ? toEdit->vertical + mods->vertical : mods->vertical;
    // *INDENT-ON*
}

void GradientParamsEdited::combine (GradientParams* toEdit, const GradientParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (enabled) toEdit->enabled = mods->enabled;
    if (degree) toEdit->degree = dontforceSet && options.baBehav[ADDSET_GRADIENT_DEGREE] ? toEdit->degree + mods->degree : mods->degree;
    if (feather) toEdit->feather = mods->feather;
    if (strength) toEdit->strength = mods->strength;
    if (centerX) toEdit->centerX = mods->centerX;
    if (centerY) toEdit->centerY = mods->centerY;
    // *INDENT-ON*
}

void PCVignetteParamsEdited::combine (PCVignetteParams* toEdit, const PCVignetteParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (enabled) toEdit->enabled = mods->enabled;
    if (strength) toEdit->strength = mods->strength;
    if (feather) toEdit->feather = mods->feather;
    if (roundness) toEdit->roundness = mods->roundness;
    // *INDENT-ON*
}

void CACorrParamsEdited::combine (CACorrParams* toEdit, const CACorrParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (red) toEdit->red = dontforceSet && options.baBehav[ADDSET_CA] ? toEdit->red + mods->red : mods->red;
    if (blue) toEdit->blue = dontforceSet && options.baBehav[ADDSET_CA] ? toEdit->blue + mods->blue : mods->blue;
    // *INDENT-ON*
}

void VignettingParamsEdited::combine (VignettingParams* toEdit, const VignettingParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (amount) toEdit->amount = dontforceSet && options.baBehav[ADDSET_VIGN_AMOUNT] ? toEdit->amount + mods->amount : mods->amount;
    if (radius) toEdit->radius = mods->radius;
    if (strength) toEdit->strength = mods->strength;
    if (centerX) toEdit->centerX = mods->centerX;
    if (centerY) toEdit->centerY = mods->centerY;
    // *INDENT-ON*
}

void ChannelMixerParamsEdited::combine (ChannelMixerParams* toEdit, const ChannelMixerParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    for (int i = 0; i < 3; i++) {
        if (red[i]) toEdit->red[i] = dontforceSet && options.baBehav[ADDSET_CHMIXER] ? toEdit->red[i] + mods->red[i] : mods->red[i];
        if (green[i]) toEdit->green[i] = dontforceSet && options.baBehav[ADDSET_CHMIXER] ? toEdit->green[i] + mods->green[i] : mods->green[i];
        if (blue[i]) toEdit->blue[i] = dontforceSet && options.baBehav[ADDSET_CHMIXER] ? toEdit->blue[i] + mods->blue[i] : mods->blue[i];
    }
    // *INDENT-ON*
}

void BlackWhiteParamsEdited::combine (BlackWhiteParams* toEdit, const BlackWhiteParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (enabled) toEdit->enabled = mods->enabled;
    if (method) toEdit->method = mods->method;
    if (luminanceCurve) toEdit->luminanceCurve = mods->luminanceCurve;
    if (autoc) toEdit->autoc = mods->autoc;
    if (setting) toEdit->setting = mods->setting;
    if (enabledcc) toEdit->enabledcc = mods->enabledcc;
    if (filter) toEdit->filter = mods->filter;
    if (mixerRed) toEdit->mixerRed = dontforceSet && options.baBehav[ADDSET_BLACKWHITE_HUES] ? toEdit->mixerRed + mods->mixerRed : mods->mixerRed;
    if (mixerOrange) toEdit->mixerOrange = dontforceSet && options.baBehav[ADDSET_BLACKWHITE_HUES] ? toEdit->mixerOrange + mods->mixerOrange : mods->mixerOrange;
    if (mixerYellow) toEdit->mixerYellow = dontforceSet && options.baBehav[ADDSET_BLACKWHITE_HUES] ? toEdit->mixerYellow + mods->mixerYellow : mods->mixerYellow;
    if (mixerGreen) toEdit->mixerGreen = dontforceSet && options.baBehav[ADDSET_BLACKWHITE_HUES] ? toEdit->mixerGreen + mods->mixerGreen : mods->mixerGreen;
    if (mixerCyan) toEdit->mixerCyan = dontforceSet && options.baBehav[ADDSET_BLACKWHITE_HUES] ? toEdit->mixerCyan + mods->mixerCyan : mods->mixerCyan;
    if (mixerBlue) toEdit->mixerBlue = dontforceSet && options.baBehav[ADDSET_BLACKWHITE_HUES] ? toEdit->mixerBlue + mods->mixerBlue : mods->mixerBlue;
    if (mixerMagenta) toEdit->mixerMagenta = dontforceSet && options.baBehav[ADDSET_BLACKWHITE_HUES] ? toEdit->mixerMagenta + mods->mixerMagenta : mods->mixerMagenta;
    if (mixerPurple) toEdit->mixerPurple = dontforceSet && options.baBehav[ADDSET_BLACKWHITE_HUES] ? toEdit->mixerPurple + mods->mixerPurple : mods->mixerPurple;
    if (gammaRed) toEdit->gammaRed = dontforceSet && options.baBehav[ADDSET_BLACKWHITE_GAMMA] ? toEdit->gammaRed + mods->gammaRed : mods->gammaRed;
    if (gammaGreen) toEdit->gammaGreen = dontforceSet && options.baBehav[ADDSET_BLACKWHITE_GAMMA] ? toEdit->gammaGreen + mods->gammaGreen : mods->gammaGreen;
    if (gammaBlue) toEdit->gammaBlue = dontforceSet && options.baBehav[ADDSET_BLACKWHITE_GAMMA] ? toEdit->gammaBlue + mods->gammaBlue : mods->gammaBlue;
    if (beforeCurve) toEdit->beforeCurve = mods->beforeCurve;
    if (beforeCurveMode) toEdit->beforeCurveMode = mods->beforeCurveMode;
    if (afterCurve) toEdit->afterCurve = mods->afterCurve;
    if (afterCurveMode) toEdit->afterCurveMode = mods->afterCurveMode;
    if (algo) toEdit->algo = mods->algo;
    // *INDENT-ON*
}

void ResizeParamsEdited::combine (ResizeParams* toEdit, const ResizeParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (scale) toEdit->scale = mods->scale;
    if (appliesTo) toEdit->appliesTo = mods->appliesTo;
    if (method) toEdit->method = mods->method;
    if (dataspec) toEdit->dataspec = mods->dataspec;
    if (width) toEdit->width = mods->width;
    if (height) toEdit->height = mods->height;
    if (enabled) toEdit->enabled = mods->enabled;
    // *INDENT-ON*
}

void ColorManagementParamsEdited::combine (ColorManagementParams* toEdit, const ColorManagementParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (input) toEdit->input = mods->input;
    if (toneCurve) toEdit->toneCurve = mods->toneCurve;
    if (applyLookTable) toEdit->applyLookTable = mods->applyLookTable;
    if (applyBaselineExposureOffset) toEdit->applyBaselineExposureOffset = mods->applyBaselineExposureOffset;
    if (applyHueSatMap) toEdit->applyHueSatMap = mods->applyHueSatMap;
    if (blendCMSMatrix) toEdit->blendCMSMatrix = mods->blendCMSMatrix;
    if (dcpIlluminant) toEdit->dcpIlluminant = mods->dcpIlluminant;
    if (working) toEdit->working = mods->working;
    if (output) toEdit->output = mods->output;
    //if (gampos) toEdit->gampos = mods->gampos;
    //if (slpos) toEdit->slpos = mods->slpos;
    if (gampos) toEdit->gampos = dontforceSet && options.baBehav[ADDSET_FREE_OUPUT_GAMMA] ? toEdit->gampos + mods->gampos : mods->gampos;
    if (slpos) toEdit->slpos = dontforceSet && options.baBehav[ADDSET_FREE_OUTPUT_SLOPE] ? toEdit->slpos + mods->slpos : mods->slpos;
    if (gamma) toEdit->gamma = mods->gamma;
    if (freegamma) toEdit->freegamma = mods->freegamma;
    // *INDENT-ON*
}

void RAWParamsEdited::combine (RAWParams* toEdit, const RAWParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (bayersensor.method) toEdit->bayersensor.method = mods->bayersensor.method;
    if (bayersensor.ccSteps) toEdit->bayersensor.ccSteps = mods->bayersensor.ccSteps;
    if (bayersensor.exBlack0) toEdit->bayersensor.black0 = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_BLACKS] ? toEdit->bayersensor.black0 + mods->bayersensor.black0 : mods->bayersensor.black0;
    if (bayersensor.exBlack1) toEdit->bayersensor.black1 = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_BLACKS] ? toEdit->bayersensor.black1 + mods->bayersensor.black1 : mods->bayersensor.black1;
    if (bayersensor.exBlack2) toEdit->bayersensor.black2 = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_BLACKS] ? toEdit->bayersensor.black2 + mods->bayersensor.black2 : mods->bayersensor.black2;
    if (bayersensor.exBlack3) toEdit->bayersensor.black3 = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_BLACKS] ? toEdit->bayersensor.black3 + mods->bayersensor.black3 : mods->bayersensor.black3;
    if (bayersensor.exTwoGreen) toEdit->bayersensor.twogreen = mods->bayersensor.twogreen;
    if (bayersensor.dcbIterations) toEdit->bayersensor.dcb_iterations = mods->bayersensor.dcb_iterations;
    if (bayersensor.dcbEnhance) toEdit->bayersensor.dcb_enhance = mods->bayersensor.dcb_enhance;
    if (bayersensor.lmmseIterations) toEdit->bayersensor.lmmse_iterations = mods->bayersensor.lmmse_iterations;
    //if (bayersensor.allEnhance) toEdit->bayersensor.all_enhance = mods->bayersensor.all_enhance;
    if (bayersensor.greenEq) toEdit->bayersensor.greenthresh = dontforceSet && options.baBehav[ADDSET_PREPROCESS_GREENEQUIL] ? toEdit->bayersensor.greenthresh + mods->bayersensor.greenthresh : mods->bayersensor.greenthresh;
    if (bayersensor.linenoise) toEdit->bayersensor.linenoise = dontforceSet && options.baBehav[ADDSET_PREPROCESS_LINEDENOISE] ? toEdit->bayersensor.linenoise + mods->bayersensor.linenoise : mods->bayersensor.linenoise;
    if (xtranssensor.method) toEdit->xtranssensor.method = mods->xtranssensor.method;
    if (xtranssensor.ccSteps) toEdit->xtranssensor.ccSteps = mods->xtranssensor.ccSteps;
    if (xtranssensor.exBlackRed) toEdit->xtranssensor.blackred = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_BLACKS] ? toEdit->xtranssensor.blackred + mods->xtranssensor.blackred : mods->xtranssensor.blackred;
    if (xtranssensor.exBlackGreen) toEdit->xtranssensor.blackgreen = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_BLACKS] ? toEdit->xtranssensor.blackgreen + mods->xtranssensor.blackgreen : mods->xtranssensor.blackgreen;
    if (xtranssensor.exBlackBlue) toEdit->xtranssensor.blackblue = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_BLACKS] ? toEdit->xtranssensor.blackblue + mods->xtranssensor.blackblue : mods->xtranssensor.blackblue;
    if (caCorrection) toEdit->ca_autocorrect = mods->ca_autocorrect;
    if (caRed) toEdit->cared = dontforceSet && options.baBehav[ADDSET_RAWCACORR] ? toEdit->cared + mods->cared : mods->cared;
    if (caBlue) toEdit->cablue = dontforceSet && options.baBehav[ADDSET_RAWCACORR] ? toEdit->cablue + mods->cablue : mods->cablue;
    if (exPos) toEdit->expos = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_LINEAR] ? toEdit->expos + mods->expos : mods->expos;
    if (exPreser) toEdit->preser = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_PRESER] ? toEdit->preser + mods->preser : mods->preser;
    if (hotPixelFilter) toEdit->hotPixelFilter = mods->hotPixelFilter;
    if (deadPixelFilter) toEdit->deadPixelFilter = mods->deadPixelFilter;
    if (hotDeadPixelThresh) toEdit->hotdeadpix_thresh = mods->hotdeadpix_thresh;
    if (darkFrame) toEdit->dark_frame = mods->dark_frame;
    if (dfAuto) toEdit->df_autoselect = mods->df_autoselect;
    if (ff_file) toEdit->ff_file = mods->ff_file;
    if (ff_AutoSelect) toEdit->ff_AutoSelect = mods->ff_AutoSelect;
    if (ff_BlurRadius) toEdit->ff_BlurRadius = mods->ff_BlurRadius;
    if (ff_BlurType) toEdit->ff_BlurType = mods->ff_BlurType;
    if (ff_AutoClipControl) toEdit->ff_AutoClipControl = mods->ff_AutoClipControl;
    if (ff_clipControl) toEdit->ff_clipControl = dontforceSet && options.baBehav[ADDSET_RAWFFCLIPCONTROL] ? toEdit->ff_clipControl + mods->ff_clipControl : mods->ff_clipControl;
    // *INDENT-ON*
}

void DirPyrEqualizerParamsEdited::combine (DirPyrEqualizerParams* toEdit, const DirPyrEqualizerParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (enabled) toEdit->enabled = mods->enabled;
    if (gamutlab) toEdit->gamutlab = mods->gamutlab;
    for(int i = 0; i < 6; i++) {
        if(mult[i]) toEdit->mult[i] = dontforceSet && options.baBehav[ADDSET_DIRPYREQ] ? toEdit->mult[i] + mods->mult[i] : mods->mult[i];
    }
    if (threshold) toEdit->threshold = dontforceSet && options.baBehav[ADDSET_DIRPYREQ_THRESHOLD] ? toEdit->threshold + mods->threshold : mods->threshold;
    if (skinprotect) toEdit->skinprotect = dontforceSet && options.baBehav[ADDSET_DIRPYREQ_SKINPROTECT] ? toEdit->skinprotect + mods->skinprotect : mods->skinprotect;
    if (hueskin) toEdit->hueskin = mods->hueskin;
    // if (algo) toEdit->algo = mods->algo;
    // *INDENT-ON*
}

void WaveletParamsEdited::combine (WaveletParams* toEdit, const WaveletParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (enabled) toEdit->enabled = mods->enabled;
    if (strength) toEdit->strength = mods->strength;
    if (balance) toEdit->balance = mods->balance;
    if (iter) toEdit->iter = mods->iter;
    if (median) toEdit->median = mods->median;
    if (medianlev) toEdit->medianlev = mods->medianlev;
    if (linkedg) toEdit->linkedg = mods->linkedg;
    if (cbenab) toEdit->cbenab = mods->cbenab;
    if (greenhigh) toEdit->greenhigh = mods->greenhigh;
    if (bluehigh) toEdit->bluehigh = mods->bluehigh;
    if (greenmed) toEdit->greenmed = mods->greenmed;
    if (bluemed) toEdit->bluemed = mods->bluemed;
    if (greenlow) toEdit->greenlow = mods->greenlow;
    if (bluelow) toEdit->bluelow = mods->bluelow;
    if (lipst) toEdit->lipst = mods->lipst;
    if (Medgreinf) toEdit->Medgreinf = mods->Medgreinf;
    if (avoid) toEdit->avoid = mods->avoid;
    if (tmr) toEdit->tmr = mods->tmr;
    if (Lmethod) toEdit->Lmethod = mods->Lmethod;
    if (CLmethod) toEdit->CLmethod = mods->CLmethod;
    if (Backmethod) toEdit->Backmethod = mods->Backmethod;
    if (Tilesmethod) toEdit->Tilesmethod = mods->Tilesmethod;
    if (daubcoeffmethod) toEdit->daubcoeffmethod = mods->daubcoeffmethod;
    if (CHmethod) toEdit->CHmethod = mods->CHmethod;
    if (CHSLmethod) toEdit->CHSLmethod = mods->CHSLmethod;
    if (EDmethod) toEdit->EDmethod = mods->EDmethod;
    if (NPmethod) toEdit->NPmethod = mods->NPmethod;
    if (BAmethod) toEdit->BAmethod = mods->BAmethod;
    if (TMmethod) toEdit->TMmethod = mods->TMmethod;
    if (HSmethod) toEdit->HSmethod = mods->HSmethod;
    if (Dirmethod) toEdit->Dirmethod = mods->Dirmethod;
    if (edgthresh) toEdit->edgthresh = mods->edgthresh;
    if (sky) toEdit->sky = dontforceSet && options.baBehav[ADDSET_WA_SKYPROTECT] ? toEdit->sky + mods->sky : mods->sky;
    if (thr) toEdit->thr = dontforceSet && options.baBehav[ADDSET_WA_THRR] ? toEdit->thr + mods->thr : mods->thr;
    if (thrH) toEdit->thrH = dontforceSet && options.baBehav[ADDSET_WA_THRRH] ? toEdit->thrH + mods->thrH : mods->thrH;
    if (sup) toEdit->sup = mods->sup;
    if (hllev) toEdit->hllev = mods->hllev;
    if (bllev) toEdit->bllev = mods->bllev;
    if (edgcont) toEdit->edgcont = mods->edgcont;
    if (level0noise) toEdit->level0noise = mods->level0noise;
    if (level1noise) toEdit->level1noise = mods->level1noise;
    if (level2noise) toEdit->level2noise = mods->level2noise;
    if (level3noise) toEdit->level3noise = mods->level3noise;
    if (pastlev) toEdit->pastlev = mods->pastlev;
    if (satlev) toEdit->satlev = mods->satlev;
    if (ccwcurve) toEdit->ccwcurve = mods->ccwcurve;
    if (opacityCurveRG) toEdit->opacityCurveRG = mods->opacityCurveRG;
    if (opacityCurveBY) toEdit->opacityCurveBY = mods->opacityCurveBY;
    if (opacityCurveW) toEdit->opacityCurveW = mods->opacityCurveW;
    if (opacityCurveWL) toEdit->opacityCurveWL = mods->opacityCurveWL;
    if (hhcurve) toEdit->hhcurve = mods->hhcurve;
    if (Chcurve) toEdit->Chcurve = mods->Chcurve;
    if (wavclCurve) toEdit->wavclCurve = mods->wavclCurve;
    //if (enacont) toEdit->enacont = mods->enacont;
    if (expcontrast) toEdit->expcontrast = mods->expcontrast;
    if (expchroma) toEdit->expchroma = mods->expchroma;
    if (expedge) toEdit->expedge = mods->expedge;
    if (expresid) toEdit->expresid = mods->expresid;
    if (expfinal) toEdit->expfinal = mods->expfinal;
    if (exptoning) toEdit->exptoning = mods->exptoning;
    if (expnoise) toEdit->expnoise = mods->expnoise;
    for(int i = 0; i < 9; ++i) {
        if(c[i]) toEdit->c[i] = dontforceSet && options.baBehav[ADDSET_WA] ? toEdit->c[i] + mods->c[i] : mods->c[i];
        if(ch[i]) toEdit->ch[i] = dontforceSet && options.baBehav[ADDSET_WA] ? toEdit->ch[i] + mods->ch[i] : mods->ch[i];
    }
    if (skinprotect) toEdit->skinprotect = dontforceSet && options.baBehav[ADDSET_WA_SKINPROTECT] ? toEdit->skinprotect + mods->skinprotect : mods->skinprotect;
    if (hueskin) toEdit->hueskin = mods->hueskin;
    if (hueskin2) toEdit->hueskin2 = mods->hueskin2;
    if (edgesensi) toEdit->edgesensi = mods->edgesensi;
    if (edgeampli) toEdit->edgeampli = mods->edgeampli;
    if (resconH) toEdit->resconH = dontforceSet && options.baBehav[ADDSET_WA_RESCONH] ? toEdit->resconH + mods->resconH : mods->resconH;
    if (reschro) toEdit->reschro = dontforceSet && options.baBehav[ADDSET_WA_RESCHRO] ? toEdit->reschro + mods->reschro : mods->reschro;
    if (tmrs) toEdit->tmrs = dontforceSet && options.baBehav[ADDSET_WA_TMRS] ? toEdit->tmrs + mods->tmrs : mods->tmrs;
    if (gamma) toEdit->gamma = dontforceSet && options.baBehav[ADDSET_WA_GAMMA] ? toEdit->gamma + mods->gamma : mods->gamma;
    if (rescon) toEdit->rescon = dontforceSet && options.baBehav[ADDSET_WA_RESCON] ? toEdit->rescon + mods->rescon : mods->rescon;
    if (thres) toEdit->thres = dontforceSet && options.baBehav[ADDSET_WA_THRES] ? toEdit->thres + mods->thres : mods->thres;
    if (threshold) toEdit->threshold = dontforceSet && options.baBehav[ADDSET_WA_THRESHOLD] ? toEdit->threshold + mods->threshold : mods->threshold;
    if (threshold2) toEdit->threshold2 = dontforceSet && options.baBehav[ADDSET_WA_THRESHOLD2] ? toEdit->threshold2 + mods->threshold2 : mods->threshold2;
    if (edgedetect) toEdit->edgedetect = dontforceSet && options.baBehav[ADDSET_WA_EDGEDETECT] ? toEdit->edgedetect + mods->edgedetect : mods->edgedetect;
    if (edgedetectthr) toEdit->edgedetectthr = dontforceSet && options.baBehav[ADDSET_WA_EDGEDETECTTHR] ? toEdit->edgedetectthr + mods->edgedetectthr : mods->edgedetectthr;
    if (edgedetectthr2) toEdit->edgedetectthr2 = dontforceSet && options.baBehav[ADDSET_WA_EDGEDETECTTHR2] ? toEdit->edgedetectthr2 + mods->edgedetectthr2 : mods->edgedetectthr2;
    if (chro) toEdit->chro = dontforceSet && options.baBehav[ADDSET_WA_CHRO] ? toEdit->chro + mods->chro : mods->chro;
    if (chroma) toEdit->chroma = dontforceSet && options.baBehav[ADDSET_WA_CHROMA] ? toEdit->chroma + mods->chroma : mods->chroma;
    if (contrast) toEdit->contrast = dontforceSet && options.baBehav[ADDSET_WA_CONTRAST] ? toEdit->contrast + mods->contrast : mods->contrast;
    if (edgrad) toEdit->edgrad = dontforceSet && options.baBehav[ADDSET_WA_EDGRAD] ? toEdit->edgrad + mods->edgrad : mods->edgrad;
    if (edgval) toEdit->edgval = dontforceSet && options.baBehav[ADDSET_WA_EDGVAL] ? toEdit->edgval + mods->edgval : mods->edgval;
    if (strength) toEdit->strength = dontforceSet && options.baBehav[ADDSET_WA_STRENGTH] ? toEdit->strength + mods->strength : mods->strength;
    // *INDENT-ON*
}

void HSVEqualizerParamsEdited::combine (HSVEqualizerParams* toEdit, const HSVEqualizerParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (hcurve) toEdit->hcurve = mods->hcurve;
    if (scurve) toEdit->scurve = mods->scurve;
    if (vcurve) toEdit->vcurve = mods->vcurve;
    // *INDENT-ON*
}

void FilmSimulationParamsEdited::combine (FilmSimulationParams* toEdit, const FilmSimulationParams* mods, bool dontforceSet)
{
    // *INDENT-OFF*
    if (enabled) toEdit->enabled = mods->enabled;
    if (clutFilename) toEdit->clutFilename = mods->clutFilename;
    if (strength) toEdit->strength = dontforceSet && options.baBehav[ADDSET_FILMSIMULATION_STRENGTH] ? toEdit->strength + mods->strength : mods->strength;
    // *INDENT-ON*
}

void ParamsEdited::combine (rtengine::procparams::ProcParams& toEdit, const rtengine::procparams::ProcParams& mods, bool forceSet)
{

    bool dontforceSet = !forceSet;

    // *INDENT-OFF*

    toneCurve.combine(&toEdit.toneCurve, &mods.toneCurve, dontforceSet);
    labCurve.combine(&toEdit.labCurve, &mods.labCurve, dontforceSet);
    rgbCurves.combine(&toEdit.rgbCurves, &mods.rgbCurves, dontforceSet);
    colorToning.combine(&toEdit.colorToning, &mods.colorToning, dontforceSet);
    sharpenEdge.combine(&toEdit.sharpenEdge, &mods.sharpenEdge, dontforceSet);
    sharpenMicro.combine(&toEdit.sharpenMicro, &mods.sharpenMicro, dontforceSet);
    sharpening.combine(&toEdit.sharpening, &mods.sharpening, dontforceSet);
    prsharpening.combine(&toEdit.prsharpening, &mods.prsharpening, dontforceSet);
    vibrance.combine(&toEdit.vibrance, &mods.vibrance, dontforceSet);
    wb.combine(&toEdit.wb, &mods.wb, dontforceSet);
    defringe.combine(&toEdit.defringe, &mods.defringe, dontforceSet);
    colorappearance.combine(&toEdit.colorappearance, &mods.colorappearance, dontforceSet);
    impulseDenoise.combine(&toEdit.impulseDenoise, &mods.impulseDenoise, dontforceSet);
    dirpyrDenoise.combine(&toEdit.dirpyrDenoise, &mods.dirpyrDenoise, dontforceSet);
    epd.combine(&toEdit.epd, &mods.epd, dontforceSet);
    sh.combine(&toEdit.sh, &mods.sh, dontforceSet);
    crop.combine(&toEdit.crop, &mods.crop, dontforceSet);
    coarse.combine(&toEdit.coarse, &mods.coarse, dontforceSet);
    commonTrans.combine(&toEdit.commonTrans, &mods.commonTrans, dontforceSet);
    rotate.combine(&toEdit.rotate, &mods.rotate, dontforceSet);
    distortion.combine(&toEdit.distortion, &mods.distortion, dontforceSet);
    lensProf.combine(&toEdit.lensProf, &mods.lensProf, dontforceSet);
    perspective.combine(&toEdit.perspective, &mods.perspective, dontforceSet);
    gradient.combine(&toEdit.gradient, &mods.gradient, dontforceSet);
    pcvignette.combine(&toEdit.pcvignette, &mods.pcvignette, dontforceSet);
    cacorrection.combine(&toEdit.cacorrection, &mods.cacorrection, dontforceSet);
    vignetting.combine(&toEdit.vignetting, &mods.vignetting, dontforceSet);
    chmixer.combine(&toEdit.chmixer, &mods.chmixer, dontforceSet);
    blackwhite.combine(&toEdit.blackwhite, &mods.blackwhite, dontforceSet);
    resize.combine(&toEdit.resize, &mods.resize, dontforceSet);
    icm.combine(&toEdit.icm, &mods.icm, dontforceSet);
    raw.combine(&toEdit.raw, &mods.raw, dontforceSet);
    wavelet.combine(&toEdit.wavelet, &mods.wavelet, dontforceSet);
    dirpyrequalizer.combine(&toEdit.dirpyrequalizer, &mods.dirpyrequalizer, dontforceSet);
    hsvequalizer.combine(&toEdit.hsvequalizer, &mods.hsvequalizer, dontforceSet);
    filmSimulation.combine(&toEdit.filmSimulation, &mods.filmSimulation, dontforceSet);
    // *INDENT-ON*

    // Exif changes are added to the existing ones
    if (exif) {
        for (procparams::ExifPairs::const_iterator i = mods.exif.begin(); i != mods.exif.end(); i++) {
            toEdit.exif[i->first] = i->second;
        }
    }

    // IPTC changes are added to the existing ones
    if (iptc) {
        for (procparams::IPTCPairs::const_iterator i = mods.iptc.begin(); i != mods.iptc.end(); i++) {
            toEdit.iptc[i->first] = i->second;
        }
    }
}

bool ParamsEdited::isTagsSet()
{
    bool retVal = general;
    return retVal;
    //return general;
}

bool ParamsEdited::isToolSet()
{
    bool retVal = toneCurve;
    retVal |= labCurve;
    retVal |= rgbCurves;
    retVal |= colorToning;
    retVal |= sharpenEdge;
    retVal |= sharpenMicro;
    retVal |= sharpening;
    retVal |= prsharpening;
    retVal |= vibrance;
    retVal |= wb;
    retVal |= defringe;
    retVal |= colorappearance;
    retVal |= impulseDenoise;
    retVal |= dirpyrDenoise;
    retVal |= epd;
    retVal |= sh;
    retVal |= crop;
    retVal |= coarse;
    retVal |= commonTrans;
    retVal |= rotate|distortion|lensProf|perspective|gradient|pcvignette|cacorrection|vignetting|chmixer|blackwhite|resize|icm|raw|wavelet|dirpyrequalizer|hsvequalizer|filmSimulation;
    return retVal;
    //return toneCurve|labCurve|rgbCurves|colorToning|sharpenEdge|sharpenMicro|sharpening|prsharpening|vibrance|wb|defringe|colorappearance|impulseDenoise|dirpyrDenoise|epd|sh|crop|coarse|commonTrans|rotate|distortion|lensProf|perspective|gradient|pcvignette|cacorrection|vignetting|chmixer|blackwhite|resize|icm|raw|wavelet|dirpyrequalizer|hsvequalizer|filmSimulation;
}

bool ParamsEdited::isExifSet()
{
    bool retVal = exif;
    return retVal;
    //return exif;
}

bool ParamsEdited::isIptcSet()
{
    bool retVal = iptc;
    return retVal;
    //return iptc;
}
