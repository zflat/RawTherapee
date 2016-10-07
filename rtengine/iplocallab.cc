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
#include <cmath>
#include <glib.h>
#include <glibmm.h>

#include "rtengine.h"
#include "improcfun.h"
#include "curves.h"
#include "gauss.h"
#include "iccstore.h"
#include "iccmatrices.h"
#include "color.h"
#include "rt_math.h"
//#define BENCHMARK
#include "StopWatch.h"

namespace rtengine
{

using namespace procparams;

extern const Settings* settings;

struct local_params {
    float yc, xc;
    float lx, ly;
    float lxL, lyT;
    int chro, cont, ligh, sens, sensh;
    double rad;
    double stren;
    int trans;
    bool inv;
    bool invrad;
    bool invret;
    float str;
};

static void calcLocalParams(int **dataspot, int oW, int oH, const LocallabParams& locallab, struct local_params& lp)
{
    int w = oW;
    int h = oH;
    double local_x = locallab.locX / 2000.0;
    double local_y = locallab.locY / 2000.0;
    double local_xL = locallab.locXL / 2000.0;
    double local_yT = locallab.locYT / 2000.0;
    double local_center_x = locallab.centerX / 2000.0 + 0.5;
    double local_center_y = locallab.centerY / 2000.0 + 0.5;
    bool inverse = locallab.invers;

    int local_chroma = locallab.chroma;
    int local_sensi = locallab.sensi;
    int local_sensih = locallab.sensih;
    int local_contrast = locallab.contrast;
    int local_lightness = locallab.lightness;
    int local_transit = locallab.transit;
    double radius = locallab.radius;
    bool inverserad = locallab.inversrad;
    bool inverseret = locallab.inversret;
    double strength = locallab.strength;
    float str = (float)locallab.str;
    lp.xc = w * local_center_x;
    lp.yc = h * local_center_y;
    lp.lx = w * local_x;
    lp.ly = h * local_y;
    lp.lxL = w * local_xL;
    lp.lyT = h * local_yT;
    lp.chro = local_chroma;
    lp.sens = local_sensi;
    lp.sensh = local_sensih;
    lp.cont = local_contrast;
    lp.ligh = local_lightness;
    lp.trans = local_transit;
    lp.rad = radius;
    lp.stren = strength;
    lp.inv = inverse;
    lp.invrad = inverserad;
    lp.invret = inverseret;
    lp.str = str;
}

inline static float calcLocalFactor(const float lox, const float loy, const float lcx, const float dx, const float lcy, const float dy, const float ach)
{
//elipse x2/a2 + y2/b2=1
//transition elipsoidal
//x==>lox y==>loy
// a==> dx  b==>dy

    float kelip = dx / dy;
    float belip = sqrt((SQR((lox - lcx) / kelip) + SQR(loy - lcy))); //determine position ellipse ==> a and b
    float aelip = belip * kelip;
    float degrad = aelip / dx;
    float ap = M_PI / (1.f - ach);
    float bp = M_PI - ap;
    return 0.5f * (1.f + xcosf(degrad * ap + bp)); //trigo cos transition

}


static void calcTransition (const float lox, const float loy, const float ach, const local_params& lp, int &zone, float &localFactor)
{
    // returns the zone (0 = outside selection, 1 = transition zone between outside and inside selection, 2 = inside selection)
    // and a factor to calculate the transition in case zone == 1

    zone = 0;

    if(lox >= lp.xc && lox < (lp.xc + lp.lx) && loy >= lp.yc && loy < lp.yc + lp.ly) {
        zone = ( (SQR(lox - lp.xc) / SQR(ach * lp.lx) + SQR(loy - lp.yc) / SQR(ach * lp.ly)) < 1.f) ? 2 : 0;

        if(!zone) {
            zone = (((SQR(lox - lp.xc) / SQR(ach * lp.lx) + SQR(loy - lp.yc) / SQR(ach * lp.ly)) > 1.f) && ((SQR(lox - lp.xc) / SQR(lp.lx) + SQR(loy - lp.yc) / SQR(lp.ly)) < 1.f)) ? 1 : 0;

            if(zone) {
                localFactor = calcLocalFactor(lox, loy, lp.xc, lp.lx, lp.yc, lp.ly, ach);
            }
        }
    } else if(lox >= lp.xc && lox < lp.xc + lp.lx && loy < lp.yc && loy > lp.yc - lp.lyT) {
        zone = (SQR(lox - lp.xc) / SQR(ach * lp.lx) + SQR(loy - lp.yc) / SQR(ach * lp.lyT)) < 1.f ? 2 : 0;

        if(!zone) {
            zone = (((SQR(lox - lp.xc) / SQR(ach * lp.lx) + SQR(loy - lp.yc) / SQR(ach * lp.lyT)) > 1.f) && ((SQR(lox - lp.xc) / SQR(lp.lx) + SQR(loy - lp.yc) / SQR(lp.lyT)) < 1.f)) ? 1 : 0;

            if(zone) {
                localFactor = calcLocalFactor(lox, loy, lp.xc, lp.lx, lp.yc, lp.lyT, ach);
            }
        }
    } else if(lox < lp.xc && lox > lp.xc - lp.lxL && loy <= lp.yc && loy > lp.yc - lp.lyT) {
        zone = (SQR(lox - lp.xc) / SQR(ach * lp.lxL) + SQR(loy - lp.yc) / SQR(ach * lp.lyT)) < 1.f ? 2 : 0;

        if(!zone) {
            zone = (((SQR(lox - lp.xc) / SQR(ach * lp.lxL) + SQR(loy - lp.yc) / SQR(ach * lp.lyT)) > 1.f) && ((SQR(lox - lp.xc) / SQR(lp.lxL) + SQR(loy - lp.yc) / SQR(lp.lyT)) < 1.f)) ? 1 : 0;

            if(zone) {
                localFactor = calcLocalFactor(lox, loy, lp.xc, lp.lxL, lp.yc, lp.lyT, ach);
            }
        }
    } else if(lox < lp.xc && lox > lp.xc - lp.lxL && loy > lp.yc && loy < lp.yc + lp.ly) {
        zone = (SQR(lox - lp.xc) / SQR(ach * lp.lxL) + SQR(loy - lp.yc) / SQR(ach * lp.ly)) < 1.f ? 2 : 0;

        if(!zone) {
            zone = (((SQR(lox - lp.xc) / SQR(ach * lp.lxL) + SQR(loy - lp.yc) / SQR(ach * lp.ly)) > 1.f) && ((SQR(lox - lp.xc) / SQR(lp.lxL) + SQR(loy - lp.yc) / SQR(lp.ly)) < 1.f)) ? 1 : 0;

            if(zone) {
                localFactor = calcLocalFactor(lox, loy, lp.xc, lp.lxL, lp.yc, lp.ly, ach);
            }
        }
    }
}


void ImProcFunctions::addGaNoise (LabImage *lab, LabImage *dst, const float mean, const float variance, const int sk)
{
    BENCHFUN
//Box-Muller method.
// add luma noise to image

    srand(1);

    const float variaFactor = SQR(variance) / sk;
    const float randFactor = 1.f / RAND_MAX;
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        float z0, z1;
        bool generate = false;
#ifdef _OPENMP
        #pragma omp for schedule(static) // static scheduling is important to avoid artefacts
#endif

        for (int y = 0; y < lab->H; y++) {
            for (int x = 0; x < lab->W; x++) {
                generate = !generate;
                float kvar = 1.f;

                if(lab->L[y][x] < 12000.f) {
                    constexpr float ah = -0.5f / 12000.f;
                    constexpr float bh = 1.5f;
                    kvar = ah * lab->L[y][x] + bh;    //increase effect for low lights < 12000.f
                } else if(lab->L[y][x] > 20000.f) {
                    constexpr float ah = -0.5f / 12768.f;
                    constexpr float bh = 1.f - 20000.f * ah;
                    kvar = ah * lab->L[y][x] + bh;    //decrease effect for high lights > 20000.f
                    kvar = kvar < 0.5f ? 0.5f : kvar;
                }

                float varia = SQR(kvar) * variaFactor;

                if(!generate) {
                    dst->L[y][x] = LIM(lab->L[y][x] + mean + varia * z1, 0.f, 32768.f);
                    continue;
                }

                int u1 = 0;
                int u2;

                while(u1 == 0) {
                    u1 = rand();
                    u2 = rand();
                }

                float u1f = u1 * randFactor;
                float u2f = u2 * randFactor;

                float2 sincosval = xsincosf(2.f * M_PI * u2f);
                float factor = sqrtf(-2.f * xlogf(u1f));
                z0 = factor * sincosval.y;
                z1 = factor * sincosval.x;

                dst->L[y][x] = LIM(lab->L[y][x] + mean + varia * z0, 0.f, 32768.f);

            }
        }
    }
}

void ImProcFunctions::BlurNoise_Local(const struct local_params& lp, LabImage* original, LabImage* transformed, const LabImage* const tmp1, int cx, int cy)
{
//local blur and noise
    BENCHFUN

    const float ach = (float)lp.trans / 100.f;

    #pragma omp parallel for schedule(dynamic,16) if (multiThread)

    for (int y = 0; y < transformed->H; y++) {
        int loy = cy + y;

        for (int x = 0; x < transformed->W; x++) {
            int lox = cx + x;

            int zone;
            float localFactor;
            calcTransition(lox, loy, ach, lp, zone, localFactor);

            switch(zone) {
                case 0: { // outside selection and outside transition zone => no effect, keep original values
                    transformed->L[y][x] = original->L[y][x];
                    transformed->a[y][x] = original->a[y][x];
                    transformed->b[y][x] = original->b[y][x];
                    break;
                }

                case 1: { // inside transition zone
                    float factorx = localFactor;
                    float difL = tmp1->L[y][x] - original->L[y][x];
                    float difa = tmp1->a[y][x] - original->a[y][x];
                    float difb = tmp1->b[y][x] - original->b[y][x];
                    difL *= factorx;
                    difa *= factorx;
                    difb *= factorx;
                    transformed->L[y][x] = original->L[y][x] + difL;
                    transformed->a[y][x] = original->a[y][x] + difa;
                    transformed->b[y][x] = original->b[y][x] + difb;
                    break;
                }

                case 2: { // inside selection => full effect, no transition
                    transformed->L[y][x] = tmp1->L[y][x];
                    transformed->a[y][x] = tmp1->a[y][x];
                    transformed->b[y][x] = tmp1->b[y][x];
                }
            }
        }
    }
}

void ImProcFunctions::InverseReti_Local(const struct local_params& lp, LabImage* original, LabImage* transformed, const LabImage* const tmp1, int cx, int cy, int chro)
{
    BENCHFUN
//inverse local retinex
    float ach = (float)lp.trans / 100.f;

    #pragma omp parallel for schedule(dynamic,16) if (multiThread)

    for (int y = 0; y < transformed->H; y++) {
        int loy = cy + y;

        for (int x = 0; x < transformed->W; x++) {
            int lox = cx + x;

            int zone;
            float localFactor;
            calcTransition (lox, loy, ach, lp, zone, localFactor);

            switch(zone) {
                case 0: { // outside selection and outside transition zone => full effect, no transition
                    if(chro == 0) {
                        transformed->L[y][x] = tmp1->L[y][x];
                    }

                    if(chro == 1) {

                        transformed->a[y][x] = tmp1->a[y][x];
                        transformed->b[y][x] = tmp1->b[y][x];
                    }

                    break;
                }

                case 1: { // inside transition zone
                    float factorx = 1.f - localFactor;

                    if(chro == 0) {
                        float difL = tmp1->L[y][x] - original->L[y][x];
                        difL *= factorx;
                        transformed->L[y][x] = original->L[y][x] + difL;
                    }

                    if(chro == 1) {
                        float difa = tmp1->a[y][x] - original->a[y][x];
                        float difb = tmp1->b[y][x] - original->b[y][x];

                        difa *= factorx;
                        difb *= factorx;

                        transformed->a[y][x] = original->a[y][x] + difa;
                        transformed->b[y][x] = original->b[y][x] + difb;
                    }

                    break;
                }

                case 2: { // inside selection => no effect, keep original values
                    if(chro == 0) {
                        transformed->L[y][x] = original->L[y][x];
                    }

                    if(chro == 1) {
                        transformed->a[y][x] = original->a[y][x];
                        transformed->b[y][x] = original->b[y][x];
                    }
                }
            }
        }
    }
}



void ImProcFunctions::Reti_Local(const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const float lumaref, const struct local_params& lp, LabImage* original, LabImage* transformed, const LabImage* const tmp1, int cx, int cy, int chro)
{

//local retinex
    BENCHFUN

    const float ach = (float)lp.trans / 100.f;

    //chroma
    constexpr float amplchsens = 2.5f;
    constexpr float achsens = (amplchsens - 1.f) / (100.f - 20.f); //20. default locallab.sensih
    constexpr float bchsens = 1.f - 20.f * achsens;
    const float multchro = lp.sensh * achsens + bchsens;

    //luma
    constexpr float ampllumsens = 2.f;
    constexpr float alumsens = (ampllumsens - 1.f) / (100.f - 20.f); //20. default locallab.sensih
    constexpr float blumsens = 1.f - 20.f * alumsens;
    const float multlum = lp.sensh * alumsens + blumsens;

    //skin
    constexpr float amplchsensskin = 1.6f;
    constexpr float achsensskin = (amplchsensskin - 1.f) / (100.f - 20.f); //20. default locallab.sensih
    constexpr float bchsensskin = 1.f - 20.f * achsensskin;
    const float multchroskin = lp.sensh * achsensskin + bchsensskin;

    //transition = difficult to avoid artifact with scope on flat area (sky...)
    float strn = lp.str / 1.f;  // we can chnage 1.f by 2 or...to chnage effect

    constexpr float delhu = 0.1f; //between 0.05 and 0.2
    const float aplus = (1.f - strn) / delhu;
    const float bplus = 1.f - aplus * hueplus;
    const float amoins = (strn - 1.f) / delhu;
    const float bmoins = 1.f - amoins * huemoins;



#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef __SSE2__
        float atan2Buffer[transformed->W] ALIGNED16;
        float sqrtBuffer[transformed->W] ALIGNED16;
        vfloat c327d68v = F2V(327.68f);
#endif

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int y = 0; y < transformed->H; y++) {
#ifdef __SSE2__
            int i = 0;

            for(; i < transformed->W - 3; i += 4) {
                vfloat av = LVFU(original->a[y][i]);
                vfloat bv = LVFU(original->b[y][i]);
                STVF(atan2Buffer[i], xatan2f(bv, av));
                STVF(sqrtBuffer[i], _mm_sqrt_ps(SQRV(bv) + SQRV(av)) / c327d68v);
            }

            for(; i < transformed->W; i++) {
                atan2Buffer[i] = xatan2f(original->b[y][i], original->a[y][i]);
                sqrtBuffer[i] = sqrt(SQR(original->b[y][i]) + SQR(original->a[y][i])) / 327.68f;
            }

#endif
            int loy = cy + y;

            for (int x = 0; x < transformed->W; x++) {
                int lox = cx + x;
#ifdef __SSE2__
                float rhue = atan2Buffer[x];
                float rchro = sqrtBuffer[x];
#else
                float rhue = xatan2f(original->b[y][x], original->a[y][x]);
                float rchro = sqrt(SQR(original->b[y][x]) + SQR(original->a[y][x])) / 327.68f;
#endif
                //    float rL = original->L[y][x] / 327.68f;

                float realstr = 1.f;
                float realstrch = 1.f;

                bool kzon = false;
                //transition = difficult to avoid artifact with scope on flat area (sky...)

                if((hueref + dhue) < M_PI && rhue < hueplus && rhue > huemoins) {//transition are good
                    if(rhue >= hueplus - delhu)  {
                        realstr = aplus * rhue + bplus;
                    } else if(rhue < huemoins + delhu)  {
                        realstr = amoins * rhue + bmoins;
                    } else {
                        realstr = strn;
                    }

                    kzon = true;
                } else if((hueref + dhue) >= M_PI && (rhue > huemoins  || rhue < hueplus )) {
                    if(rhue >= hueplus - delhu  && rhue < hueplus)  {
                        realstr = aplus * rhue + bplus;
                    } else if(rhue >= huemoins && rhue < huemoins + delhu)  {
                        realstr = amoins * rhue + bmoins;
                    } else {
                        realstr = strn;
                    }

                    kzon = true;
                }

                if((hueref - dhue) > -M_PI && rhue < hueplus && rhue > huemoins) {
                    if(rhue >= hueplus - delhu  && rhue < hueplus)  {
                        realstr = aplus * rhue + bplus;
                    } else if(rhue >= huemoins && rhue < huemoins + delhu)  {
                        realstr = amoins * rhue + bmoins;
                    } else {
                        realstr = strn;
                    }

                    kzon = true;
                } else if((hueref - dhue) <= -M_PI && (rhue > huemoins  || rhue < hueplus )) {
                    if(rhue >= hueplus - delhu  && rhue < hueplus)  {
                        realstr = aplus * rhue + bplus;
                    } else if(rhue >= huemoins && rhue < huemoins + delhu)  {
                        realstr = amoins * rhue + bmoins;
                    } else {
                        realstr = strn;
                    }

                    kzon = true;
                }

                //    float kLinf = rL / (100.f);
                //    float kLsup = kLinf;

                //    float kdiff = 0.f;
                // I add these functions...perhaps not good
                if(kzon) {
                    if(lp.sensh < 60.f) { //arbitrary value
                        if(hueref < -1.1f && hueref > -2.8f) {// detect blue sky
                            if(chromaref > 0.f && chromaref < 35.f * multchro) { // detect blue sky
                                if( (rhue > -2.79f && rhue < -1.11f) && (rchro < 35.f * multchro)) {
                                    realstr *= 0.9f;
                                } else {
                                    realstr = 1.f;
                                }
                            }
                        } else {
                            realstr = strn;
                        }

                        if(lp.sensh < 50.f) {//&& lp.chro > 0.f
                            if(hueref > -0.1f && hueref < 1.6f) {// detect skin
                                if(chromaref > 0.f && chromaref < 55.f * multchroskin) { // detect skin
                                    if( (rhue > -0.09f && rhue < 1.59f) && (rchro < 55.f * multchroskin)) {
                                        realstr *= 0.7f;
                                    } else {
                                        realstr = 1.f;
                                    }
                                }
                            } else {
                                realstr = strn;
                            }
                        }
                    }

                }

                int zone;
                float localFactor;
                calcTransition (lox, loy, ach, lp, zone, localFactor);

                switch(zone) {
                    case 0: { // outside selection and outside transition zone => no effect, keep original values
                        if(chro == 0) {
                            transformed->L[y][x] = original->L[y][x];
                        }

                        if(chro == 1) {
                            transformed->a[y][x] = original->a[y][x];
                            transformed->b[y][x] = original->b[y][x];
                        }

                        break;
                    }

                    case 1: { // inside transition zone
                        float factorx = localFactor;

                        if(chro == 0) {
                            float difL = (tmp1->L[y][x]) - original->L[y][x];
                            difL *= factorx * (100.f + realstr * (1.f - factorx)) / 100.f;
                            transformed->L[y][x] = original->L[y][x] + difL;
                        }

                        if(chro == 1) {
                            float difa = tmp1->a[y][x] - original->a[y][x];
                            float difb = tmp1->b[y][x] - original->b[y][x];
                            difa *= factorx * (100.f + realstr * (1.f - factorx)) / 100.f;
                            difb *= factorx * (100.f + realstr * (1.f - factorx)) / 100.f;
                            transformed->a[y][x] = original->a[y][x] + difa;
                            transformed->b[y][x] = original->b[y][x] + difb;
                        }

                        break;

                    }

                    case 2: { // inside selection => full effect, no transition
                        if(chro == 0) {
                            transformed->L[y][x] = tmp1->L[y][x];
                        }

                        if(chro == 1) {
                            transformed->a[y][x] = tmp1->a[y][x];
                            transformed->b[y][x] = tmp1->b[y][x];
                        }
                    }
                }
            }
        }
    }


}

void ImProcFunctions::InverseBlurNoise_Local(const struct local_params& lp, LabImage* original, LabImage* transformed, const LabImage* const tmp1, int cx, int cy)
{
    BENCHFUN
//inverse local blur and noise
    float ach = (float)lp.trans / 100.f;

    #pragma omp parallel for schedule(dynamic,16) if (multiThread)

    for (int y = 0; y < transformed->H; y++) {
        int loy = cy + y;

        for (int x = 0; x < transformed->W; x++) {
            int lox = cx + x;

            int zone;
            float localFactor;
            calcTransition (lox, loy, ach, lp, zone, localFactor);

            switch(zone) {
                case 0: { // outside selection and outside transition zone => full effect, no transition
                    transformed->L[y][x] = tmp1->L[y][x];
                    transformed->a[y][x] = tmp1->a[y][x];
                    transformed->b[y][x] = tmp1->b[y][x];
                    break;
                }

                case 1: { // inside transition zone
                    float difL = tmp1->L[y][x] - original->L[y][x];
                    float difa = tmp1->a[y][x] - original->a[y][x];
                    float difb = tmp1->b[y][x] - original->b[y][x];

                    float factorx = 1.f - localFactor;
                    difL *= factorx;
                    difa *= factorx;
                    difb *= factorx;

                    transformed->L[y][x] = original->L[y][x] + difL;
                    transformed->a[y][x] = original->a[y][x] + difa;
                    transformed->b[y][x] = original->b[y][x] + difb;
                    break;
                }

                case 2: { // inside selection => no effect, keep original values
                    transformed->L[y][x] = original->L[y][x];
                    transformed->a[y][x] = original->a[y][x];
                    transformed->b[y][x] = original->b[y][x];
                }
            }
        }
    }
}

struct local_contra {
    float alsup, blsup;
    float alsup2, blsup2;
    float alsup3, blsup3;
    float alinf;
    float aDY;
    float aa;
    float bb;
    float aaa, bbb;
    float ccc;
//    float DY;
    float dx, dy;
    float ah, bh;
    float al, bl;
};

void ImProcFunctions::Contrast_Local(const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, float pm, struct local_contra &lco, float lumaref, float av, const struct local_params& lp, LabImage* original, LabImage* transformed, int cx, int cy)
{
    BENCHFUN
// contrast - perhaps for 4 areas   if need
// I tried shmap adaptaed to Lab, but no real gain and artifacts
    const float localtype = lumaref; // always spot area
    const float ach = (float)lp.trans / 100.f;
    float reducac;

    if(lp.sens < 30.f) {
        reducac = 0.2f * (lp.sens / 100.f);
    } else {
        float areduc = 0.6285714f; //0.44f/0.7f;
        float breduc = 0.5f - areduc;
        reducac = areduc * (lp.sens / 100.f) + breduc;
    }

    const float realcox = lco.dx, realcoy = lco.dy;

    lco.alsup = (-realcox) / (localtype / 2.f);
    lco.blsup = -lco.alsup * localtype;
    lco.alsup2 = (realcoy) / (50.f - localtype / 2.f);
    lco.blsup2 = -lco.alsup2 * localtype;
    lco.alsup3 = (realcoy) / (localtype / 2.f - 50.f);
    lco.blsup3 = -lco.alsup3 * 100.f;
    lco.aDY = realcoy;



    constexpr float delhu = 0.1f; //between 0.05 and 0.2

    const float apl = (-1.f) / delhu;
    const float bpl = - apl * hueplus;
    const float amo = 1.f / delhu;
    const float bmo = - amo * huemoins;


    /*  const float aplus = (1.f - lp.sensh) / delhu;
      const float bplus = 1.f - aplus * hueplus;
      const float amoins = (lp.sensh - 1.f) / delhu;
      const float bmoins = 1.f - amoins * huemoins;
    */
    const float pb = 4.f;
    const float pa = (1.f - pb) / 40.f;

    const float ahu = 1.f / (lp.sensh - 280.f);
    const float bhu = 1.f - ahu * lp.sensh;

    lco.alinf = realcox / (localtype / 2.f);
    const float vi = (localtype / 2.f) / 100.f;
    const float vinf = (50.f + localtype / 2.f) / 100.f;
    ImProcFunctions::secondeg_begin (reducac, vi, lco.aa, lco.bb);//parabolic
    ImProcFunctions::secondeg_end (reducac, vinf, lco.aaa, lco.bbb, lco.ccc);//parabolic

    printf("huref=%f huplus=%f huemoins=%f dhue=%f\n", hueref, hueplus, huemoins, dhue);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif

    for (int y = 0; y < transformed->H; y++) {
        int loy = cy + y;

        for (int x = 0; x < transformed->W; x++) {
            int lox = cx + x;

            int zone;
            float localFactor;
            calcTransition (lox, loy, ach, lp, zone, localFactor);
            float khu = 0.f;
            float kch = 1.f;
            bool kzon = false;
            float fach = 1.f;

            //   float rhue = xatan2f(original->b[y][x], original->a[y][x]);
            float rchro = sqrt(SQR(original->b[y][x]) + SQR(original->a[y][x])) / 327.68f;
            float deltachro = fabs(rchro - chromaref);
            //   float deltahue = fabs(rhue - hueref);
            //   float deltaE = 20.f * deltahue + deltachro; //between 0 and 280

            if(deltachro < 160.f * SQR(lp.sens / 100.f)) {
                kch = 1.f;
            } else {
                float ck = 160.f * SQR(lp.sens / 100.f);
                float ak = 1.f / (ck - 160.f);
                float bk = -160.f * ak;
                kch = ak * deltachro + bk;
            }

            if(lp.sens < 40.f ) {
                kch = pow(kch, pa * lp.sens + pb);    //increase under 40
            }


            /*    // algo with detection of hue ==> artifacts
                        if(lp.sensh >= 8.f) {//to try...
                                //if(deltaE < lp.sensh) khu = 1.f;
                            //  else khu = ahu*deltaE + bhu;

                            if((hueref + dhue) < M_PI && rhue < hueplus && rhue > huemoins) {//transition are good
                                if(rhue >= hueplus - delhu )  {
                                    khu  = apl * rhue + bpl;
                                } else if(rhue < huemoins + delhu)  {
                                    khu = amo * rhue + bmo;
                                } else {
                                    khu = 1.f;
                                }


                                kzon = true;
                            } else if((hueref + dhue) >= M_PI && (rhue > huemoins  || rhue < hueplus )) {
                                if(rhue >= hueplus - delhu  && rhue < hueplus)  {
                                    khu  = apl * rhue + bpl;
                                } else if(rhue >= huemoins && rhue < huemoins + delhu)  {
                                    khu = amo * rhue + bmo;
                                } else {
                                    khu = 1.f;
                                }

                                kzon = true;
                            }

                            if((hueref - dhue) > -M_PI && rhue < hueplus && rhue > huemoins ) {
                                if(rhue >= hueplus - delhu  && rhue < hueplus)  {
                                    khu  = apl * rhue + bpl;
                                } else if(rhue >= huemoins && rhue < huemoins + delhu)  {
                                    khu = amo * rhue + bmo;
                                } else {
                                    khu = 1.f;
                                }

                                kzon = true;
                            } else if((hueref - dhue) <= -M_PI && (rhue > huemoins  || rhue < hueplus )) {
                                if(rhue >= hueplus - delhu  && rhue < hueplus)  {
                                    khu  = apl * rhue + bpl;
                                } else if(rhue >= huemoins && rhue < huemoins + delhu)  {
                                     khu = amo * rhue + bmo;
                               } else {
                                    khu = 1.f;
                                }

                                kzon = true;
                            }
                                    if(deltaE < lp.sensh) fach = khu;
                                    else fach = khu*(ahu*deltaE + bhu);

                            //fach = khu ;

                        }
            */
            switch(zone) {
                case 0: { // outside selection and outside transition zone => no effect, keep original values
                    transformed->L[y][x] = original->L[y][x];
                    transformed->a[y][x] = original->a[y][x];
                    transformed->b[y][x] = original->b[y][x];
                    break;
                }

                case 1: { // inside transition zone
                    if(original->L[y][x] < 32768.f) {
                        float factorx = localFactor;
                        float prov100 = original->L[y][x] / 32768.f;
                        float prov = prov100 * 100.f;
                        bool contin = true;


                        if(contin) {



                            if(prov > localtype) {
                                if(prov >= localtype && prov < 50.f + localtype / 2.f) {
                                    float core = (lco.alsup2 * prov + lco.blsup2) ;
                                    core *= factorx;

                                    transformed->L[y][x] = 327.68f * (prov + pm * (prov - localtype) * (core) * kch * fach);
                                } else {
                                    float core = lco.aDY * (lco.aaa * prov100 * prov100 + lco.bbb * prov100 + lco.ccc);

                                    core *= factorx;
                                    transformed->L[y][x] = 327.68f * (prov + pm * (prov - localtype) * (core) * kch * fach);
                                }
                            } else  { //inferior
                                if(2.f * prov > localtype && prov < localtype)  {
                                    float core = (lco.alsup * prov + lco.blsup) ;
                                    core *= factorx;
                                    transformed->L[y][x] = 327.68f * (prov - pm * (localtype - prov) * core * kch * fach);
                                } else if(2.f * prov <= localtype) {
                                    float core = prov * lco.alinf * (lco.aa * prov100 * prov100 + lco.bb * prov100);

                                    core *= factorx;
                                    transformed->L[y][x] = 327.68f * (prov - pm * (localtype - prov) * core * kch * fach);
                                }
                            }
                        }
                    } else {
                        transformed->L[y][x] = original->L[y][x];
                    }

                    break;
                }

                case 2: { // inside selection => full effect, no transition
                    if(original->L[y][x] < 32768.f) {
                        float prov100 = original->L[y][x] / 32768.f;
                        float prov = prov100 * 100.f;

                        bool contin = true;

                        if(contin) {



                            if(prov > localtype  ) {
                                if(prov >= localtype && prov < 50.f + localtype / 2.f) {
                                    float core = (lco.alsup2 * prov + lco.blsup2) ;
                                    transformed->L[y][x] = 327.68f * (prov + pm * (prov - localtype) * core * kch * fach);
                                } else {
                                    float core = lco.aDY * (lco.aaa * prov100 * prov100 + lco.bbb * prov100 + lco.ccc);
                                    transformed->L[y][x] = 327.68f * (prov + pm * (prov - localtype) * core * kch * fach);
                                }
                            } else  { //inferior
                                if(2.f * prov > localtype && prov < localtype)  {
                                    float core = (lco.alsup * prov + lco.blsup) ;
                                    transformed->L[y][x] = 327.68f * (prov - pm * (localtype - prov) * core * kch * fach);
                                } else if(2.f * prov <= localtype) {
                                    float core = prov * lco.alinf * (lco.aa * prov100 * prov100 + lco.bb * prov100);
                                    transformed->L[y][x] = 327.68f * (prov - pm * (localtype - prov) * core * kch * fach);
                                }
                            }
                        }
                    } else {
                        transformed->L[y][x] = original->L[y][x];
                    }
                }
            }
        }
    }
}


void ImProcFunctions::InverseContrast_Local(float ave, const local_contra& lco, const struct local_params& lp, LabImage* original, LabImage* transformed, int cx, int cy)
{
    BENCHFUN
    float ach = (float)lp.trans / 100.f;

    #pragma omp parallel for schedule(dynamic,16) if (multiThread)

    for (int y = 0; y < transformed->H; y++) {
        int loy = cy + y;

        for (int x = 0; x < transformed->W; x++) {
            int lox = cx + x;

            int zone;
            float localFactor;
            calcTransition (lox, loy, ach, lp, zone, localFactor);

            switch(zone) {
                case 0: { // outside selection and outside transition zone => full effect, no transition
                    if(original->L[y][x] < 32768.f) {
                        float prov = original->L[y][x];

                        if(original->L[y][x] > ave) {
                            float kh = lco.ah * (original->L[y][x] / 327.68f) + lco.bh;
                            original->L[y][x] = ave + kh * (original->L[y][x] - ave);
                        } else  {
                            float kl = lco.al * (original->L[y][x] / 327.68f) + 1.f;
                            original->L[y][x] = ave - kl * (ave - original->L[y][x]);
                        }

                        float diflc = original->L[y][x] - prov;
                        transformed->L[y][x] =  prov + diflc;
                    } else {
                        transformed->L[y][x] = original->L[y][x];
                    }

                    break;
                }

                case 1: { // inside transition zone
                    if(original->L[y][x] < 32768.f) {
                        float factorx = localFactor;
                        factorx = 1.f - factorx;
                        float prov = original->L[y][x];

                        if(original->L[y][x] > ave) {
                            float kh = lco.ah * (original->L[y][x] / 327.68f) + lco.bh;
                            original->L[y][x] = ave + kh * (original->L[y][x] - ave);
                        } else  {
                            float kl = lco.al * (original->L[y][x] / 327.68f) + 1.f;
                            original->L[y][x] = ave - kl * (ave - original->L[y][x]);
                        }

                        float diflc = original->L[y][x] - prov;
                        diflc *= factorx;
                        transformed->L[y][x] =  prov + diflc;

                    } else {
                        transformed->L[y][x] = original->L[y][x];
                    }

                    break;
                }

                case 2: { // inside selection => no effect, keep original values
                    transformed->L[y][x] = original->L[y][x];
                }
            }
        }
    }
}

static void calclight (float lum, int  koef, float &lumnew)
{
    if(koef > 0) {
        lumnew = lum * (1.f + (float)koef / 100.f);

        if(lumnew > 32768.f) {
            float kc = 32768.f / lumnew;
            lumnew = lum * (1.f + kc * (float)koef / 100.f);

        }
    }

    if(koef < 0) {
        lumnew = lum * (1.f + (float)koef / 100.f);

        if(lumnew < 0.f) {
            float kc = lum / (lum - lumnew);
            lumnew = lum * (1.f + kc * (float)koef / 100.f);

        }
    }


}


void ImProcFunctions::ColorLight_Local(const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const float lumaref, const local_params& lp, LabImage* original, LabImage* transformed, int cx, int cy, const LUTf & localcurve)
{
    BENCHFUN
// chroma and lightness

    const float ach = (float)lp.trans / 100.f;

    //chroma
    constexpr float amplchsens = 2.5f;
    constexpr float achsens = (amplchsens - 1.f) / (100.f - 20.f); //20. default locallab.sensi
    constexpr float bchsens = 1.f - 20.f * achsens;
    const float multchro = lp.sens * achsens + bchsens;

    //luma
    constexpr float ampllumsens = 2.f;
    constexpr float alumsens = (ampllumsens - 1.f) / (100.f - 20.f); //20. default locallab.sensi
    constexpr float blumsens = 1.f - 20.f * alumsens;
    const float multlum = lp.sens * alumsens + blumsens;

    //skin
    constexpr float amplchsensskin = 1.6f;
    constexpr float achsensskin = (amplchsensskin - 1.f) / (100.f - 20.f); //20. default locallab.sensi
    constexpr float bchsensskin = 1.f - 20.f * achsensskin;
    const float multchroskin = lp.sens * achsensskin + bchsensskin;

    //transition = difficult to avoid artifact with scope on flat area (sky...)
    constexpr float delhu = 0.1f; //between 0.05 and 0.2
    const float aplus = (1.f - lp.chro) / delhu;
    const float bplus = 1.f - aplus * hueplus;
    const float amoins = (lp.chro - 1.f) / delhu;
    const float bmoins = 1.f - amoins * huemoins;

    const float pb = 4.f;
    const float pa = (1.f - pb) / 40.f;

    //luma
    constexpr float lumdelta = 11.f; //11
    float modlum = lumdelta * multlum;

    if(lumaref + modlum >= 100.f) {
        modlum = (100.f - lumaref) / 2.f;
    }

    if(lumaref - modlum <= 0.f) {
        modlum = (lumaref) / 2.f;
    }

    float alu = 1.f / (lumaref + modlum - 100.f); //linear
    float aa, bb, aaa, bbb, ccc;
    float reducac = 0.05f;

    float vinf = (lumaref + modlum) / 100.f;
    float vi = (lumaref - modlum) / 100.f;
    ImProcFunctions::secondeg_begin (reducac, vi, aa, bb);//parabolic
    ImProcFunctions::secondeg_end (reducac, vinf, aaa, bbb, ccc);//parabolic

    float vinf2 = (lumaref + modlum) / 100.f;
    float vi2 = (lumaref - modlum) / 100.f;
    float aaaa, bbbb, cccc, aO, bO;
    ImProcFunctions::secondeg_end (0.0001f, vinf2, aaaa, bbbb, cccc);//parabolic
    ImProcFunctions::secondeg_begin (0.0001f, vi2, aO, bO);//parabolic

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef __SSE2__
        float atan2Buffer[transformed->W] ALIGNED16;
        float sqrtBuffer[transformed->W] ALIGNED16;
        vfloat c327d68v = F2V(327.68f);
#endif

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int y = 0; y < transformed->H; y++) {
#ifdef __SSE2__
            int i = 0;

            for(; i < transformed->W - 3; i += 4) {
                vfloat av = LVFU(original->a[y][i]);
                vfloat bv = LVFU(original->b[y][i]);
                STVF(atan2Buffer[i], xatan2f(bv, av));
                STVF(sqrtBuffer[i], _mm_sqrt_ps(SQRV(bv) + SQRV(av)) / c327d68v);
            }

            for(; i < transformed->W; i++) {
                atan2Buffer[i] = xatan2f(original->b[y][i], original->a[y][i]);
                sqrtBuffer[i] = sqrt(SQR(original->b[y][i]) + SQR(original->a[y][i])) / 327.68f;
            }

#endif

            int loy = cy + y;

            for (int x = 0; x < transformed->W; x++) {
                int lox = cx + x;
#ifdef __SSE2__
                float rhue = atan2Buffer[x];
                float rchro = sqrtBuffer[x];
#else
                float rhue = xatan2f(original->b[y][x], original->a[y][x]);
                float rchro = sqrt(SQR(original->b[y][x]) + SQR(original->a[y][x])) / 327.68f;
#endif
                float rL = original->L[y][x] / 327.68f;

                float realchro = 1.f;
                float deltachro = fabs(rchro - chromaref);
                float deltahue = fabs(rhue - hueref);
                float kch = 1.f;
                float khu = 1.f;


                if(deltachro < 160.f * SQR(lp.sens / 100.f)) {
                    kch = 1.f;
                } else {
                    float ck = 160.f * SQR(lp.sens / 100.f);
                    float ak = 1.f / (ck - 160.f);
                    float bk = -160.f * ak;
                    kch = ak * deltachro + bk;
                }

                if(lp.sens < 40.f ) {
                    kch = pow(kch, pa * lp.sens + pb);    //increase under 40
                }


                bool kzon = false;
                //transition = difficult to avoid artifact with scope on flat area (sky...)

                if((hueref + dhue) < M_PI && rhue < hueplus && rhue > huemoins) {//transition are good
                    if(rhue >= hueplus - delhu)  {
                        realchro = aplus * rhue + bplus;

                    } else if(rhue < huemoins + delhu)  {
                        realchro = amoins * rhue + bmoins;

                    } else {
                        realchro = lp.chro;

                    }

                    kzon = true;
                } else if((hueref + dhue) >= M_PI && (rhue > huemoins  || rhue < hueplus )) {
                    if(rhue >= hueplus - delhu  && rhue < hueplus)  {
                        realchro = aplus * rhue + bplus;

                    } else if(rhue >= huemoins && rhue < huemoins + delhu)  {
                        realchro = amoins * rhue + bmoins;

                    } else {
                        realchro = lp.chro;

                    }

                    kzon = true;
                }

                if((hueref - dhue) > -M_PI && rhue < hueplus && rhue > huemoins) {
                    if(rhue >= hueplus - delhu  && rhue < hueplus)  {
                        realchro = aplus * rhue + bplus;

                    } else if(rhue >= huemoins && rhue < huemoins + delhu)  {
                        realchro = amoins * rhue + bmoins;

                    } else {
                        realchro = lp.chro;

                    }

                    kzon = true;
                } else if((hueref - dhue) <= -M_PI && (rhue > huemoins  || rhue < hueplus )) {
                    if(rhue >= hueplus - delhu  && rhue < hueplus)  {
                        realchro = aplus * rhue + bplus;

                    } else if(rhue >= huemoins && rhue < huemoins + delhu)  {
                        realchro = amoins * rhue + bmoins;

                    } else {
                        realchro = lp.chro;
                    }

                    kzon = true;
                }

                float fach = 1.f; //khu / lp.sensh;

                if(kzon) {
                    if(lp.sens < 60.f) { //arbitrary value
                        if(hueref < -1.1f && hueref > -2.8f) {// detect blue sky
                            if(chromaref > 0.f && chromaref < 35.f * multchro) { // detect blue sky
                                if( (rhue > -2.79f && rhue < -1.11f) && (rchro < 35.f * multchro)) {
                                    realchro *= 0.9f;
                                } else {
                                    realchro = 1.f;
                                }
                            }
                        } else {
                            realchro = lp.chro;
                        }

                        if(lp.sens < 50.f && lp.chro > 0.f) {
                            if(hueref > -0.1f && hueref < 1.6f) {// detect skin
                                if(chromaref > 0.f && chromaref < 55.f * multchroskin) { // detect skin
                                    if( (rhue > -0.09f && rhue < 1.59f) && (rchro < 55.f * multchroskin)) {
                                        realchro *= 0.9f;
                                    } else {
                                        realchro = 1.f;
                                    }
                                }
                            } else {
                                realchro = lp.chro;
                            }
                        }
                    }

                }

                float kLinf = rL / (100.f);
                float kLsup = kLinf;

                float kdiff = 0.f;

                if(kzon/*rhue < hueplus && rhue > huemoins*/) {

                    if( (rL > (lumaref - modlum) && rL < (lumaref + modlum))) {
                        kdiff = 1.f;
                    } else if (rL > 0.f && rL <= (lumaref - modlum)) {
                        kdiff = aa * kLinf * kLinf + bb * kLinf;    //parabolic
                    } else if (rL <= 100.f && rL >= (lumaref + modlum)) {
                        kdiff = aaa * kLsup * kLsup + bbb * kLsup + ccc;    //parabolic
                    }

                    //end luma
                } else {
                    if( (rL > (lumaref - modlum) && rL < (lumaref + modlum))) {
                        kdiff = 0.9f;
                    } else if (rL > 0.f && rL <= (lumaref - modlum)) {
                        kdiff = 0.9f * (aO * kLinf * kLinf + bO * kLinf);    //parabolic
                    } else if (rL <= 100.f && rL >= (lumaref + modlum)) {
                        kdiff = 0.9f * (aaaa * kLsup * kLsup + bbbb * kLsup + cccc);    //parabolic
                    }

                }

                int zone;
                float localFactor;
                calcTransition (lox, loy, ach, lp, zone, localFactor);

                switch(zone) {
                    case 0: { // outside selection and outside transition zone => no effect, keep original values
                        transformed->L[y][x] = original->L[y][x];
                        transformed->a[y][x] = original->a[y][x];
                        transformed->b[y][x] = original->b[y][x];
                        break;
                    }

                    case 1: { // inside transition zone
                        float lumnew = original->L[y][x];

                        if(lp.ligh != 0) {
                            calclight (original->L[y][x], lp.ligh , lumnew);
                        }

                        //    float lightcont = localcurve[original->L[y][x]]; //apply lightness
                        float lightcont = lumnew ; //original->L[y][x] + (lp.ligh /100.f)*original->L[y][x] ; //apply lightness
                        float factorx = localFactor;
                        float fac = (100.f + factorx * realchro) / 100.f; //chroma factor transition
                        float diflc = lightcont - original->L[y][x];
                        diflc *= kdiff * kch * khu;

                        diflc *= factorx; //transition lightess

                        transformed->L[y][x] = original->L[y][x] + diflc;
                        transformed->a[y][x] = original->a[y][x] * fac;
                        transformed->b[y][x] = original->b[y][x] * fac;
                        break;
                    }

                    case 2: { // inside selection => full effect, no transition
                        //  float lightcont = localcurve[original->L[y][x]]; //apply lightness
                        //  float lightcont = original->L[y][x] + (lp.ligh /100.f)*original->L[y][x]; //apply lightness
                        float lumnew = original->L[y][x];

                        if(lp.ligh != 0) {
                            calclight (original->L[y][x], lp.ligh , lumnew);
                        }

                        //    float lightcont = localcurve[original->L[y][x]]; //apply lightness
                        float lightcont = lumnew ; //original->L[y][x] + (lp.ligh /100.f)*original->L[y][x] ; //apply lightness

                        float fac = (100.f + realchro) / 100.f; //chroma factor transition
                        float diflc = lightcont - original->L[y][x];
                        diflc *= kdiff * kch * khu;
                        transformed->L[y][x] = original->L[y][x] + diflc;
                        transformed->a[y][x] = original->a[y][x] * fac;
                        transformed->b[y][x] = original->b[y][x] * fac;

                    }
                }
            }
        }
    }
}

void ImProcFunctions::InverseColorLight_Local(const struct local_params& lp, LabImage* original, LabImage* transformed, int cx, int cy, const LUTf & localcurve)
{
    BENCHFUN
    float ach = (float)lp.trans / 100.f;
    const float facc = (100.f + lp.chro) / 100.f; //chroma factor transition

    #pragma omp parallel for schedule(dynamic,16) if (multiThread)

    for (int y = 0; y < transformed->H; y++) {
        int loy = cy + y;

        for (int x = 0; x < transformed->W; x++) {
            int lox = cx + x;

            int zone;
            float localFactor;
            calcTransition (lox, loy, ach, lp, zone, localFactor);

            switch(zone) {
                case 0: { // outside selection and outside transition zone => no effect, keep original values
                    float lumnew = original->L[y][x];

                    if(lp.ligh != 0) {
                        calclight (original->L[y][x], lp.ligh , lumnew);
                    }

                    //    float lightcont = localcurve[original->L[y][x]]; //apply lightness
                    float lightcont = lumnew ; //original->L[y][x] + (lp.ligh /100.f)*original->L[y][x] ; //apply lightness



                    transformed->L[y][x] = lightcont; //localcurve[original->L[y][x]];  //apply lightness
                    transformed->a[y][x] = original->a[y][x] * facc;
                    transformed->b[y][x] = original->b[y][x] * facc;
                    break;
                }

                case 1: { // inside transition zone
                    float factorx = 1.f - localFactor;
                    float fac = (100.f + factorx * lp.chro) / 100.f; //chroma factor transition
                    float lumnew = original->L[y][x];

                    if(lp.ligh != 0) {
                        calclight (original->L[y][x], lp.ligh , lumnew);
                    }

                    //    float lightcont = localcurve[original->L[y][x]]; //apply lightness
                    float lightcont = lumnew ; //original->L[y][x] + (lp.ligh /100.f)*original->L[y][x] ; //apply lightness

                    //float lightcont = localcurve[original->L[y][x]]; //apply lightness
                    float diflc = lightcont - original->L[y][x];
                    diflc *= factorx;
                    transformed->L[y][x] = original->L[y][x] + diflc;
                    transformed->a[y][x] = original->a[y][x] * fac;
                    transformed->b[y][x] = original->b[y][x] * fac;
                    break;
                }

                case 2: { // inside selection => no effect, keep original values
                    transformed->L[y][x] = original->L[y][x];
                    transformed->a[y][x] = original->a[y][x];
                    transformed->b[y][x] = original->b[y][x];
                }
            }
        }
    }

}


void ImProcFunctions::Lab_Local(int **dataspot, LabImage* original, LabImage* transformed, int sx, int sy, int cx, int cy, int oW, int oH,  int fw, int fh,  LUTf & localcurve, bool locutili, int sk, const LocretigainCurve & locRETgainCcurve, double &hueref, double &chromaref, double &lumaref)
{
    if(params->locallab.enabled) {
        BENCHFUN
#ifdef _DEBUG
        MyTime t1e, t2e;
        t1e.set();
        // init variables to display Munsell corrections
        MunsellDebugInfo* MunsDebugInfo = new MunsellDebugInfo();
#endif

        struct local_params lp;
        calcLocalParams(dataspot, oW, oH, params->locallab, lp);

        const float radius = lp.rad / (sk * 2.f); //0 to 50 ==> see skip

        if(radius >= GAUSS_SKIP || lp.stren > 0.1) { // radius < GAUSS_SKIP means no gauss, just copy of original image
            LabImage *tmp1 = new LabImage(transformed->W, transformed->H);;
            int GW = transformed->W;
            int GH = transformed->H;

#ifdef _OPENMP
            #pragma omp parallel
#endif
            {
                gaussianBlur (original->L, tmp1->L, GW, GH, radius);
                gaussianBlur (original->a, tmp1->a, GW, GH, radius);
                gaussianBlur (original->b, tmp1->b, GW, GH, radius);

            }

            if(lp.stren > 0.1f) {
                float mean = 0.f;//0 best result
                float variance = lp.stren ; //(double) SQR(lp.stren)/sk;
                addGaNoise (tmp1, tmp1, mean, variance, sk) ;
            }

            if(!lp.invrad) { //blur and noise (center)
                BlurNoise_Local(lp, original, transformed, tmp1, cx, cy);
            } else {
                InverseBlurNoise_Local(lp, original, transformed, tmp1, cx, cy);
            }

            delete tmp1;
        }

        /*       if(lp.str > 0.f) {
                   int GW = transformed->W;
                   int GH = transformed->H;

                   float *orig[GH] ALIGNED16;
                   float *origBuffer = new float[GH * GW];

                   for (int i = 0; i < GH; i++) {
                       orig[i] = &origBuffer[i * GW];
                   }

                   float *orig1[GH] ALIGNED16;
                   float *origBuffer1 = new float[GH * GW];

                   for (int i = 0; i < GH; i++) {
                       orig1[i] = &origBuffer1[i * GW];
                   }


                   LabImage *tmpl = new LabImage(transformed->W, transformed->H);

        #ifdef _OPENMP
                   #pragma omp parallel for schedule(dynamic,16)
        #endif

                   for(int ir = 0; ir < GH; ir += 1)
                       for(int jr = 0; jr < GW; jr += 1) {
                           orig[ir][jr] = original->L[ir][jr];
                           orig1[ir][jr] = transformed->L[ir][jr];
                       }

                   float minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax;
                   ImProcFunctions::MSRLocal(orig, tmpl->L, orig1, GW, GH, params->locallab, sk, locRETgainCcurve, 0, 4, 0.8f, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax);
        #ifdef _OPENMP
                   #pragma omp parallel for schedule(dynamic,16)
        #endif

                   for(int ir = 0; ir < GH; ir += 1)
                       for(int jr = 0; jr < GW; jr += 1) {
                           tmpl->L[ir][jr] = orig[ir][jr];
                       }

                   if(!lp.invret) {
                       Reti_Local(lp, original, transformed, tmpl, cx, cy, 0);
                   } else {
                       InverseReti_Local(lp, original, transformed, tmpl, cx, cy, 0);
                   }

                   if(params->locallab.chrrt > 0) {
        #ifdef _OPENMP
                       #pragma omp parallel for schedule(dynamic,16)
        #endif

                       for(int ir = 0; ir < GH; ir += 1)
                           for(int jr = 0; jr < GW; jr += 1) {
                               orig[ir][jr] = sqrt(SQR(original->a[ir][jr]) + SQR(original->b[ir][jr]));
                               orig1[ir][jr] = sqrt(SQR(transformed->a[ir][jr]) + SQR(transformed->b[ir][jr]));
                           }

                       ImProcFunctions::MSRLocal(orig, tmpl->L, orig1, GW, GH, params->locallab, sk, locRETgainCcurve, 1, 4, 0.8f, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax);
        #ifdef _OPENMP
                       #pragma omp parallel for schedule(dynamic,16)
        #endif

                       for(int ir = 0; ir < GH; ir += 1)
                           for(int jr = 0; jr < GW; jr += 1) {
                               float Chprov = orig1[ir][jr];
                               float2 sincosval;
                               sincosval.y = Chprov == 0.0f ? 1.f : transformed->a[ir][jr] / Chprov;
                               sincosval.x = Chprov == 0.0f ? 0.f : transformed->b[ir][jr] / Chprov;
                               tmpl->a[ir][jr] = orig[ir][jr] * sincosval.y;
                               tmpl->b[ir][jr] = orig[ir][jr] * sincosval.x;

                           }

                       if(!lp.invret) {
                           Reti_Local(lp, original, transformed, tmpl, cx, cy, 1);
                       } else {
                           InverseReti_Local(lp, original, transformed, tmpl, cx, cy, 1);
                       }

                   }

                   delete tmpl;
                   delete [] origBuffer;
                   delete [] origBuffer1;

               }
        */
        //begin contrast and evalue hue
        // double precision for large summations
        double ave = 0.;
        double aveA = 0.;
        double aveB = 0.;
        double aveL = 0.;
        double aveChro = 0.;
        // int precision for the counters
        int n = 0;
        int nab = 0;
        // single precision for the result
        float av, avA, avB, avL;
        //evauate mean luminance for contrast : actually one area
        // evaluate also hue


        if ((!lp.inv  || !lp.invret)  && hueref == INFINITY && chromaref == INFINITY && lumaref == INFINITY) {
            //evaluate hue, chroma, luma in center spot
            int spotSize = max(1, 18 / sk);

            // very small region, don't use omp here
            for (int y = max(cy, (int)(lp.yc - spotSize)); y < min(transformed->H + cy, (int)(lp.yc + spotSize + 1)); y++) {
                for (int x = max(cx, (int)(lp.xc - spotSize)); x < min(transformed->W + cx, (int)(lp.xc + spotSize + 1)); x++) {
                    aveL += original->L[y - cy][x - cx];
                    aveA += original->a[y - cy][x - cx];
                    aveB += original->b[y - cy][x - cx];
                    aveChro += sqrtf(SQR(original->b[y - cy][x - cx]) + SQR(original->a[y - cy][x - cx]));

                    nab++;
                }
            }

        } else if (lp.inv || lp.invret) { //exterior
            ave = 0.f;
            n = 0;

            #pragma omp parallel for reduction(+:ave,n)

            for (int y = 0; y < transformed->H; y++) {
                for (int x = 0; x < transformed->W; x++) {
                    int lox = cx + x;
                    int loy = cy + y;

                    if(lox >= lp.xc && lox < lp.xc + lp.lx && loy >= lp.yc && loy < lp.yc + lp.ly) {
                    } else if(lox >= lp.xc && lox < lp.xc + lp.lx && loy < lp.yc && loy > lp.yc - lp.lyT) {
                    } else if(lox < lp.xc && lox > lp.xc - lp.lxL && loy <= lp.yc && loy > lp.yc - lp.lyT) {
                    } else if(lox < lp.xc && lox > lp.xc - lp.lxL && loy > lp.yc && loy < lp.yc + lp.ly) {
                    } else {
                        ave += original->L[y][x];
                        n++;
                    }
                }
            }

            if(n == 0) {
                ave = 15000.f;
                n = 1;
            }

            ave = ave / n;
            av = ave / 327.68f;
        }

        aveL = aveL / nab;
        aveA = aveA / nab;
        aveB = aveB / nab;
        aveChro = aveChro / nab;
        aveChro /= 327.68f;
        avA = aveA / 327.68f;
        avB = aveB / 327.68f;
        avL = aveL / 327.68f;

        if(hueref == INFINITY) {
            hueref = xatan2f(avB, avA);    //mean hue
        }

        if(chromaref == INFINITY) {
            chromaref = aveChro;
        }

        if(lumaref == INFINITY) {
            lumaref = avL;
        }

        struct local_contra lco;

        // we must here detect : general case, skin, sky,...foliages ???
        // delta dhue, luminance and chroma
        constexpr float ared = (M_PI - 0.05f) / 100.f;

        constexpr float bred = 0.05f;

        float dhue = ared * lp.sens + bred; //delta hue lght chroma

        float dhueret = ared * lp.sensh + bred; //delta hue retinex

        constexpr float maxh = 4.f; //amplification contrast above mean

        constexpr float maxl = 3.f; //reductio contrast under mean

        float multh = (float) fabs(lp.cont) * (maxh - 1.f) / 100.f + 1.f;

        float mult = (float)fabs(lp.cont) * (maxl - 1.f) / 100.f + 1.f;

        lco.dx = 1.f - 1.f / mult;

        lco.dy = 1.f - 1.f / multh;



        if(!lp.inv) {   //contrast interior ellipse
            const float pm = lp.cont < 0.f ? -1.f : 1.f;
            float hueplus = hueref + dhue;
            float huemoins = hueref - dhue;

            if(hueplus > M_PI) {
                hueplus = hueref + dhue - 2.f * M_PI;
            }

            if(huemoins < -M_PI) {
                huemoins = hueref - dhue + 2.f * M_PI;
            }

            Contrast_Local(hueplus, huemoins, hueref, dhue, chromaref, pm, lco, lumaref, av, lp, original, transformed, cx, cy);
        } else if(lp.inv) {

            float multL = (float)lp.cont * (maxl - 1.f) / 100.f + 1.f;
            float multH = (float) lp.cont * (maxh - 1.f) / 100.f + 1.f;

            lco.ah = (multH - 1.f) / (av - 100.f); //av ==> lumaref
            lco.bh = 1.f - 100.f * lco.ah;
            lco.al = (multL - 1.f) / av;
            lco.bl = 1.f;
            InverseContrast_Local(ave, lco, lp, original, transformed, cx, cy);
        }

// end contrast interior and exterior

        if(!lp.inv) { //interior ellipse renforced lightness and chroma
            float hueplus = hueref + dhue;
            float huemoins = hueref - dhue;

            if(hueplus > M_PI) {
                hueplus = hueref + dhue - 2.f * M_PI;
            }

            if(huemoins < -M_PI) {
                huemoins = hueref - dhue + 2.f * M_PI;
            }

            ColorLight_Local(hueplus, huemoins, hueref, dhue, chromaref, lumaref, lp, original, transformed, cx, cy, localcurve);
        }
        //inverse
        else if(lp.inv) {
            InverseColorLight_Local(lp, original, transformed, cx, cy, localcurve);
        }



        if(lp.str > 0.f) {
            int GW = transformed->W;
            int GH = transformed->H;
            float hueplus = hueref + dhueret;
            float huemoins = hueref - dhueret;

            if(hueplus > M_PI) {
                hueplus = hueref + dhueret - 2.f * M_PI;
            }

            if(huemoins < -M_PI) {
                huemoins = hueref - dhueret + 2.f * M_PI;
            }

            float *orig[GH] ALIGNED16;
            float *origBuffer = new float[GH * GW];

            for (int i = 0; i < GH; i++) {
                orig[i] = &origBuffer[i * GW];
            }

            float *orig1[GH] ALIGNED16;
            float *origBuffer1 = new float[GH * GW];

            for (int i = 0; i < GH; i++) {
                orig1[i] = &origBuffer1[i * GW];
            }


            LabImage *tmpl = new LabImage(transformed->W, transformed->H);

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16)
#endif

            for(int ir = 0; ir < GH; ir += 1)
                for(int jr = 0; jr < GW; jr += 1) {
                    orig[ir][jr] = original->L[ir][jr];
                    orig1[ir][jr] = transformed->L[ir][jr];
                }

            float minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax;
            ImProcFunctions::MSRLocal(orig, tmpl->L, orig1, GW, GH, params->locallab, sk, locRETgainCcurve, 0, 4, 0.8f, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax);
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16)
#endif

            for(int ir = 0; ir < GH; ir += 1)
                for(int jr = 0; jr < GW; jr += 1) {
                    tmpl->L[ir][jr] = orig[ir][jr];
                }

            if(!lp.invret) {
                Reti_Local(hueplus, huemoins, hueref, dhueret, chromaref, lumaref, lp, original, transformed, tmpl, cx, cy, 0);
            } else {
                InverseReti_Local(lp, original, transformed, tmpl, cx, cy, 0);
            }

            if(params->locallab.chrrt > 0) {
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for(int ir = 0; ir < GH; ir += 1)
                    for(int jr = 0; jr < GW; jr += 1) {
                        orig[ir][jr] = sqrt(SQR(original->a[ir][jr]) + SQR(original->b[ir][jr]));
                        orig1[ir][jr] = sqrt(SQR(transformed->a[ir][jr]) + SQR(transformed->b[ir][jr]));
                    }

                ImProcFunctions::MSRLocal(orig, tmpl->L, orig1, GW, GH, params->locallab, sk, locRETgainCcurve, 1, 4, 0.8f, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax);
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for(int ir = 0; ir < GH; ir += 1)
                    for(int jr = 0; jr < GW; jr += 1) {
                        float Chprov = orig1[ir][jr];
                        float2 sincosval;
                        sincosval.y = Chprov == 0.0f ? 1.f : transformed->a[ir][jr] / Chprov;
                        sincosval.x = Chprov == 0.0f ? 0.f : transformed->b[ir][jr] / Chprov;
                        tmpl->a[ir][jr] = orig[ir][jr] * sincosval.y;
                        tmpl->b[ir][jr] = orig[ir][jr] * sincosval.x;

                    }

                if(!lp.invret) {
                    Reti_Local(hueplus, huemoins, hueref, dhueret, chromaref, lumaref, lp, original, transformed, tmpl, cx, cy, 1);
                } else {
                    InverseReti_Local(lp, original, transformed, tmpl, cx, cy, 1);
                }

            }

            delete tmpl;
            delete [] origBuffer;
            delete [] origBuffer1;

        }


// Gamut and Munsell control
        if(params->locallab.avoid) {
            TMatrix wiprof = iccStore->workingSpaceInverseMatrix (params->icm.working);
            float wip[3][3] = {
                {static_cast<float>(wiprof[0][0]), static_cast<float>(wiprof[0][1]), static_cast<float>(wiprof[0][2])},
                {static_cast<float>(wiprof[1][0]), static_cast<float>(wiprof[1][1]), static_cast<float>(wiprof[1][2])},
                {static_cast<float>(wiprof[2][0]), static_cast<float>(wiprof[2][1]), static_cast<float>(wiprof[2][2])}
            };
            const bool highlight = params->toneCurve.hrenabled;
            const bool needHH = (lp.chro != 0.f);
#ifdef _OPENMP
            #pragma omp parallel if (multiThread)
#endif
            {
#ifdef __SSE2__
                float atan2Buffer[transformed->W] ALIGNED16;
                float sqrtBuffer[transformed->W] ALIGNED16;
                float sincosyBuffer[transformed->W] ALIGNED16;
                float sincosxBuffer[transformed->W] ALIGNED16;
                vfloat c327d68v = F2V(327.68f);
                vfloat onev = F2V(1.f);
#endif

#ifdef _OPENMP
#ifdef _DEBUG
                #pragma omp for schedule(dynamic,16) firstprivate(MunsDebugInfo)
#else
                #pragma omp for schedule(dynamic,16)
#endif
#endif

                for (int y = 0; y < transformed->H; y++) {
#ifdef __SSE2__
                    int i = 0;

                    for(; i < transformed->W - 3; i += 4) {
                        vfloat av = LVFU(transformed->a[y][i]);
                        vfloat bv = LVFU(transformed->b[y][i]);

                        if(needHH) { // only do expensive atan2 calculation if needed
                            STVF(atan2Buffer[i], xatan2f(bv, av));
                        }

                        vfloat Chprov1v = vsqrtf(SQRV(bv) + SQRV(av));
                        STVF(sqrtBuffer[i], Chprov1v / c327d68v);
                        vfloat sincosyv = av / Chprov1v;
                        vfloat sincosxv = bv / Chprov1v;
                        vmask selmask = vmaskf_eq(Chprov1v, ZEROV);
                        sincosyv = vself(selmask, onev, sincosyv);
                        sincosxv = vselfnotzero(selmask, sincosxv);
                        STVF(sincosyBuffer[i], sincosyv);
                        STVF(sincosxBuffer[i], sincosxv);
                    }

                    for(; i < transformed->W; i++) {
                        float aa = transformed->a[y][i];
                        float bb = transformed->b[y][i];

                        if(needHH) { // only do expensive atan2 calculation if needed
                            atan2Buffer[i] = xatan2f(bb, aa);
                        }

                        float Chprov1 = sqrtf(SQR(bb) + SQR(aa));
                        sqrtBuffer[i] = Chprov1 / 327.68f;

                        if(Chprov1 == 0.0f) {
                            sincosyBuffer[i] = 1.f;
                            sincosxBuffer[i] = 0.0f;
                        } else {
                            sincosyBuffer[i] = aa / Chprov1;
                            sincosxBuffer[i] = bb / Chprov1;
                        }

                    }

#endif

                    for (int x = 0; x < transformed->W; x++) {
                        float Lprov1 = transformed->L[y][x] / 327.68f;
                        float2 sincosval;
#ifdef __SSE2__
                        float HH = atan2Buffer[x]; // reading HH from line buffer even if line buffer is not filled is faster than branching
                        float Chprov1 = sqrtBuffer[x];
                        sincosval.y = sincosyBuffer[x];
                        sincosval.x = sincosxBuffer[x];
#else
                        float aa = transformed->a[y][x];
                        float bb = transformed->b[y][x];
                        float HH;

                        if(needHH) { // only do expensive atan2 calculation if needed
                            HH = xatan2f(bb, aa);
                        }

                        float Chprov1 = sqrtf(SQR(aa) + SQR(bb)) / 327.68f;

                        if(Chprov1 == 0.0f) {
                            sincosval.y = 1.f;
                            sincosval.x = 0.0f;
                        } else {
                            sincosval.y = aa / (Chprov1 * 327.68f);
                            sincosval.x = bb / (Chprov1 * 327.68f);
                        }

#endif

#ifdef _DEBUG
                        bool neg = false;
                        bool more_rgb = false;
                        Color::gamutLchonly(sincosval, Lprov1, Chprov1, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
                        Color::gamutLchonly(sincosval, Lprov1, Chprov1, wip, highlight, 0.15f, 0.96f);
#endif

                        transformed->L[y][x] = Lprov1 * 327.68f;
                        transformed->a[y][x] = 327.68f * Chprov1 * sincosval.y;
                        transformed->b[y][x] = 327.68f * Chprov1 * sincosval.x;

                        if (needHH) {
                            float Lprov2 = original->L[y][x] / 327.68f;
                            float correctionHue = 0.f; // Munsell's correction
                            float correctlum = 0.f;
                            float memChprov = sqrtf(SQR(original->a[y][x]) + SQR(original->b[y][x])) / 327.68f;
                            float Chprov = sqrtf(SQR(transformed->a[y][x]) + SQR(transformed->b[y][x])) / 327.68f;
#ifdef _DEBUG
                            Color::AllMunsellLch(true, Lprov1, Lprov2, HH, Chprov, memChprov, correctionHue, correctlum, MunsDebugInfo);
#else
                            Color::AllMunsellLch(true, Lprov1, Lprov2, HH, Chprov, memChprov, correctionHue, correctlum);
#endif

                            if(fabs(correctionHue) < 0.015f) {
                                HH += correctlum;    // correct only if correct Munsell chroma very little.
                            }

                            float2 sincosval = xsincosf(HH + correctionHue);

                            transformed->a[y][x] = 327.68f * Chprov * sincosval.y; // apply Munsell
                            transformed->b[y][x] = 327.68f * Chprov * sincosval.x;
                        }
                    }
                }
            }
        }

#ifdef _DEBUG

        if (settings->verbose) {
            t2e.set();
            printf("Color::AllMunsellLch (correction performed in %d usec):\n", t2e.etime(t1e));
            //  printf("   Munsell chrominance: MaxBP=%1.2frad MaxRY=%1.2frad MaxGY=%1.2frad MaxRP=%1.2frad  dep=%i\n", MunsDebugInfo->maxdhue[0],    MunsDebugInfo->maxdhue[1],    MunsDebugInfo->maxdhue[2],    MunsDebugInfo->maxdhue[3],    MunsDebugInfo->depass);
            //  printf("   Munsell luminance  : MaxBP=%1.2frad MaxRY=%1.2frad MaxGY=%1.2frad MaxRP=%1.2frad  dep=%i\n", MunsDebugInfo->maxdhuelum[0], MunsDebugInfo->maxdhuelum[1], MunsDebugInfo->maxdhuelum[2], MunsDebugInfo->maxdhuelum[3], MunsDebugInfo->depassLum);
        }

        delete MunsDebugInfo;
#endif

    }

}

}
