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
 *  2016 Jacques Desmis <jdesmis@gmail.com>
 *  2016 Ingo Weyrich <heckflosse@i-weyrich.de>

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
#ifdef _DEBUG
#include "mytime.h"
#endif

#include "cplx_wavelet_dec.h"

//#define BENCHMARK
//#include "StopWatch.h"

#define cliploc( val, minv, maxv )    (( val = (val < minv ? minv : val ) ) > maxv ? maxv : val )

#define CLIPC(a) ((a)>-42000?((a)<42000?(a):42000):-42000)  // limit a and b  to 130 probably enough ?
#define CLIPL(x) LIM(x,0.f,40000.f) // limit L to about L=120 probably enough ?
namespace rtengine
{
using namespace procparams;

extern const Settings* settings;

struct local_params {
    float yc, xc;
    float lx, ly;
    float lxL, lyT;
    float dxx, dyy;
    int cir;
    float thr;
    int prox;
    int chro, cont, ligh, sens, sensh;
    int shamo, shdamp, shiter, senssha;
    double shrad;
    double rad;
    double stren;
    int trans;
    bool inv;
    bool invrad;
    bool invret;
    bool invshar;
    bool actsp;
    float str;
    int qualmet;
    float noiself;
    float noiselc;
    float noisecf;
    float noisecc;
};

static void calcLocalParams(int oW, int oH, const LocallabParams& locallab, struct local_params& lp)
{
    int w = oW;
    int h = oH;
    int circr = locallab.circrad;
    float thre = locallab.thres / 100.f;
    double local_x = locallab.locX / 2000.0;
    double local_y = locallab.locY / 2000.0;
    double local_xL = locallab.locXL / 2000.0;
    double local_yT = locallab.locYT / 2000.0;
    double local_center_x = locallab.centerX / 2000.0 + 0.5;
    double local_center_y = locallab.centerY / 2000.0 + 0.5;
    double local_dxx = locallab.proxi / 8000.0;//for proxi = 2==> # 1 pixel
    double local_dyy = locallab.proxi / 8000.0;

    if(locallab.qualityMethod == "std") {
        lp.qualmet = 0;
    } else if(locallab.qualityMethod == "enh") {
        lp.qualmet = 1;
    }

    float local_noiself = locallab.noiselumf;
    float local_noiselc = locallab.noiselumc;
    float local_noisecf = locallab.noisechrof;
    float local_noisecc = locallab.noisechroc;

    int local_chroma = locallab.chroma;
    int local_sensi = locallab.sensi;
    int local_sensih = locallab.sensih;
    int local_contrast = locallab.contrast;
    int local_lightness = locallab.lightness;
    int local_transit = locallab.transit;
    double radius = locallab.radius;
    double sharradius = ((double) locallab.sharradius) / 100. ;
    int local_sensisha = locallab.sensisha;
    int local_sharamount = locallab.sharamount;
    int local_shardamping = locallab.shardamping;
    int local_shariter = locallab.shariter;
    bool inverse = locallab.invers;
    bool acti = locallab.activsp;

    bool inverserad = locallab.inversrad;
    bool inverseret = locallab.inversret;
    bool inversesha = locallab.inverssha;
    double strength = locallab.strength;
    float str = (float)locallab.str;
    lp.cir = circr;
    lp.actsp = acti;
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
    lp.invshar = inversesha;
    lp.str = str;
    lp.shrad = sharradius;
    lp.senssha = local_sensisha;
    lp.shamo = local_sharamount;
    lp.shdamp = local_shardamping;
    lp.shiter = local_shariter;
    lp.dxx = w * local_dxx;
    lp.dyy = h * local_dyy;
    lp.thr = thre;
    lp.noiself = local_noiself;
    lp.noiselc = local_noiselc;
    lp.noisecf = local_noisecf;
    lp.noisecc = local_noisecc;


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
//   BENCHFUN
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

void ImProcFunctions::DeNoise_Local(int call, const struct local_params& lp, LabImage* original, LabImage* transformed, const LabImage* const tmp1, int cx, int cy)
{
    // local denoise
    // BENCHFUN
    const float ach = (float)lp.trans / 100.f;

    #pragma omp parallel for schedule(dynamic,16) if (multiThread)

    for (int y = 0; y < transformed->H; y++) {
        int loy = cy + y;

        for (int x = 0; x < transformed->W; x++) {
            int lox = cx + x;

            int zone;
            float localFactor;
            calcTransition(lox, loy, ach, lp, zone, localFactor);
            int begx = int(lp.xc - lp.lxL);
            int begy = int(lp.yc - lp.lyT);

            switch(zone) {
                case 0: { // outside selection and outside transition zone => no effect, keep original values
                    transformed->L[y][x] = original->L[y][x];
                    transformed->a[y][x] = original->a[y][x];
                    transformed->b[y][x] = original->b[y][x];
                    break;
                }

                case 1: { // inside transition zone
                    float factorx = localFactor;
                    float difL, difa, difb;

                    if(call == 2) {//simpleprocess
                        if(lox >= (lp.xc - lp.lxL) && lox < (lp.xc + lp.lx) && loy >= (lp.yc - lp.lyT) && loy < (lp.yc + lp.ly)) {
                            difL = tmp1->L[loy - begy - 1][lox - begx - 1] - original->L[y][x];
                            difa = tmp1->a[loy - begy - 1][lox - begx - 1] - original->a[y][x];
                            difb = tmp1->b[loy - begy - 1][lox - begx - 1] - original->b[y][x];
                        }
                    } else if(call == 1) {//dcrop
                        difL = tmp1->L[y][x] - original->L[y][x];
                        difa = tmp1->a[y][x] - original->a[y][x];
                        difb = tmp1->b[y][x] - original->b[y][x];

                    }

                    difL *= factorx;
                    difa *= factorx;
                    difb *= factorx;
                    transformed->L[y][x] = original->L[y][x] + difL;
                    transformed->a[y][x] = original->a[y][x] + difa;
                    transformed->b[y][x] = original->b[y][x] + difb;
                    break;
                }

                case 2: { // inside selection => full effect, no transition
                    float difL, difa, difb;

                    if(call == 2) {//simpleprocess

                        if(lox >= (lp.xc - lp.lxL) && lox < (lp.xc + lp.lx) && loy >= (lp.yc - lp.lyT) && loy < (lp.yc + lp.ly)) {
                            difL = tmp1->L[loy - begy - 1][lox - begx - 1] - original->L[y][x];
                            difa = tmp1->a[loy - begy - 1][lox - begx - 1] - original->a[y][x];
                            difb = tmp1->b[loy - begy - 1][lox - begx - 1] - original->b[y][x];
                        }
                    } else if(call == 1) {//dcrop
                        difL = tmp1->L[y][x] - original->L[y][x];
                        difa = tmp1->a[y][x] - original->a[y][x];
                        difb = tmp1->b[y][x] - original->b[y][x];

                    }

                    transformed->L[y][x] = original->L[y][x] + difL;
                    transformed->a[y][x] = original->a[y][x] + difa;
                    transformed->b[y][x] = original->b[y][x] + difb;
                }
            }
        }
    }



}

void ImProcFunctions::BlurNoise_Local(const struct local_params& lp, LabImage* original, LabImage* transformed, const LabImage* const tmp1, int cx, int cy)
{
//local blur and noise
    //BENCHFUN

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
    // BENCHFUN
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



void ImProcFunctions::Reti_Local(int call, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const float lumaref, const struct local_params& lp, float **deltE, LabImage* original, LabImage* transformed, const LabImage* const tmp1, int cx, int cy, int chro)
{

//local retinex
    //BENCHFUN

    {
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

        const float apl = (-1.f) / delhu;
        const float bpl = - apl * hueplus;
        const float amo = 1.f / delhu;
        const float bmo = - amo * huemoins;


        const float pb = 4.f;
        const float pa = (1.f - pb) / 40.f;

        const float ahu = 1.f / (2.8f * lp.sensh - 280.f);
        const float bhu = 1.f - ahu * 2.8f * lp.sensh;

        const float alum = 1.f / (lp.sensh - 100.f);
        const float blum = 1.f - alum * lp.sensh;


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
                    float eps = 0.f;

                    if(fabs(original->b[y][x]) < 0.001f) {
                        eps = 0.01f;
                    }

                    float kab = original->a[y][x] / (original->b[y][x] + eps);

                    float realstr = 1.f;
                    float realstrch = 1.f;
                    //prepare shape detection
                    float deltachro = fabs(rchro - chromaref);
                    float deltahue = fabs(rhue - hueref);

                    if(deltahue > M_PI) {
                        deltahue = -(deltahue - 2.f * M_PI);
                    }

                    float deltaE = 20.f * deltahue + deltachro; //between 0 and 280
                    float deltaL = fabs (lumaref - rL); //between 0 and 100

                    float kch = 1.f;
                    float khu = 0.f;
                    float fach = 1.f;
                    float falu = 1.f;

                    if(deltachro < 160.f * SQR(lp.sensh / 100.f)) {
                        kch = 1.f;
                    } else {
                        float ck = 160.f * SQR(lp.sensh / 100.f);
                        float ak = 1.f / (ck - 160.f);
                        float bk = -160.f * ak;
                        kch = ak * deltachro + bk;
                    }

                    if(lp.sensh < 40.f ) {
                        kch = pow(kch, pa * lp.sensh + pb);    //increase under 40
                    }

                    bool kzon = false;

                    //transition = difficult to avoid artifact with scope on flat area (sky...)
                    //hue detection
                    if((hueref + dhue) < M_PI && rhue < hueplus && rhue > huemoins) {//transition are good
                        if(rhue >= hueplus - delhu)  {
                            realstr = aplus * rhue + bplus;
                            khu  = apl * rhue + bpl;

                        } else if(rhue < huemoins + delhu)  {
                            realstr = amoins * rhue + bmoins;
                            khu = amo * rhue + bmo;

                        } else {
                            realstr = strn;
                            khu = 1.f;

                        }

                        kzon = true;
                    } else if((hueref + dhue) >= M_PI && (rhue > huemoins  || rhue < hueplus )) {
                        if(rhue >= hueplus - delhu  && rhue < hueplus)  {
                            realstr = aplus * rhue + bplus;
                            khu  = apl * rhue + bpl;

                        } else if(rhue >= huemoins && rhue < huemoins + delhu)  {
                            realstr = amoins * rhue + bmoins;
                            khu = amo * rhue + bmo;

                        } else {
                            realstr = strn;
                            khu = 1.f;

                        }

                        kzon = true;
                    }

                    if((hueref - dhue) > -M_PI && rhue < hueplus && rhue > huemoins) {
                        if(rhue >= hueplus - delhu  && rhue < hueplus)  {
                            realstr = aplus * rhue + bplus;
                            khu  = apl * rhue + bpl;

                        } else if(rhue >= huemoins && rhue < huemoins + delhu)  {
                            realstr = amoins * rhue + bmoins;
                            khu = amo * rhue + bmo;

                        } else {
                            realstr = strn;
                            khu = 1.f;

                        }

                        kzon = true;
                    } else if((hueref - dhue) <= -M_PI && (rhue > huemoins  || rhue < hueplus )) {
                        if(rhue >= hueplus - delhu  && rhue < hueplus)  {
                            realstr = aplus * rhue + bplus;
                            khu  = apl * rhue + bpl;

                        } else if(rhue >= huemoins && rhue < huemoins + delhu)  {
                            realstr = amoins * rhue + bmoins;
                            khu = amo * rhue + bmo;

                        } else {
                            realstr = strn;
                            khu = 1.f;

                        }

                        kzon = true;
                    }

                    //shape detection for hue chroma and luma
                    if(lp.sensh <= 20.f) {//to try...

                        if(deltaE <  2.8f * lp.sensh) {
                            fach = khu;
                        } else {
                            fach = khu * (ahu * deltaE + bhu);
                        }

                        float kcr = 10.f;

                        if (rchro < kcr) {
                            fach *= (1.f / (kcr * kcr)) * rchro * rchro;
                        }

                        if(lp.qualmet == 1) {
                            if(deltE[y][x] > lp.thr) {
                                fach = 1.f;
                            }
                        } else {
                            fach = 1.f;
                        }

                        if(deltaL <  lp.sensh) {
                            falu = 1.f;
                        } else {
                            falu = alum * deltaL + blum;
                        }

                    }


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
                    int begx = int(lp.xc - lp.lxL);
                    int begy = int(lp.yc - lp.lyT);

                    if(rL > 0.1f) {//to avoid crash with very low gamut in rare cases ex : L=0.01 a=0.5 b=-0.9
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
                                    float difL;

                                    if(call == 2) {
                                        if(lox >= (lp.xc - lp.lxL) && lox < (lp.xc + lp.lx) && loy >= (lp.yc - lp.lyT) && loy < (lp.yc + lp.ly)) {
                                            difL = tmp1->L[loy - begy - 1][lox - begx - 1] - original->L[y][x];
                                            difL *= factorx * (100.f + realstr * (1.f - factorx)) / 100.f;
                                            difL *= kch * fach;

                                            transformed->L[y][x] = original->L[y][x] + difL;
                                        }
                                    } else {
                                        difL = (tmp1->L[y][x]) - original->L[y][x];
                                        difL *= factorx * (100.f + realstr * (1.f - factorx)) / 100.f;
                                        difL *= kch * fach;
                                        transformed->L[y][x] = original->L[y][x] + difL;
                                    }
                                }

                                if(chro == 1) {
                                    float difa, difb;

                                    if(call == 2) {
                                        if(lox >= (lp.xc - lp.lxL) && lox < (lp.xc + lp.lx) && loy >= (lp.yc - lp.lyT) && loy < (lp.yc + lp.ly)) {
                                            difa = tmp1->a[loy - begy - 1][lox - begx - 1] - original->a[y][x];
                                            difb = tmp1->b[loy - begy - 1][lox - begx - 1] - original->b[y][x];
                                            difa *= factorx * (100.f + realstr * falu * (1.f - factorx)) / 100.f;
                                            difb *= factorx * (100.f + realstr * falu * (1.f - factorx)) / 100.f;
                                            transformed->a[y][x] = CLIPC(original->a[y][x] + difa);
                                            transformed->b[y][x] = CLIPC(original->b[y][x] + difb);
                                        }

                                    } else {
                                        difa = tmp1->a[y][x] - original->a[y][x];
                                        difb = tmp1->b[y][x] - original->b[y][x];
                                        difa *= factorx * (100.f + realstr * falu * (1.f - factorx)) / 100.f;
                                        difb *= factorx * (100.f + realstr * falu * (1.f - factorx)) / 100.f;
                                        transformed->a[y][x] = CLIPC(original->a[y][x] + difa);
                                        transformed->b[y][x] = CLIPC(original->b[y][x] + difb);
                                    }
                                }

                                break;

                            }

                            case 2: { // inside selection => full effect, no transition
                                if(chro == 0) {
                                    float difL;

                                    if(call == 2) {
                                        if(lox >= (lp.xc - lp.lxL) && lox < (lp.xc + lp.lx) && loy >= (lp.yc - lp.lyT) && loy < (lp.yc + lp.ly)) {
                                            difL = tmp1->L[loy - begy - 1][lox - begx - 1] - original->L[y][x];
                                            difL *= (100.f + realstr) / 100.f;
                                            difL *= kch * fach;
                                            transformed->L[y][x] = original->L[y][x] + difL;
                                        }

                                    } else {
                                        difL = tmp1->L[y][x] - original->L[y][x];
                                        difL *= (100.f + realstr) / 100.f;
                                        difL *= kch * fach;
                                        transformed->L[y][x] = original->L[y][x] + difL;
                                    }
                                }

                                if(chro == 1) {
                                    float difa, difb;

                                    if(call == 2) {
                                        if(lox >= (lp.xc - lp.lxL) && lox < (lp.xc + lp.lx) && loy >= (lp.yc - lp.lyT) && loy < (lp.yc + lp.ly)) {
                                            difa = tmp1->a[loy - begy - 1][lox - begx - 1] - original->a[y][x];
                                            difb = tmp1->b[loy - begy - 1][lox - begx - 1] - original->b[y][x];
                                            difa *= (100.f + realstr * falu) / 100.f;
                                            difb *= (100.f + realstr * falu) / 100.f;
                                            transformed->a[y][x] = CLIPC(original->a[y][x] + difa);
                                            transformed->b[y][x] = CLIPC(original->b[y][x] + difb);
                                        }
                                    } else {
                                        difa = tmp1->a[y][x] - original->a[y][x];
                                        difb = tmp1->b[y][x] - original->b[y][x];
                                        difa *= (100.f + realstr * falu) / 100.f;
                                        difb *= (100.f + realstr * falu) / 100.f;
                                        transformed->a[y][x] = CLIPC(original->a[y][x] + difa);
                                        transformed->b[y][x] = CLIPC(original->b[y][x] + difb);
                                    }

                                }
                            }
                        }
                    }
                }
            }
        }

    }
}

void ImProcFunctions::InverseBlurNoise_Local(const struct local_params & lp, LabImage * original, LabImage * transformed, const LabImage * const tmp1, int cx, int cy)
{
    // BENCHFUN
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

void ImProcFunctions::Contrast_Local(float moy, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, float pm, struct local_contra & lco, float lumaref, float av, const struct local_params & lp, float **deltE, LabImage * original, LabImage * transformed, int cx, int cy)
{
    // BENCHFUN
// contrast - perhaps for 4 areas   if need
// I tried shmap adaptaed to Lab, but no real gain and artifacts
    const float localtype = lumaref; // always spot area
    const float ach = (float)lp.trans / 100.f;
    float reducac;

    //constant and variable to prepare shape detection
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


    const float pb = 4.f;
    const float pa = (1.f - pb) / 40.f;

    const float ahu = 1.f / (2.8f * lp.sens - 280.f);
    const float bhu = 1.f - ahu * 2.8f * lp.sens;

    lco.alinf = realcox / (localtype / 2.f);
    const float vi = (localtype / 2.f) / 100.f;
    const float vinf = (50.f + localtype / 2.f) / 100.f;
    ImProcFunctions::secondeg_begin (reducac, vi, lco.aa, lco.bb);//parabolic
    ImProcFunctions::secondeg_end (reducac, vinf, lco.aaa, lco.bbb, lco.ccc);//parabolic

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
                int zone;
                float localFactor = 1.f;
                calcTransition (lox, loy, ach, lp, zone, localFactor);
                //prepare shape detection
                float khu = 0.f;
                float kch = 1.f;
                bool kzon = false;
                float fach = 1.f;
                float deltachro = fabs(rchro - chromaref);
                float deltahue = fabs(rhue - hueref);

                if(deltahue > M_PI) {
                    deltahue = -(deltahue - 2.f * M_PI);
                }

                float deltaE = 20.f * deltahue + deltachro; //pseudo deltaE between 0 and 280
                float rL = original->L[y][x] / 327.68f;

                //kch to modulate action with chroma
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


                // algo with detection of hue ==> artifacts for noisy images  ==> denoise before
                if(lp.sens < 20.f) {//to try...
                    //hue detection
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

                    if(deltaE <  2.8f * lp.sens) {
                        fach = khu;
                    } else {
                        fach = khu * (ahu * deltaE + bhu);
                    }

                    float kcr = 10.f;

                    if (rchro < kcr) {
                        fach *= (1.f / (kcr * kcr)) * rchro * rchro;
                    }

                    if(lp.qualmet == 1) {
                        if(deltE[y][x] > lp.thr) {
                            fach = 1.f;
                        }
                    } else {
                        fach = 1.f;
                    }

                    //fach = khu ;

                } else {
                    /*
                        float kcr = 8.f;
                        if(lp.sensh > 30.f){
                        if (rchro < kcr) {
                            fach *= (1.f / (kcr)) * rchro;

                        }
                        }
                        */
                }

                if(rL > 0.01f) {//to avoid crash with very low gamut in rare cases ex : L=0.01 a=0.5 b=-0.9

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
    }
}


void ImProcFunctions::InverseContrast_Local(float ave, const local_contra & lco, const struct local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy)
{
    //  BENCHFUN
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

static void calclight (float lum, int  koef, float & lumnew)
//replace L-curve that does not work in local
{
    if(koef > 0) {
        lumnew = lum + 0.2f * (33000.f - lum) * (float) koef / 100.f;

        if(lumnew > 32500.f) {
            float kc = 32500.f / lumnew;
            lumnew = lum + kc * 0.2f * (33000.f - lum) * (float)koef / 100.f;

        }
    }

    if(koef < 0) {
        lumnew = lum + 0.2f * lum * (float)koef / 100.f;

        if(lumnew < 0.f) {
            float kc = lum / (lum - lumnew);
            lumnew = lum + kc * 0.2f * lum * (float)koef / 100.f;

        }
    }


}

void ImProcFunctions::InverseSharp_Local(int sp, float **loctemp, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const float lumaref, const local_params & lp, float **deltE, LabImage * original, LabImage * transformed, int cx, int cy)
{
//local blur and noise
    //  BENCHFUN
    const float localtype = lumaref; // always spot area
    const float ach = (float)lp.trans / 100.f;
    float reducac;

    //constant and variable to prepare shape detection
    if(lp.senssha < 30.f) {
        reducac = 0.2f * (lp.senssha / 100.f);
    } else {
        float areduc = 0.6285714f; //0.44f/0.7f;
        float breduc = 0.5f - areduc;
        reducac = areduc * (lp.senssha / 100.f) + breduc;
    }



    constexpr float delhu = 0.1f; //between 0.05 and 0.2

    const float apl = (-1.f) / delhu;
    const float bpl = - apl * hueplus;
    const float amo = 1.f / delhu;
    const float bmo = - amo * huemoins;


    const float pb = 4.f;
    const float pa = (1.f - pb) / 40.f;

    const float ahu = 1.f / (2.8f * lp.senssha - 280.f);
    const float bhu = 1.f - ahu * 2.8f * lp.senssha;

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
                int zone;
                float localFactor = 1.f;
                calcTransition (lox, loy, ach, lp, zone, localFactor);
                //prepare shape detection
                float khu = 0.f;
                float kch = 1.f;
                bool kzon = false;
                float fach = 1.f;
                float deltachro = fabs(rchro - chromaref);
                float deltahue = fabs(rhue - hueref);

                if(deltahue > M_PI) {
                    deltahue = -(deltahue - 2.f * M_PI);
                }

                float deltaE = 20.f * deltahue + deltachro; //pseudo deltaE between 0 and 280

                //kch to modulate action with chroma
                if(deltachro < 160.f * SQR(lp.senssha / 100.f)) {
                    kch = 1.f;
                } else {
                    float ck = 160.f * SQR(lp.senssha / 100.f);
                    float ak = 1.f / (ck - 160.f);
                    float bk = -160.f * ak;
                    kch = ak * deltachro + bk;
                }

                if(lp.senssha < 40.f ) {
                    kch = pow(kch, pa * lp.senssha + pb);    //increase under 40
                }


                // algo with detection of hue ==> artifacts for noisy images  ==> denoise before
                if(lp.senssha < 20.f) {//to try...
                    //hue detection
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

                    if(deltaE <  2.8f * lp.senssha) {
                        fach = khu;
                    } else {
                        fach = khu * (ahu * deltaE + bhu);
                    }


                    float kcr = 10.f;

                    if (rchro < kcr) {
                        fach *= (1.f / (kcr * kcr)) * rchro * rchro;
                    }

                    if(lp.qualmet == 1) {
                        if(deltE[y][x] > lp.thr) {
                            fach = 1.f;
                        }
                    } else {
                        fach = 1.f;
                    }

                    //fach = khu ;

                } else {
                    /*
                        float kcr = 8.f;
                        if(lp.senssha > 30.f){
                        if (rchro < kcr) {
                            fach *= (1.f / (kcr)) * rchro;

                        }
                        }
                        */
                }



                switch(zone) {
                    case 0: { // outside selection and outside transition zone => full effect, no transition
                        float difL = loctemp[y][x] - original->L[y][x];
                        transformed->L[y][x] = original->L[y][x] + difL * kch * fach;

                        break;
                    }

                    case 1: { // inside transition zone
                        float difL = loctemp[y][x] - original->L[y][x];

                        float factorx = 1.f - localFactor;
                        difL *= factorx;

                        transformed->L[y][x] = original->L[y][x] + difL * kch * fach;
                        break;
                    }

                    case 2: { // inside selection => no effect, keep original values
                        transformed->L[y][x] = original->L[y][x];
                    }
                }
            }
        }
    }
}


void ImProcFunctions::Sharp_Local(int call, int sp, float **loctemp, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const float lumaref, const local_params & lp, float **deltE, LabImage * original, LabImage * transformed, int cx, int cy)
{
//local blur and noise
    // BENCHFUN
    const float localtype = lumaref; // always spot area
    const float ach = (float)lp.trans / 100.f;
    float reducac;

    //constant and variable to prepare shape detection
    if(lp.senssha < 30.f) {
        reducac = 0.2f * (lp.senssha / 100.f);
    } else {
        float areduc = 0.6285714f; //0.44f/0.7f;
        float breduc = 0.5f - areduc;
        reducac = areduc * (lp.senssha / 100.f) + breduc;
    }

    //printf("call=%i\n", call);

    constexpr float delhu = 0.1f; //between 0.05 and 0.2

    const float apl = (-1.f) / delhu;
    const float bpl = - apl * hueplus;
    const float amo = 1.f / delhu;
    const float bmo = - amo * huemoins;


    const float pb = 4.f;
    const float pa = (1.f - pb) / 40.f;

    const float ahu = 1.f / (2.8f * lp.senssha - 280.f);
    const float bhu = 1.f - ahu * 2.8f * lp.senssha;

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
                int zone;
                float localFactor = 1.f;
                calcTransition (lox, loy, ach, lp, zone, localFactor);
                //prepare shape detection
                float khu = 0.f;
                float kch = 1.f;
                bool kzon = false;
                float fach = 1.f;
                float deltachro = fabs(rchro - chromaref);
                float deltahue = fabs(rhue - hueref);

                if(deltahue > M_PI) {
                    deltahue = -(deltahue - 2.f * M_PI);
                }

                float deltaE = 20.f * deltahue + deltachro; //pseudo deltaE between 0 and 280

                //kch to modulate action with chroma
                if(deltachro < 160.f * SQR(lp.senssha / 100.f)) {
                    kch = 1.f;
                } else {
                    float ck = 160.f * SQR(lp.senssha / 100.f);
                    float ak = 1.f / (ck - 160.f);
                    float bk = -160.f * ak;
                    kch = ak * deltachro + bk;
                }

                if(lp.senssha < 40.f ) {
                    kch = pow(kch, pa * lp.senssha + pb);    //increase under 40
                }


                // algo with detection of hue ==> artifacts for noisy images  ==> denoise before
                if(lp.senssha < 20.f) {//to try...
                    //hue detection
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

                    if(deltaE <  2.8f * lp.senssha) {
                        fach = khu;
                    } else {
                        fach = khu * (ahu * deltaE + bhu);
                    }


                    float kcr = 10.f;

                    if (rchro < kcr) {
                        fach *= (1.f / (kcr * kcr)) * rchro * rchro;
                    }

                    if(lp.qualmet == 1) {
                        if(deltE[y][x] > lp.thr) {
                            fach = 1.f;
                        }
                    } else {
                        fach = 1.f;
                    }

                    //fach = khu ;

                } else {
                    /*
                        float kcr = 8.f;
                        if(lp.senssha > 30.f){
                        if (rchro < kcr) {
                            fach *= (1.f / (kcr)) * rchro;

                        }
                        }
                        */
                }

                int begx = int(lp.xc - lp.lxL);
                int begy = int(lp.yc - lp.lyT);

                switch(zone) {
                    case 0: { // outside selection and outside transition zone => no effect, keep original values
                        transformed->L[y][x] = original->L[y][x];
                        break;
                    }

                    case 1: { // inside transition zone
                        float factorx = localFactor;
                        float difL;

                        if(call == 2) {
                            if(lox >= (lp.xc - lp.lxL) && lox < (lp.xc + lp.lx) && loy >= (lp.yc - lp.lyT) && loy < (lp.yc + lp.ly)) {
                                difL = loctemp[loy - begy - 1][lox - begx - 1] - original->L[y][x];
                            }
                        } else if(call == 1) {
                            difL = loctemp[y][x] - original->L[y][x];

                        }

                        //float difL = loctemp[y][x] - original->L[y][x];
                        difL *= factorx;
                        transformed->L[y][x] = original->L[y][x] + difL * kch * fach;

                        break;
                    }

                    case 2: { // inside selection => full effect, no transition
                        // float difL = loctemp[y][x] - original->L[y][x];
                        float difL;

                        if(call == 2) {

                            if(lox >= (lp.xc - lp.lxL) && lox < (lp.xc + lp.lx) && loy >= (lp.yc - lp.lyT) && loy < (lp.yc + lp.ly)) {
                                //       bufsh[loy - begy - 1][lox - begx - 1]
                                difL = loctemp[loy - begy - 1][lox - begx - 1] - original->L[y][x];
                            }
                        } else if(call == 1) {
                            difL = loctemp[y][x] - original->L[y][x];
                        }

                        transformed->L[y][x] = original->L[y][x] + difL * kch * fach;
                    }
                }
            }
        }
    }
}



void ImProcFunctions::ColorLight_Local(int sp, float moy, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const float lumaref, const local_params & lp, float ** deltE, LabImage * original, LabImage * transformed, int cx, int cy)
{
    // BENCHFUN
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


    const float apl = (-1.f) / delhu;
    const float bpl = - apl * hueplus;
    const float amo = 1.f / delhu;
    const float bmo = - amo * huemoins;


    const float pb = 4.f;
    const float pa = (1.f - pb) / 40.f;

    const float ahu = 1.f / (2.8f * lp.sens - 280.f);
    const float bhu = 1.f - ahu * 2.8f * lp.sens;

    const float alum = 1.f / (lp.sens - 100.f);
    const float blum = 1.f - alum * lp.sens;

    //luma
    constexpr float lumdelta = 11.f; //11
    float modlum = lumdelta * multlum;
//  printf("multlum=%f modlum=%f\n", multlum, modlum);

    // constant and varaibles to prepare shape detection
    if(lumaref + modlum >= 100.f) {
        modlum = (100.f - lumaref) / 2.f;
    }

    if(lumaref - modlum <= 0.f) {
        modlum = (lumaref) / 2.f;
    }

    float alu = 1.f / (lumaref + modlum - 100.f); //linear
    float aa, bb, aaa, bbb, ccc;
    float reducac = settings->reduchigh;//0.85f;
    float reducac2 = settings->reduclow;//0.2f;

    float vinf = (lumaref + modlum) / 100.f;
    float vi = (lumaref - modlum) / 100.f;
    ImProcFunctions::secondeg_begin (reducac, vi, aa, bb);//parabolic
    ImProcFunctions::secondeg_end (reducac, vinf, aaa, bbb, ccc);//parabolic
//  printf("vi=%f aa=%f bb=%f vinf=%f aaa=%f bbb=%f ccc=%f\n", vi,aa,bb, vinf, aaa, bbb, ccc);
    float vinf2 = (lumaref + modlum) / 100.f;
    float vi2 = (lumaref - modlum) / 100.f;
    float aaaa, bbbb, cccc, aO, bO;
    ImProcFunctions::secondeg_end (reducac2, vinf2, aaaa, bbbb, cccc);//parabolic
    ImProcFunctions::secondeg_begin (reducac2, vi2, aO, bO);//parabolic

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef __SSE2__
        float atan2Buffer[transformed->W] ALIGNED16;
        float sqrtBuffer[transformed->W] ALIGNED16;
        vfloat c327d68v = F2V(327.68f);
#endif

        float maxl = -100000.f;
        float maxa = -100000.f;
        float maxb = -100000.f;
        float minl = 100000.f;
        float mina = 100000.f;
        float minb = 100000.f;
        float maxrl = -100000.f;
        float minrl = 100000.f;

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
                float rLL = original->L[y][x] / 327.68f;

                if(fabs(original->b[y][x]) < 0.01f) {
                    original->b[y][x] = 0.01f;
                }

                float eps = 0.f;

                if(fabs(original->b[y][x]) < 0.001f) {
                    eps = 0.01f;
                }

                float kab = (original->a[y][x] / (original->b[y][x] + eps));
                //prepare shape detection
                float realchro = 1.f;
                float deltachro = fabs(rchro - chromaref);

                float deltahue = fabs(rhue - hueref);

                if(deltahue > M_PI) {
                    deltahue = -(deltahue - 2.f * M_PI);
                }

                float deltaE = 20.f * deltahue + deltachro; //pseudo deltaE between 0 and 280
                float deltaL = fabs (lumaref - rL); //between 0 and 100

                float kch = 1.f;
                float khu = 0.f;
                float fach = 1.f;
                float falu = 1.f;

                //kch acts on luma
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
                //hue detection
                if((hueref + dhue) < M_PI && rhue < hueplus && rhue > huemoins) {//transition are good
                    if(rhue >= hueplus - delhu)  {
                        realchro = aplus * rhue + bplus;
                        khu  = apl * rhue + bpl;

                    } else if(rhue < huemoins + delhu)  {
                        realchro = amoins * rhue + bmoins;
                        khu = amo * rhue + bmo;

                    } else {
                        realchro = lp.chro;
                        khu = 1.f;

                    }

                    kzon = true;
                } else if((hueref + dhue) >= M_PI && (rhue > huemoins  || rhue < hueplus )) {
                    if(rhue >= hueplus - delhu  && rhue < hueplus)  {
                        realchro = aplus * rhue + bplus;
                        khu  = apl * rhue + bpl;

                    } else if(rhue >= huemoins && rhue < huemoins + delhu)  {
                        realchro = amoins * rhue + bmoins;
                        khu = amo * rhue + bmo;

                    } else {
                        realchro = lp.chro;
                        khu = 1.f;

                    }

                    kzon = true;
                }

                if((hueref - dhue) > -M_PI && rhue < hueplus && rhue > huemoins) {
                    if(rhue >= hueplus - delhu  && rhue < hueplus)  {
                        realchro = aplus * rhue + bplus;
                        khu  = apl * rhue + bpl;

                    } else if(rhue >= huemoins && rhue < huemoins + delhu)  {
                        realchro = amoins * rhue + bmoins;
                        khu = amo * rhue + bmo;

                    } else {
                        realchro = lp.chro;
                        khu = 1.f;

                    }

                    kzon = true;
                } else if((hueref - dhue) <= -M_PI && (rhue > huemoins  || rhue < hueplus )) {
                    if(rhue >= hueplus - delhu  && rhue < hueplus)  {
                        realchro = aplus * rhue + bplus;
                        khu  = apl * rhue + bpl;

                    } else if(rhue >= huemoins && rhue < huemoins + delhu)  {
                        realchro = amoins * rhue + bmoins;
                        khu = amo * rhue + bmo;

                    } else {
                        realchro = lp.chro;
                        khu = 1.f;

                    }

                    kzon = true;
                }


                //detection of deltaE and deltaL
                if(lp.sens <= 20.f) {//to try...
                    //fach and kch acts on luma
                    if(deltaE <  2.8f * lp.sens) {
                        fach = khu;
                    } else {
                        fach = khu * (ahu * deltaE + bhu);
                    }

                    float kcr = 10.f;

                    if (rchro < kcr) {
                        fach *= (1.f / (kcr * kcr)) * rchro * rchro;
                    }

                    //fach = 1.f;//to avoid artifacts in some cases
                    //can be probably improved
                    if(lp.qualmet == 1) {
                        if(deltE[y][x] > lp.thr) {
                            fach = 1.f;
                        }
                    } else {
                        fach = 1.f;
                    }

                    //falu acts on chroma
                    if(deltaL <  lp.sens) {
                        falu = 1.f;
                    } else {
                        falu = 1.f;// alum * deltaL + blum;
                    }

                }

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

                float kLinf = rLL / (100.f);
                float kLsup = kLinf;

                float kdiff = 1.f;

                if(kzon) {///rhue < hueplus && rhue > huemoins

                    if( (rLL > (lumaref - modlum) && rLL < (lumaref + modlum))) {
                        kdiff = 1.f;
                    } else if (rLL > 0.f && rLL <= (lumaref - modlum)) {
                        kdiff = (aa * kLinf * kLinf + bb * kLinf);   //parabolic

                        if(kdiff < 0.01f) {
                            kdiff = 0.01f;
                        }
                    } else if (rLL <= 100.f && rLL >= (lumaref + modlum)) {

                        kdiff = (aaa * kLsup * kLsup + bbb * kLsup + ccc);   //parabolic

                        if(kdiff < 0.01f) {
                            kdiff = 0.01f;
                        }

                    }

                    //end luma
                } else {
                    float ktes = 1.f;

                    if( (rLL > (lumaref - modlum) && rLL < (lumaref + modlum))) {
                        kdiff = ktes;
                    } else if (rLL > 0.f && rLL <= (lumaref - modlum)) {

                        kdiff = (ktes * (aO * kLinf * kLinf + bO * kLinf));    //parabolic

                        if(kdiff < 0.01f) {
                            kdiff = 0.01f;
                        }

                    } else if (rLL <= 100.f && rLL >= (lumaref + modlum)) {

                        kdiff = (ktes * (aaaa * kLsup * kLsup + bbbb * kLsup + cccc));    //parabolic

                        if(kdiff < 0.01f) {
                            kdiff = 0.01f;
                        }

                    }

                }

                int zone;
                float localFactor;
                calcTransition (lox, loy, ach, lp, zone, localFactor);
                float th_r = 0.01f;

                /*
                                if(moy < 50.f) {
                                    th_r = 0.3f;
                                } else {
                                    th_r = 20.f;
                                }
                */
                if(rL > th_r) {//to avoid crash with very low gamut in rare cases ex : L=0.01 a=0.5 b=-0.9
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
                                calclight (original->L[y][x], lp.ligh , lumnew);//replace L-curve
                            }

                            float lightcont = lumnew ; //original->L[y][x] + (lp.ligh /100.f)*original->L[y][x] ; //apply lightness
                            float factorx = localFactor;
                            float fac = (100.f + factorx * realchro * falu) / 100.f; //chroma factor transition
                            float diflc = lightcont - original->L[y][x];
                            kdiff *= fach * kch;
                            diflc *= kdiff ;

                            diflc *= factorx; //transition lightess
                            transformed->L[y][x] = CLIPL(original->L[y][x] + diflc);

                            if(fabs(kab) > 1.f) {
                                transformed->a[y][x] = CLIPC(original->a[y][x] * fac) ;
                                transformed->b[y][x] = CLIPC(original->a[y][x] * fac) / kab;
                            } else {
                                transformed->b[y][x] = CLIPC(original->b[y][x] * fac);
                                transformed->a[y][x] = CLIPC(original->b[y][x] * fac) * kab ;

                            }

                            break;
                        }

                        case 2: { // inside selection => full effect, no transition
                            //  float lightcont = original->L[y][x] + (lp.ligh /100.f)*original->L[y][x]; //apply lightness
                            float lumnew = original->L[y][x];

                            if(lp.ligh != 0) {
                                calclight (original->L[y][x], lp.ligh , lumnew);

                            }


                            //    float lightcont = localcurve[original->L[y][x]]; //apply lightness
                            float lightcont = lumnew ; //original->L[y][x] + (lp.ligh /100.f)*original->L[y][x] ; //apply lightness

                            float fac = (100.f + realchro * falu) / 100.f; //chroma factor transition
                            float diflc = lightcont - original->L[y][x];
                            kdiff *= fach * kch;
                            diflc *= kdiff ;

                            transformed->L[y][x] = CLIPL(original->L[y][x] + diflc);

                            if(fabs(kab) > 1.f) {
                                transformed->a[y][x] = CLIPC(original->a[y][x] * fac) ;
                                transformed->b[y][x] = CLIPC(original->a[y][x] * fac) / kab;
                            } else {
                                transformed->b[y][x] = CLIPC(original->b[y][x] * fac);
                                transformed->a[y][x] = CLIPC(original->b[y][x] * fac) * kab;

                            }


                        }

                            //      if(rL > maxrl) maxrl = rL;
                            //      if(rL < minrl) minrl = rL;


                            /*
                                                if(transformed->L[y][x] > maxl) {
                                                    maxl = transformed->L[y][x];
                                                }

                                                if(transformed->L[y][x] < minl) {
                                                    minl = transformed->L[y][x];
                                                }

                                                if(transformed->a[y][x] > maxa) {
                                                    maxa = transformed->a[y][x];
                                                }

                                                if(transformed->a[y][x] < mina) {
                                                    mina = transformed->a[y][x];
                                                }

                                                if(transformed->b[y][x] > maxb) {
                                                    maxb = transformed->b[y][x];
                                                }

                                                if(transformed->b[y][x] < minb) {
                                                    minb = transformed->b[y][x];
                                                }
                            */
                    }
                }
            }
        }

        //  printf("maxRL=%f minRL=%f\n", maxrl, minrl);
        //     printf("sp=%i ML=%f ml=%f Ma=%f ma=%f Mb=%f mb=%f\n", sp, maxl / 327.f, minl / 327.f, maxa / 327.f, mina / 327.f, maxb / 327.f, minb / 327.f);
    }
}

void ImProcFunctions::InverseColorLight_Local(const struct local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy)
{
    // BENCHFUN
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


void ImProcFunctions::Lab_Local(int call, int sp, float** shbuffer, LabImage * original, LabImage * transformed, int sx, int sy, int cx, int cy, int oW, int oH,  int fw, int fh, bool locutili, int sk, const LocretigainCurve & locRETgainCcurve, double & hueref, double & chromaref, double & lumaref)
{
    //general call of others functions : important return hueref, chromaref, lumaref
    if(params->locallab.enabled) {
        // BENCHFUN
#ifdef _DEBUG
        MyTime t1e, t2e;
        t1e.set();
        // init variables to display Munsell corrections
        MunsellDebugInfo* MunsDebugInfo = new MunsellDebugInfo();
#endif
        //ImProcFunctions ipf;



        int GW = transformed->W;
        int GH = transformed->H;
        float moy = 0.f;
        float maxmad = -10000.f;
        float minmad = 1000000.f;

        struct local_params lp;
        calcLocalParams(oW, oH, params->locallab, lp);

        const float radius = lp.rad / (sk * 2.f); //0 to 50 ==> see skip
        GW = transformed->W;
        GH = transformed->H;
        float **deltE;

        if(lp.qualmet == 1) {

            deltE   = new float*[GH];

            for (int i = 0; i < GH; i++) {
                deltE[i] = new float[GW];
            }

            for(int ir = 0; ir < GH; ir++)
                for(int jr = 0; jr < GW; jr++) {
                    deltE[ir][jr] = 0.f;
                }
        }

        if(radius >= GAUSS_SKIP || lp.stren > 0.1) { // radius < GAUSS_SKIP means no gauss, just copy of original image
            LabImage *tmp1 = new LabImage(transformed->W, transformed->H);;

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
            int spotSize = 0.88623f * max(1,  lp.cir / sk); //18
            //O.88623 = sqrt(PI / 4) ==> sqare equal to circle

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

        //local denoise

        if(lp.noiself > 0.f || lp.noiselc > 0.f || lp.noisecf > 0.f || lp.noisecc > 0.f  && call < 3) {
            if(call == 1) {
                LabImage *tmp1 = new LabImage(transformed->W, transformed->H);
                int GW = transformed->W;
                int GH = transformed->H;


                for(int ir = 0; ir < GH; ir++)
                    for(int jr = 0; jr < GW; jr++) {
                        tmp1->L[ir][jr] = original->L[ir][jr];
                        tmp1->a[ir][jr] = original->a[ir][jr];
                        tmp1->b[ir][jr] = original->b[ir][jr];
                    }

                int DaubLen = 6;
                int wavNestedLevels = 1;

                int levwavL = 7;
                int skip = 1;
                wavelet_decomposition* Ldecomp = new wavelet_decomposition (tmp1->L[0], tmp1->W, tmp1->H, levwavL, 1, skip, max(1, wavNestedLevels), DaubLen);
                wavelet_decomposition* adecomp = new wavelet_decomposition (tmp1->a[0], tmp1->W, tmp1->H, levwavL, 1, skip, max(1, wavNestedLevels), DaubLen);
                wavelet_decomposition* bdecomp = new wavelet_decomposition (tmp1->b[0], tmp1->W, tmp1->H, levwavL, 1, skip, max(1, wavNestedLevels), DaubLen);

                float madL[8][3];
                float madab[8][3];

                if(!Ldecomp->memoryAllocationFailed) {

                    for (int lvl = 0; lvl < 7; lvl++) {
                        for (int dir = 1; dir < 4; dir++) {
                            int Wlvl_L = Ldecomp->level_W(lvl);
                            int Hlvl_L = Ldecomp->level_H(lvl);

                            float ** WavCoeffs_L = Ldecomp->level_coeffs(lvl);

                            madL[lvl][dir - 1] = SQR(Mad(WavCoeffs_L[dir], Wlvl_L * Hlvl_L));
                        }
                    }


                    int ind = 0;

                    float vari[7];

                    vari[0] = 8.f * SQR((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                    vari[1] = 8.f * SQR((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                    vari[2] = 8.f * SQR((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));

                    vari[3] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                    vari[4] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                    vari[5] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                    vari[6] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));

                    if(( lp.noiself > 0.1f ||  lp.noiselc > 0.1f)) {
                        int edge = 1;
                        vari[0] = max(0.0001f, vari[0]);
                        vari[1] = max(0.0001f, vari[1]);
                        vari[2] = max(0.0001f, vari[2]);
                        vari[3] = max(0.0001f, vari[3]);
                        vari[4] = max(0.0001f, vari[4]);
                        vari[5] = max(0.0001f, vari[5]);
                        vari[6] = max(0.0001f, vari[6]);

                        float* noisevarlum = nullptr;  // we need a dummy to pass it to WaveletDenoiseAllL

                        WaveletDenoiseAllL(*Ldecomp, noisevarlum, madL, vari, 2);
                    }
                }

                if(!adecomp->memoryAllocationFailed && !bdecomp->memoryAllocationFailed) {

                    float variC[7];

                    variC[0] = lp.noisecf / 10.0;
                    variC[1] = lp.noisecf / 10.0;
                    variC[2] = lp.noisecf / 10.0;

                    variC[3] = lp.noisecf / 10.0;
                    variC[4] = lp.noisecf / 10.0;
                    variC[5] = lp.noisecc / 10.0;
                    variC[6] = lp.noisecc / 10.0;

                    if(( lp.noisecf > 0.1f ||  lp.noisecc > 0.1f)) {
                        variC[0] = max(0.0001f, variC[0]);
                        variC[1] = max(0.0001f, variC[1]);
                        variC[2] = max(0.0001f, variC[2]);
                        variC[3] = max(0.0001f, variC[3]);
                        variC[4] = max(0.0001f, variC[4]);
                        variC[5] = max(0.0001f, variC[5]);
                        variC[6] = max(0.0001f, variC[6]);

                        float* noisevarchrom = new float[GH * GW];

                        for(int q = 0; q < GH * GW; q++) {
                            noisevarchrom[q] = 1.f;
                        }

                        float noisevarab_r = 100.f; //SQR(lp.noisecc / 10.0);
                        WaveletDenoiseAllAB(*Ldecomp, *adecomp, noisevarchrom, madL, variC, 2, noisevarab_r, false, false, false);
                        WaveletDenoiseAllAB(*Ldecomp, *bdecomp, noisevarchrom, madL, variC, 2, noisevarab_r, false, false, false);
                        delete[] noisevarchrom;

                    }
                }

                if(!Ldecomp->memoryAllocationFailed) {

                    Ldecomp->reconstruct(tmp1->L[0]);
                }

                if(!adecomp->memoryAllocationFailed) {

                    adecomp->reconstruct(tmp1->a[0]);
                }

                if(!bdecomp->memoryAllocationFailed) {

                    bdecomp->reconstruct(tmp1->b[0]);
                }

                DeNoise_Local(call, lp, original, transformed, tmp1, cx, cy);
                delete tmp1;
                delete Ldecomp;
                delete adecomp;
                delete bdecomp;
            }

            if(call == 2) { //simpleprocess
                LabImage *bufwv;
                int GW = transformed->W;
                int GH = transformed->H;
                int bfh = int(lp.ly + lp.lyT) + 1;//bfw bfh real size of square zone
                int bfw = int(lp.lx + lp.lxL) + 1;
                bufwv = new LabImage(bfw, bfh);
#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for(int ir = 0; ir < bfh; ir++) //fill with 0
                    for(int jr = 0; jr < bfw; jr++) {
                        bufwv->L[ir][jr] = 0.f;
                        bufwv->a[ir][jr] = 0.f;
                        bufwv->b[ir][jr] = 0.f;

                    }

#ifdef _OPENMP
//           #pragma omp parallel for
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;
                        int begx = int(lp.xc - lp.lxL);
                        int begy = int(lp.yc - lp.lyT);

                        if(lox >= (lp.xc - lp.lxL) && lox < (lp.xc + lp.lx) && loy >= (lp.yc - lp.lyT) && loy < (lp.yc + lp.ly)) {
                            bufwv->L[loy - begy - 1][lox - begx - 1] = original->L[y][x];//fill square buffer with datas
                            bufwv->a[loy - begy - 1][lox - begx - 1] = original->a[y][x];//fill square buffer with datas
                            bufwv->b[loy - begy - 1][lox - begx - 1] = original->b[y][x];//fill square buffer with datas
                        }
                    }

                int DaubLen = 6;
                int wavNestedLevels = 1;

                int levwavL = 7;
                int skip = 1;
                wavelet_decomposition* Ldecomp = new wavelet_decomposition (bufwv->L[0], bufwv->W, bufwv->H, levwavL, 1, skip, max(1, wavNestedLevels), DaubLen);
                wavelet_decomposition* adecomp = new wavelet_decomposition (bufwv->a[0], bufwv->W, bufwv->H, levwavL, 1, skip, max(1, wavNestedLevels), DaubLen);
                wavelet_decomposition* bdecomp = new wavelet_decomposition (bufwv->b[0], bufwv->W, bufwv->H, levwavL, 1, skip, max(1, wavNestedLevels), DaubLen);

                float madL[8][3];
                float madab[8][3];

//           if(!Ldecomp->memoryAllocationFailed) {
                if(!Ldecomp->memoryAllocationFailed) {

                    for (int lvl = 0; lvl < 7; lvl++) {
                        for (int dir = 1; dir < 4; dir++) {
                            int Wlvl_L = Ldecomp->level_W(lvl);
                            int Hlvl_L = Ldecomp->level_H(lvl);

                            float ** WavCoeffs_L = Ldecomp->level_coeffs(lvl);

                            madL[lvl][dir - 1] = SQR(Mad(WavCoeffs_L[dir], Wlvl_L * Hlvl_L));
                        }
                    }


                    int ind = 0;

                    float vari[7];
                    vari[0] = 8.f * SQR((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                    vari[1] = 8.f * SQR((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                    vari[2] = 8.f * SQR((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));

                    vari[3] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                    vari[4] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                    vari[5] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                    vari[6] = 8.f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));

                    if(( lp.noiself > 0.1f ||  lp.noiselc > 0.1f)) {
                        int edge = 1;
                        vari[0] = max(0.0001f, vari[0]);
                        vari[1] = max(0.0001f, vari[1]);
                        vari[2] = max(0.0001f, vari[2]);
                        vari[3] = max(0.0001f, vari[3]);
                        vari[4] = max(0.0001f, vari[4]);
                        vari[5] = max(0.0001f, vari[5]);
                        vari[6] = max(0.0001f, vari[6]);

                        float* noisevarlum = nullptr;  // we need a dummy to pass it to WaveletDenoiseAllL

                        WaveletDenoiseAllL(*Ldecomp, noisevarlum, madL, vari, 2);
                    }
                }

                if(!adecomp->memoryAllocationFailed && !bdecomp->memoryAllocationFailed) {

                    float variC[7];

                    variC[0] = lp.noisecf / 10.0;
                    variC[1] = lp.noisecf / 10.0;
                    variC[2] = lp.noisecf / 10.0;

                    variC[3] = lp.noisecf / 10.0;
                    variC[4] = lp.noisecf / 10.0;
                    variC[5] = lp.noisecc / 10.0;
                    variC[6] = lp.noisecc / 10.0;

                    if(( lp.noisecf > 0.1f ||  lp.noisecc > 0.1f)) {
                        variC[0] = max(0.0001f, variC[0]);
                        variC[1] = max(0.0001f, variC[1]);
                        variC[2] = max(0.0001f, variC[2]);
                        variC[3] = max(0.0001f, variC[3]);
                        variC[4] = max(0.0001f, variC[4]);
                        variC[5] = max(0.0001f, variC[5]);
                        variC[6] = max(0.0001f, variC[6]);

                        float* noisevarchrom = new float[bfh * bfw];

                        for(int q = 0; q < bfh * bfw; q++) {
                            noisevarchrom[q] = 1.f;
                        }

                        float noisevarab_r = 100.f; //SQR(lp.noisecc / 10.0);
                        WaveletDenoiseAllAB(*Ldecomp, *adecomp, noisevarchrom, madL, variC, 2, noisevarab_r, false, false, false);
                        WaveletDenoiseAllAB(*Ldecomp, *bdecomp, noisevarchrom, madL, variC, 2, noisevarab_r, false, false, false);
                        delete[] noisevarchrom;

                    }
                }

                if(!Ldecomp->memoryAllocationFailed) {

                    Ldecomp->reconstruct(bufwv->L[0]);
                }

                if(!adecomp->memoryAllocationFailed) {

                    adecomp->reconstruct(bufwv->a[0]);
                }

                if(!bdecomp->memoryAllocationFailed) {

                    bdecomp->reconstruct(bufwv->b[0]);
                }

                DeNoise_Local(call, lp, original, transformed, bufwv, cx, cy);
                delete bufwv;
                delete Ldecomp;
                delete adecomp;
                delete bdecomp;


            }

        }




        if(!lp.inv  && (lp.chro != 0 || lp.ligh != 0)) { //interior ellipse renforced lightness and chroma
            float maxhur = -10.f;
            float minhur = 10.f;
            float hueplus = hueref + dhue;
            float huemoins = hueref - dhue;

            if(hueplus > M_PI) {
                hueplus = hueref + dhue - 2.f * M_PI;
            }

            if(huemoins < -M_PI) {
                huemoins = hueref - dhue + 2.f * M_PI;
            }


            ColorLight_Local(sp, moy, hueplus, huemoins, hueref, dhue, chromaref, lumaref, lp, deltE, original, transformed, cx, cy);

        }
        //inverse
        else if(lp.inv  && (lp.chro != 0 || lp.ligh != 0)) {

            InverseColorLight_Local(lp, original, transformed, cx, cy);
        }


        if(!lp.inv  && lp.cont != 0) {   //contrast interior ellipse
            const float pm = lp.cont < 0.f ? -1.f : 1.f;
            float hueplus = hueref + dhue;
            float huemoins = hueref - dhue;

            if(hueplus > M_PI) {
                hueplus = hueref + dhue - 2.f * M_PI;
            }

            if(huemoins < -M_PI) {
                huemoins = hueref - dhue + 2.f * M_PI;
            }

            Contrast_Local(moy, hueplus, huemoins, hueref, dhue, chromaref, pm, lco, lumaref, av, lp, deltE, original, transformed, cx, cy);
        } else if(lp.inv && lp.cont != 0) {

            float multL = (float)lp.cont * (maxl - 1.f) / 100.f + 1.f;
            float multH = (float) lp.cont * (maxh - 1.f) / 100.f + 1.f;

            lco.ah = (multH - 1.f) / (av - 100.f); //av ==> lumaref
            lco.bh = 1.f - 100.f * lco.ah;
            lco.al = (multL - 1.f) / av;
            lco.bl = 1.f;

            InverseContrast_Local(ave, lco, lp, original, transformed, cx, cy);
        }

// end contrast interior and exterior


        if(!lp.invshar && lp.shrad > 0.42 && call < 3) { //interior ellipse for sharpening, call = 1 and 2 only with Dcrop and simpleprocess

            int GW = original->W;
            int GH = original->H;
            float **bufsh;//buffer por square zone
            float **loctemp;
            float **hbuffer;
            int bfh = int(lp.ly + lp.lyT) + 1;//bfw bfh real size of square zone
            int bfw = int(lp.lx + lp.lxL) + 1;

            if(call == 2) { //call from simpleprocess
                //printf("GW=%i GH=%i\n", GW, GH);
                //    int bfh = int(lp.ly + lp.lyT) + 1;//bfw bfh real size of square zone
                //    int bfw = int(lp.lx + lp.lxL) + 1;
                //printf("bw=%i bh=%i\n", bfw, bfh);
                bufsh   = new float*[bfh];

                for (int i = 0; i < bfh; i++) {
                    bufsh[i] = new float[bfw];
                }

#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for(int ir = 0; ir < bfh; ir++) //fill with 0
                    for(int jr = 0; jr < bfw; jr++) {
                        bufsh[ir][jr] = 0.f;
                    }

                //    int begx = int(lp.xc - lp.lxL);
                //    int begy = int(lp.yc - lp.lyT);


#ifdef _OPENMP
//           #pragma omp parallel for
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;
                        int begx = int(lp.xc - lp.lxL);
                        int begy = int(lp.yc - lp.lyT);

                        if(lox >= (lp.xc - lp.lxL) && lox < (lp.xc + lp.lx) && loy >= (lp.yc - lp.lyT) && loy < (lp.yc + lp.ly)) {
                            bufsh[loy - begy - 1][lox - begx - 1] = original->L[y][x];//fill square buffer with datas
                        }
                    }

                loctemp = new float*[bfh];//allocate temp

                for (int i = 0; i < bfh; i++) {
                    loctemp[i] = new float[bfw];
                }

                hbuffer = new float*[bfh];//allocate buffer for sharp

                for (int i = 0; i < bfh; i++) {
                    hbuffer[i] = new float[bfw];
                }



                //      ImProcFunctions::deconvsharpeningloc(original->L, shbuffer, GW, GH, loctemp, params->locallab.shardamping, (double)params->locallab.sharradius / 100., params->locallab.shariter, params->locallab.sharamount);
                //sharpen only square area instaed of all image
                ImProcFunctions::deconvsharpeningloc(bufsh, hbuffer, bfw, bfh, loctemp, params->locallab.shardamping, (double)params->locallab.sharradius / 100., params->locallab.shariter, params->locallab.sharamount);
            } else { //call from dcrop.cc
                loctemp = new float*[GH];//allocate temp

                for (int i = 0; i < GH; i++) {
                    loctemp[i] = new float[GW];
                }

                ImProcFunctions::deconvsharpeningloc(original->L, shbuffer, GW, GH, loctemp, params->locallab.shardamping, (double)params->locallab.sharradius / 100., params->locallab.shariter, params->locallab.sharamount);

            }

            float hueplus = hueref + dhue;
            float huemoins = hueref - dhue;

            if(hueplus > M_PI) {
                hueplus = hueref + dhue - 2.f * M_PI;
            }

            if(huemoins < -M_PI) {
                huemoins = hueref - dhue + 2.f * M_PI;
            }

            //sharpen ellipse and transition
            Sharp_Local(call, sp, loctemp, hueplus, huemoins, hueref, dhue, chromaref, lumaref, lp, deltE, original, transformed, cx, cy);

            //claen all
            if(call == 2) {
                for (int i = 0; i < bfh; i++) {
                    delete [] loctemp[i];
                }

                delete [] loctemp;

                for (int i = 0; i < bfh; i++) {
                    delete [] bufsh[i];
                }

                delete [] bufsh;

                for (int i = 0; i < bfh; i++) {
                    delete [] hbuffer[i];
                }

                delete [] hbuffer;
            } else {
                for (int i = 0; i < GH; i++) {
                    delete [] loctemp[i];
                }

                delete [] loctemp;

            }

            /*            for (int i = 0; i < GH; i++) {
                            delete [] hbuffer[i];
                        }

                        delete [] hbuffer;
            */

        } else if(lp.invshar && lp.shrad > 0.42 && call < 3) {
            int GW = original->W;
            int GH = original->H;

            float **loctemp = new float*[GH];

            for (int i = 0; i < GH; i++) {
                loctemp[i] = new float[GW];
            }

            ImProcFunctions::deconvsharpeningloc(original->L, shbuffer, GW, GH, loctemp, params->locallab.shardamping, (double)params->locallab.sharradius / 100., params->locallab.shariter, params->locallab.sharamount);

            float hueplus = hueref + dhue;
            float huemoins = hueref - dhue;

            if(hueplus > M_PI) {
                hueplus = hueref + dhue - 2.f * M_PI;
            }

            if(huemoins < -M_PI) {
                huemoins = hueref - dhue + 2.f * M_PI;
            }

            InverseSharp_Local(sp, loctemp, hueplus, huemoins, hueref, dhue, chromaref, lumaref, lp, deltE, original, transformed, cx, cy);
            // InverseSharp_Local(sp, loctemp, dataspot, lp, original, transformed, cx, cy);

            for (int i = 0; i < GH; i++) {
                delete [] loctemp[i];
            }

            delete [] loctemp;

        }


        if(lp.str > 0.f) {
            int GW = transformed->W;
            int GH = transformed->H;

            // float **bufreti;//buffer por square zone
            LabImage *bufreti;

            float **loctemp;
            float **hbuffer;
            int bfh = int(lp.ly + lp.lyT) + 1;//bfw bfh real size of square zone
            int bfw = int(lp.lx + lp.lxL) + 1;

            float hueplus = hueref + dhueret;
            float huemoins = hueref - dhueret;

            if(hueplus > M_PI) {
                hueplus = hueref + dhueret - 2.f * M_PI;
            }

            if(huemoins < -M_PI) {
                huemoins = hueref - dhueret + 2.f * M_PI;
            }

            int Hd, Wd;
            Hd = GH;
            Wd = GW;

            if(!lp.invret && call == 2) {
                Hd = bfh;
                Wd = bfw;
                bufreti = new LabImage(bfw, bfh);

                /*                bufreti   = new float*[bfh];

                                for (int i = 0; i < bfh; i++) {
                                    bufreti[i] = new float[bfw];
                                }
                */
#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for(int ir = 0; ir < bfh; ir++) //fill with 0
                    for(int jr = 0; jr < bfw; jr++) {
                        bufreti->L[ir][jr] = 0.f;
                        bufreti->a[ir][jr] = 0.f;
                        bufreti->b[ir][jr] = 0.f;
                    }

                //    int begx = int(lp.xc - lp.lxL);
                //    int begy = int(lp.yc - lp.lyT);


#ifdef _OPENMP
//           #pragma omp parallel for
#endif

                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;
                        int begx = int(lp.xc - lp.lxL);
                        int begy = int(lp.yc - lp.lyT);

                        if(lox >= (lp.xc - lp.lxL) && lox < (lp.xc + lp.lx) && loy >= (lp.yc - lp.lyT) && loy < (lp.yc + lp.ly)) {
                            bufreti->L[loy - begy - 1][lox - begx - 1] = original->L[y][x];//fill square buffer with datas
                            bufreti->a[loy - begy - 1][lox - begx - 1] = original->a[y][x];//fill square buffer with datas
                            bufreti->b[loy - begy - 1][lox - begx - 1] = original->b[y][x];//fill square buffer with datas
                        }
                    }



            }

            float *orig[Hd] ALIGNED16;
            float *origBuffer = new float[Hd * Wd];

            for (int i = 0; i < Hd; i++) {
                orig[i] = &origBuffer[i * Wd];
            }

            float *orig1[Hd] ALIGNED16;
            float *origBuffer1 = new float[Hd * Wd];

            for (int i = 0; i < Hd; i++) {
                orig1[i] = &origBuffer1[i * Wd];
            }


            LabImage *tmpl;

            //    LabImage *tmpl = new LabImage(Wd, Hd);
            if(!lp.invret && call == 2) {

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for(int ir = 0; ir < Hd; ir += 1)
                    for(int jr = 0; jr < Wd; jr += 1) {
                        orig[ir][jr] = bufreti->L[ir][jr];
                        orig1[ir][jr] = bufreti->L[ir][jr];
                    }

                tmpl = new LabImage(Wd, Hd);

            } else {
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for(int ir = 0; ir < Hd; ir += 1)
                    for(int jr = 0; jr < Wd; jr += 1) {
                        orig[ir][jr] = original->L[ir][jr];
                        orig1[ir][jr] = transformed->L[ir][jr];
                    }

                tmpl = new LabImage(transformed->W, transformed->H);


            }

            float minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax;
            ImProcFunctions::MSRLocal(orig, tmpl->L, orig1, Wd, Hd, params->locallab, sk, locRETgainCcurve, 0, 4, 0.8f, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax);
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16)
#endif

            for(int ir = 0; ir < Hd; ir += 1)
                for(int jr = 0; jr < Wd; jr += 1) {
                    tmpl->L[ir][jr] = orig[ir][jr];
                }

            if(!lp.invret) {
                Reti_Local(call, hueplus, huemoins, hueref, dhueret, chromaref, lumaref, lp, deltE, original, transformed, tmpl, cx, cy, 0);
            } else {
                InverseReti_Local(lp, original, transformed, tmpl, cx, cy, 0);
            }

            if(params->locallab.chrrt > 0) {

                if(!lp.invret && call == 2) {
#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16)
#endif

                    for(int ir = 0; ir < Hd; ir += 1)
                        for(int jr = 0; jr < Wd; jr += 1) {

                            orig[ir][jr] = sqrt(SQR(bufreti->a[ir][jr]) + SQR(bufreti->b[ir][jr]));
                            orig1[ir][jr] = sqrt(SQR(bufreti->a[ir][jr]) + SQR(bufreti->b[ir][jr]));
                        }

                } else {

#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16)
#endif

                    for(int ir = 0; ir < GH; ir += 1)
                        for(int jr = 0; jr < GW; jr += 1) {
                            orig[ir][jr] = sqrt(SQR(original->a[ir][jr]) + SQR(original->b[ir][jr]));
                            orig1[ir][jr] = sqrt(SQR(transformed->a[ir][jr]) + SQR(transformed->b[ir][jr]));
                        }
                }

                ImProcFunctions::MSRLocal(orig, tmpl->L, orig1, Wd, Hd, params->locallab, sk, locRETgainCcurve, 1, 4, 0.8f, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax);

                if(!lp.invret && call == 2) {
#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16)
#endif

                    for(int ir = 0; ir < Hd; ir += 1)
                        for(int jr = 0; jr < Wd; jr += 1) {
                            float Chprov = orig1[ir][jr];
                            float2 sincosval;
                            sincosval.y = Chprov == 0.0f ? 1.f : bufreti->a[ir][jr] / Chprov;
                            sincosval.x = Chprov == 0.0f ? 0.f : bufreti->b[ir][jr] / Chprov;
                            tmpl->a[ir][jr] = orig[ir][jr] * sincosval.y;
                            tmpl->b[ir][jr] = orig[ir][jr] * sincosval.x;

                        }


                } else {
#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16)
#endif

                    for(int ir = 0; ir < Hd; ir += 1)
                        for(int jr = 0; jr < Wd; jr += 1) {
                            float Chprov = orig1[ir][jr];
                            float2 sincosval;
                            sincosval.y = Chprov == 0.0f ? 1.f : transformed->a[ir][jr] / Chprov;
                            sincosval.x = Chprov == 0.0f ? 0.f : transformed->b[ir][jr] / Chprov;
                            tmpl->a[ir][jr] = orig[ir][jr] * sincosval.y;
                            tmpl->b[ir][jr] = orig[ir][jr] * sincosval.x;

                        }
                }

                if(!lp.invret) {
                    Reti_Local(call, hueplus, huemoins, hueref, dhueret, chromaref, lumaref, lp, deltE, original, transformed, tmpl, cx, cy, 1);
                } else {
                    InverseReti_Local(lp, original, transformed, tmpl, cx, cy, 1);
                }

            }

            delete tmpl;
            delete [] origBuffer;
            delete [] origBuffer1;

            if(!lp.invret && call == 2) {

                delete  bufreti;
            }
        }


// Gamut and Munsell control - very important do not desactivated to avoid crash
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
                        float chr;

#else
                        float aa = transformed->a[y][x];
                        float bb = transformed->b[y][x];
                        float HH, chr;

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
                        // Color::pregamutlab (Lprov1, HH, chr);
                        Chprov1 = min(Chprov1, chr);

                        Color::gamutLchonly(sincosval, Lprov1, Chprov1, wip, highlight, 0.15f, 0.92f, neg, more_rgb);
#else
                        Color::pregamutlab (Lprov1, HH, chr);
                        Chprov1 = min(Chprov1, chr);
                        Color::gamutLchonly(sincosval, Lprov1, Chprov1, wip, highlight, 0.15f, 0.92f);
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

        if(lp.qualmet == 1) {

            for (int i = 0; i < GH; i++) {
                delete [] deltE[i];
            }

            delete [] deltE;
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
