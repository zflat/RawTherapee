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
#include "mytime.h"
#include "iccstore.h"
#include "iccmatrices.h"
#include "color.h"
#include "rt_math.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace rtengine
{

using namespace procparams;

#undef ABS
#undef CLIPS
#undef CLIPC

#define ABS(a) ((a)<0?-(a):(a))
#define CLIPS(a) ((a)>-32768?((a)<32767?(a):32767):-32768)
#define CLIPC(a) ((a)>-32000?((a)<32000?(a):32000):-32000)
#define CLIP2(a) ((a)<MAXVAL ? a : MAXVAL )
#define FCLIP(a) ((a)>0.0?((a)<65535.5?(a):65535.5):0.0)


extern const Settings* settings;

struct local_params {
    float ta, yc, xc;
    float lx, ly;
    float lxL, lyT;
    int h;
    int chro, cont, ligh, sens;
    double rad;
    double stren;
    int trans;
    bool inv;
    bool invrad;
};

static void calcLocalParams(int oW, int oH, const LocallabParams& locallab, struct local_params& lp)
{
    int w = oW;
    int h = oH;
    double local_x = locallab.locX / 200.0;
    double local_y = locallab.locY / 200.0;
    double local_xL = locallab.locXL / 200.0;
    double local_yT = locallab.locYT / 200.0;
    double local_center_x = locallab.centerX / 200.0 + 0.5;
    double local_center_y = locallab.centerY / 200.0 + 0.5;
    bool inverse = locallab.invers;
    int local_chroma = locallab.chroma;
    int local_sensi = locallab.sensi;
    int local_contrast = locallab.contrast;
    int local_lightness = locallab.lightness;
    int local_transit = locallab.transit;
    double radius = locallab.radius;
    bool inverserad = locallab.inversrad;
    double strength = locallab.strength;

    lp.xc = w * local_center_x;
    lp.yc = h * local_center_y;
    lp.lx = w * local_x;
    lp.ly = h * local_y;
    lp.lxL = w * local_xL;
    lp.lyT = h * local_yT;
    lp.chro = local_chroma;
    lp.sens = local_sensi;
    lp.cont = local_contrast;
    lp.ligh = local_lightness;
    lp.trans = local_transit;
    lp.rad = radius;
    lp.stren = strength;
    lp.inv = inverse;
    lp.invrad = inverserad;
}

static float calcLocalFactor(const struct local_params& lp, int lox, int loy, int lcx, int dx, int lcy, int dy, double &factorx, double &factory, int trans)
{
//elipse x2/a2 + y2/b2=1
//transition elipsoidal
//x==>lox y==>loy
// a==> dx  b==>dy
    float ach = (float)trans / 100.f; //transition

    if( (SQR(lox - lcx) / SQR(ach * dx) + SQR(loy - lcy) / SQR(ach * dy)) > 1.f && (SQR(lox - lcx) / SQR(dx) + SQR(loy - lcy) / SQR(dy)) <= 1.f ) {
        float kelip = (float)dx / (float)dy;
        float belip = sqrt((SQR((lox - lcx) / kelip) + SQR(loy - lcy))); //determin position ellipse ==> a and b
        float aelip = belip * kelip;
        float degrad = aelip / dx;
        float aa = 1.f / (ach - 1.f);
        float bb = -aa;
        float ap = M_PI / (1.f - ach);
        float bp = M_PI - ap;
        factorx = 0.5f * (1.f + cos(degrad * ap + bp)); //trigo cos transition
    }
}

static void calc_fourarea (float lox, float loy, float ach, struct local_params& lp, bool &qu1, bool &qu1nt, bool &qu1wt, bool &qu2, bool &qu2nt, bool &qu2wt, bool &qu3, bool &qu3nt, bool &qu3wt, bool &qu4, bool &qu4nt, bool &qu4wt)
{
    //boolean calculation for the four ellipse area
    //one can replace by other "curve" or LUT
    qu1 = (lox >= lp.xc && lox < (lp.xc + lp.lx)) && (loy >= lp.yc && loy < (lp.yc + lp.ly));
    qu1nt = ( (SQR(lox - lp.xc) / SQR(ach * lp.lx) + SQR(loy - lp.yc) / SQR(ach * lp.ly)) < 1.f);
    qu1wt = (((SQR(lox - lp.xc) / SQR(ach * lp.lx) + SQR(loy - lp.yc) / SQR(ach * lp.ly)) > 1.f)  && ((SQR(lox - lp.xc) / SQR(lp.lx) + SQR(loy - lp.yc) / SQR(lp.ly)) < 1.f));

    qu2 = (lox >= lp.xc && lox < (lp.xc + lp.lx)) && (loy < lp.yc && loy > (lp.yc - lp.lyT));
    qu2nt = (SQR(lox - lp.xc) / SQR(ach * lp.lx) + SQR(loy - lp.yc) / SQR(ach * lp.lyT)) < 1.f;
    qu2wt = (((SQR(lox - lp.xc) / SQR(ach * lp.lx) + SQR(loy - lp.yc) / SQR(ach * lp.lyT)) > 1.f)  && ((SQR(lox - lp.xc) / SQR(lp.lx) + SQR(loy - lp.yc) / SQR(lp.lyT)) < 1.f));

    qu3 = ((lox < lp.xc && lox > (lp.xc - lp.lxL)) && (loy <= lp.yc && loy > (lp.yc - lp.lyT)));
    qu3nt = (SQR(lox - lp.xc) / SQR(ach * lp.lxL) + SQR(loy - lp.yc) / SQR(ach * lp.lyT)) < 1.f;
    qu3wt = (((SQR(lox - lp.xc) / SQR(ach * lp.lxL) + SQR(loy - lp.yc) / SQR(ach * lp.lyT)) > 1.f)  && ((SQR(lox - lp.xc) / SQR(lp.lxL) + SQR(loy - lp.yc) / SQR(lp.lyT)) < 1.f));

    qu4 = ((lox < lp.xc && lox > (lp.xc - lp.lxL)) && (loy > lp.yc && loy < (lp.yc + lp.ly)));
    qu4nt = (SQR(lox - lp.xc) / SQR(ach * lp.lxL) + SQR(loy - lp.yc) / SQR(ach * lp.ly)) < 1.f;
    qu4wt = (((SQR(lox - lp.xc) / SQR(ach * lp.lxL) + SQR(loy - lp.yc) / SQR(ach * lp.ly)) > 1.f)  && ((SQR(lox - lp.xc) / SQR(lp.lxL) + SQR(loy - lp.yc) / SQR(lp.ly)) < 1.f));
}

void ImProcFunctions::addGaNoise (LabImage *lab, LabImage *dst, double mean, double variance, int chroma, int sk)
{
//Box-Muller method.
// add noise to image : luma and chroma
// actually chroma disabled
    srand(1);
#define unit_rand() (1.0*rand()/RAND_MAX)
    double u1, u2;

    u1 = 0.0;

    for (int y = 1; y < lab->H; y++) {
        for (int x = 1; x < lab->W; x++) {
            while (u1 == 0.0) {
                u1 = unit_rand();
            }

            //  double X=sqrt(-2.0*variance*log(u1)) * cos(2.*M_PI*u2);//sqrt does not run in debug
            double kvar = 1.;
            double varia;

            if(lab->L[y][x] < 12000.f) {
                kvar = (double)(-0.5f * (lab->L[y][x]) / 12000.f + 1.5f);    //increase effect for low lights < 12000.f
            }

            float ah = -0.5f / 12768.f;
            float bh = 1.f - 20000.f * ah;

            if(lab->L[y][x] > 20000.f) {
                kvar = (double)(ah * (lab->L[y][x]) + bh);    //decrease effect for high lights > 20000.f
            }

            if(kvar < 0.5) {
                kvar = 0.5;
            }

            varia = variance * kvar;
            varia = SQR(varia) / sk;
            double X = pow(-2.0 * varia * log(u1), 0.5) * cos(2.*M_PI * u2);
            double Y = pow(-2.0 * varia * log(u1), 0.5) * sin(2.*M_PI * u2);
            double Xp = mean + sqrt(varia) * X ;
            double Yp = mean + sqrt(varia) * Y ;
            float tempf1 = lab->L[y][x] + Xp;
            float tempa1;
            float tempb1;
            float tempab1;
            float HH;

            float k = 1.f;

            if(chroma == 0) {
                k = 2.f;    //increase chroma noise
            }

            float redPi = 5000.f; //about 32768 / 2*Pi

            if (chroma <= 2) {
                tempab1 = sqrt(SQR(lab->a[y][x]) + SQR(lab->b[y][x])) + k * Xp;

                if(chroma <= 2) {
                    HH = xatan2f(lab->b[y][x], lab->a[y][x]) + k * Xp / redPi;
                } else {
                    HH = xatan2f(lab->b[y][x], lab->a[y][x]);
                }
            }

            if(tempf1 > 32768.f) {
                dst->L[y][x] = 32768.f;
            } else if (tempf1 < 0.f) {
                dst->L[y][x] = 0.f;
            } else {
                dst->L[y][x] = tempf1;
            }

            if(chroma <= 2) {
                if(tempab1 < 0) {
                    dst->a[y][x] = 0.f;
                    dst->b[y][x] = 0.f;
                } else {
                    dst->a[y][x] = tempab1 * cos(HH);
                    dst->b[y][x] = tempab1 * sin(HH);
                }
            }

            float tempf2 = lab->L[y][x] + Yp;
            float tempab2;

            if(chroma <= 2) {
                tempab2 = sqrt(SQR(lab->a[y][x]) + SQR(lab->b[y][x])) + k * Yp;

                if(chroma <= 2) {
                    HH = xatan2f(lab->b[y][x], lab->a[y][x]) + k * Yp / redPi;
                } else {
                    HH = xatan2f(lab->b[y][x], lab->a[y][x]);
                }

            }

            if(tempf2 > 32768.f) {
                dst->L[y - 1][x - 1] = 32768.f;
            } else if (tempf2 < 0.f) {
                dst->L[y - 1][x - 1] = 0.f;
            } else {
                dst->L[y - 1][x - 1] = tempf2;
            }

            if(chroma <= 2) {
                if(tempab2 < 0) {
                    dst->a[y - 1][x - 1] = 0.f;
                    dst->b[y - 1][x - 1] = 0.f;
                } else {
                    dst->a[y - 1][x - 1] = tempab2 * cos(HH);
                    dst->b[y - 1][x - 1] = tempab2 * sin(HH);
                }
            }

            u1 = 0.0;
        }
    }
}

void ImProcFunctions::BlurNoise_Local(struct local_params& lp, LabImage* original, LabImage* transformed, LabImage* tmp1, int sx, int sy, int cx, int cy, int oW, int oH,  int fw, int fh,  LUTf & localcurve, bool locutili, int sk)
{
//local blurr and noise
    #pragma omp parallel for schedule(dynamic,16) if (multiThread)
    for (int y = 0; y < transformed->H; y++) {
        for (int x = 0; x < transformed->W; x++) {
            double factor = 1.;
            double factorx = 1.0;
            double factory = 1.0;
            double fac = 1.f;
            int lox = cx + x;
            int loy = cy + y;
            factorx = 1.f;
            factory = 1.f;
            float ach = (float)lp.trans / 100.f;

            bool qu1, qu1nt, qu1wt, qu2, qu2nt, qu2wt, qu3, qu3nt, qu3wt, qu4, qu4nt, qu4wt;
            calc_fourarea (lox, loy, ach, lp, qu1, qu1nt, qu1wt, qu2, qu2nt, qu2wt, qu3, qu3nt, qu3wt, qu4, qu4nt, qu4wt);

            float l_x, l_y;

            if(qu1) {
                l_x = lp.lx;
                l_y = lp.ly;
            }

            if(qu2) {
                l_x = lp.lx;
                l_y = lp.lyT;
            }

            if(qu3) {
                l_x = lp.lxL;
                l_y = lp.lyT;
            }

            if(qu4) {
                l_x = lp.lxL;
                l_y = lp.ly;
            }

            if(qu1 || qu2  || qu3  || qu4) {
                factorx = 1.f;
                factory = 1.f;
                float difL = tmp1->L[y][x] - original->L[y][x];
                float difa = tmp1->a[y][x] - original->a[y][x];
                float difb = tmp1->b[y][x] - original->b[y][x];
                difL *= factorx;
                difa *= factorx;
                difb *= factorx;

                if(qu1nt || qu2nt  || qu3nt  || qu4nt) {//interior ellipse
                    factorx = 1.f;
                    factory = 1.f;
                    difL *= factorx; //in case of..
                    difa *= factorx;
                    difb *= factorx;
                    transformed->L[y][x] = original->L[y][x] + difL;
                    transformed->a[y][x] = original->a[y][x] + difa;
                    transformed->b[y][x] = original->b[y][x] + difb;
                }
                //intermed transition
                else if (qu1wt || qu2wt  || qu3wt  || qu4wt) {
                    calcLocalFactor(lp, lox, loy, lp.xc, l_x, lp.yc, l_y, factorx, factory, lp.trans);
                    fac = (100.f + factorx * lp.chro) / 100.f; //chroma factor transition
                    difL *= factorx;
                    difa *= factorx;
                    difb *= factorx;
                    transformed->L[y][x] = original->L[y][x] + difL;
                    transformed->a[y][x] = original->a[y][x] + difa;
                    transformed->b[y][x] = original->b[y][x] + difb;
                } else { //exter ellipse full borders no transition
                    factorx = 1.f;
                    factory = 1.f;
                    difL *= factorx;
                    difa *= factorx;
                    difb *= factorx;
                    transformed->L[y][x] = original->L[y][x];
                    transformed->a[y][x] = original->a[y][x];
                    transformed->b[y][x] = original->b[y][x];
                }
            } else {
                transformed->L[y][x] = original->L[y][x];
                transformed->a[y][x] = original->a[y][x];
                transformed->b[y][x] = original->b[y][x];
            }
        }
    }


}

void ImProcFunctions::InverseBlurNoise_Local(struct local_params& lp, LabImage* original, LabImage* transformed, LabImage* tmp1, int sx, int sy, int cx, int cy, int oW, int oH,  int fw, int fh,  LUTf & localcurve, bool locutili, int sk)
{
//inverse local blurr and noise
    #pragma omp parallel for schedule(dynamic,16) if (multiThread)
    for (int y = 0; y < transformed->H; y++) {
        for (int x = 0; x < transformed->W; x++) {
            double factor = 1.;
            double factorx = 1.0;
            double factory = 1.0;
            double fac = 1.f;
            int lox = cx + x;
            int loy = cy + y;
            factorx = 1.f;
            factory = 1.f;
            float difL = tmp1->L[y][x] - original->L[y][x];
            float difa = tmp1->a[y][x] - original->a[y][x];
            float difb = tmp1->b[y][x] - original->b[y][x];
            difL *= factorx;
            difa *= factorx;
            difb *= factorx;
            float ach = (float)lp.trans / 100.f;

            bool qu1, qu1nt, qu1wt, qu2, qu2nt, qu2wt, qu3, qu3nt, qu3wt, qu4, qu4nt, qu4wt;
            calc_fourarea (lox, loy, ach, lp, qu1, qu1nt, qu1wt, qu2, qu2nt, qu2wt, qu3, qu3nt, qu3wt, qu4, qu4nt, qu4wt);

            float l_x, l_y;

            if(qu1) {
                l_x = lp.lx;
                l_y = lp.ly;
            }

            if(qu2) {
                l_x = lp.lx;
                l_y = lp.lyT;
            }

            if(qu3) {
                l_x = lp.lxL;
                l_y = lp.lyT;
            }

            if(qu4) {
                l_x = lp.lxL;
                l_y = lp.ly;
            }


            if(qu1 || qu2  || qu3  || qu4) {
                factorx = 1.f;
                factory = 1.f;
                float difL = tmp1->L[y][x] - original->L[y][x];
                float difa = tmp1->a[y][x] - original->a[y][x];
                float difb = tmp1->b[y][x] - original->b[y][x];
                difL *= factorx;
                difa *= factorx;
                difb *= factorx;

                if(qu1nt || qu2nt  || qu3nt  || qu4nt) {//interior ellipse with max transit no application
                    transformed->L[y][x] = original->L[y][x];
                    transformed->a[y][x] = original->a[y][x];
                    transformed->b[y][x] = original->b[y][x];
                }
                //intermed transition
                else if (qu1wt || qu2wt  || qu3wt  || qu4wt) {
                    calcLocalFactor(lp, lox, loy, lp.xc, l_x, lp.yc, l_y, factorx, factory, lp.trans);
                    factorx = 1.f - factorx;
                    fac = (100.f + factorx * lp.chro) / 100.f; //chroma factor transition
                    difL *= factorx;
                    difa *= factorx;
                    difb *= factorx;

                    transformed->L[y][x] = original->L[y][x] + difL;
                    transformed->a[y][x] = original->a[y][x] + difa;
                    transformed->b[y][x] = original->b[y][x] + difb;
                } else { //exter ellipse full borders no transition
                    factorx = 1.f;
                    factory = 1.f;
                    difL *= factorx;
                    difa *= factorx;
                    difb *= factorx;
                    transformed->L[y][x] = original->L[y][x] + difL;
                    transformed->a[y][x] = original->a[y][x] + difa;
                    transformed->b[y][x] = original->b[y][x] + difb;
                }
            } else {
                transformed->L[y][x] = original->L[y][x] + difL;
                transformed->a[y][x] = original->a[y][x] + difa;
                transformed->b[y][x] = original->b[y][x] + difb;
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
    float DY;
    float dx, dy;
    float ah, bh;
    float al, bl;
};


void ImProcFunctions::Contrast_Local(float pm, bool locL, struct local_contra& lco, float hueref, float dhue, float lumaref, float av, struct local_params& lp, LabImage* original, LabImage* transformed, SHMap* shmap, int sx, int sy, int cx, int cy, int oW, int oH,  int fw, int fh,  LUTf & localcurve, bool locutili, int sk)
{
// contrast - perhaps for 4 areas   if need
// I tried shmap adaptaed to Lab, but no real gain and artifacts
// locL : spot or mean
    #pragma omp parallel for schedule(dynamic,16) if (multiThread)
    for (int y = 0; y < transformed->H; y++) {
        for (int x = 0; x < transformed->W; x++) {
            int lox = cx + x;
            double factorx = 1.0;
            double factory = 1.0;
            int loy = cy + y;
            float conl;
            float fac;
            factorx = 1.f;
            factory = 1.f;
            float diflc;
            float localtype = lumaref; //spot

            if(!locL) {
                localtype = av;    //mean area
            }

            float ach = (float)lp.trans / 100.f;
            float rhue = xatan2f(original->b[y][x], original->a[y][x]);
            float rchro = sqrt(SQR(original->b[y][x]) + SQR(original->a[y][x])) / 327.68f;
            float rL = (original->L[y][x]) / 327.68f;
            float ra = original->a[y][x] / 327.68f;
            float rb = original->b[y][x] / 327.68f;
            float dU = 0.1f; //for transition action
            float realchro = 1.f;
            float dUx = 1.f;
            float hueplus, huemoins, huemoinsred, hueplusred;
            hueplus = hueref + dhue;
            huemoins = hueref - dhue;
            hueplusred = hueref + 0.9f * dhue;
            huemoinsred = hueref - 0.9f * dhue;

            if(hueplusred > M_PI) {
                hueplusred = hueref + 0.9f * dhue - 2.f * M_PI;
            }

            if(huemoinsred < -M_PI) {
                huemoinsred = hueref - 0.9f * dhue + 2.f * M_PI;
            }

            if(hueplus > M_PI) {
                hueplus = hueref + dhue - 2.f * M_PI;
            }

            if(huemoins < -M_PI) {
                huemoins = hueref - dhue + 2.f * M_PI;
            }

            int kzon = 0;
            float kk = 1.f;
            float reducac;

            if(lp.sens < 30.f) {
                reducac = 0.2f * (lp.sens / 100.f);
            } else {
                float areduc = 0.6285714f; //0.44f/0.7f;
                float breduc = 0.5f - areduc;
                reducac = areduc * (lp.sens / 100.f) + breduc;
            }

            float delhu = 0.1f; //between 0.05 and 0.2
            float aplusx = (1.f - lco.dx) / delhu;
            float bplusx = 1.f - aplusx * hueplus;
            float amoinsx = (lco.dx - 1.f) / delhu;
            float bmoinsx = 1.f - amoinsx * huemoins;
            float aplusy = (1.f - lco.dy) / delhu;
            float bplusy = 1.f - aplusy * hueplus;
            float amoinsy = (lco.dy - 1.f) / delhu;
            float bmoinsy = 1.f - amoinsy * huemoins;

            float realcox = lco.dx, realcoy = lco.dy;

            if((hueref + dhue) < M_PI && rhue < hueplus && rhue > huemoins) {//transition leeds to artifacts
                //  if(rhue >= hueplus-delhu  && rhue < hueplus)  {realcox=aplusx*rhue+bplusx;realcoy=aplusy*rhue+bplusy;}//artifacts
                //  else if(rhue >= huemoins && rhue < huemoins+delhu)  {realcox=amoinsx*rhue+bmoinsx;realcoy=amoinsy*rhue+bmoinsy;}
                //  else
                {
                    realcox = lco.dx;
                    realcoy = lco.dy;
                }

                kzon = 1;
            } else if((hueref + dhue) >= M_PI && (rhue > huemoins  || rhue < hueplus )) {
                //  if(rhue >= hueplus-delhu  && rhue < hueplus)  {realcox=aplusx*rhue+bplusx;realcoy=aplusy*rhue+bplusy;}
                //  else if(rhue >= huemoins && rhue < huemoins+delhu)  {realcox=amoinsx*rhue+bmoinsx;realcoy=amoinsy*rhue+bmoinsy;}
                //  else
                {
                    realcox = lco.dx;
                    realcoy = lco.dy;
                }

                kzon = 2;
            }

            if((hueref - dhue) > -M_PI && rhue < hueplus && rhue > huemoins) {
                //  if(rhue >= hueplus-delhu  && rhue < hueplus)  {realcox=aplusx*rhue+bplusx;realcoy=aplusy*rhue+bplusy;}
                //  else if(rhue >= huemoins && rhue < huemoins+delhu)  {realcox=amoinsx*rhue+bmoinsx;realcoy=amoinsy*rhue+bmoinsy;}
                //  else
                {
                    realcox = lco.dx;
                    realcoy = lco.dy;
                }
                kzon = 3;
            } else if((hueref - dhue) <= -M_PI && (rhue > huemoins  || rhue < hueplus )) {
                //  if(rhue >= hueplus-delhu  && rhue < hueplus)  {realcox=aplusx*rhue+bplusx;realcoy=aplusy*rhue+bplusy;}
                //  else if(rhue >= huemoins && rhue < huemoins+delhu)  {realcox=amoinsx*rhue+bmoinsx;realcoy=amoinsy*rhue+bmoinsy;}
                //  else
                {
                    realcox = lco.dx;
                    realcoy = lco.dy;
                }
                kzon = 4;
            }

            lco.alsup = (-realcox) / (localtype / 2.f);
            lco.blsup = -lco.alsup * localtype;
            lco.alsup2 = (realcoy) / (50.f - localtype / 2.f);
            lco.blsup2 = -lco.alsup2 * localtype;
            lco.alsup3 = (realcoy) / (localtype / 2.f - 50.f);
            lco.blsup3 = -lco.alsup3 * 100.f;
            lco.aDY = realcoy;

            lco.alinf = realcox / (localtype / 2.f);
            float vi = (localtype / 2.f) / 100.f;
            float vinf = (50.f + localtype / 2.f) / 100.f;
            ImProcFunctions::secondeg_begin (reducac, vi, lco.aa, lco.bb);//parabolic
            ImProcFunctions::secondeg_end (reducac, vinf, lco.aaa, lco.bbb, lco.ccc);//parabolic

            bool qu1, qu1nt, qu1wt, qu2, qu2nt, qu2wt, qu3, qu3nt, qu3wt, qu4, qu4nt, qu4wt;
            calc_fourarea (lox, loy, ach, lp, qu1, qu1nt, qu1wt, qu2, qu2nt, qu2wt, qu3, qu3nt, qu3wt, qu4, qu4nt, qu4wt);

            float l_x, l_y;

            if(qu1) {
                l_x = lp.lx;
                l_y = lp.ly;
            }

            if(qu2) {
                l_x = lp.lx;
                l_y = lp.lyT;
            }

            if(qu3) {
                l_x = lp.lxL;
                l_y = lp.lyT;
            }

            if(qu4) {
                l_x = lp.lxL;
                l_y = lp.ly;
            }

            if(qu1 || qu2  || qu3  || qu4) {
                factorx = 1.f;
                float ach = (float)lp.trans / 100.f;
                float core;
                float prov100;

                if(qu1nt || qu2nt  || qu3nt  || qu4nt) {//interior ellipse with max transit no application
                    if(original->L[y][x] < 32768.f) {
                        float prov;
                        prov = original->L[y][x] / 327.68f;
                        //  prov=shmap->map[y][x]/327.68f;
                        prov100 = prov / 100.f;

                        //  prov100=shmap->map[y][x]/32768.f;
                        if(original->L[y][x] / 327.68f > localtype) {
                            if(prov >= localtype && prov < 50.f + localtype / 2.f) {
                                core = (lco.alsup2 * prov + lco.blsup2) ;
                                transformed->L[y][x] = 327.68f * (prov + pm * (prov - localtype) * core);
                            } else {
                                //core=(lco.alsup3*prov+lco.blsup3) ;
                                core = lco.aDY * (lco.aaa * prov100 * prov100 + lco.bbb * prov100 + lco.ccc);
                                transformed->L[y][x] = 327.68f * (prov + pm * (prov - localtype) * core);
                            }
                        } else  { //inferior
                            if(prov > localtype / 2.f && prov < localtype)  {
                                core = (lco.alsup * prov + lco.blsup) ;
                                transformed->L[y][x] = 327.68f * (prov - pm * (localtype - prov) * core);
                            } else if(prov <= localtype / 2.f) {
                                core = prov * lco.alinf * (lco.aa * prov100 * prov100 + lco.bb * prov100);
                                //core=lco.alinf*prov;

                                transformed->L[y][x] = 327.68f * (prov - pm * (localtype - prov) * core);
                            }
                        }
                    }
                }
                //intermed transition
                else if (qu1wt || qu2wt  || qu3wt  || qu4wt) {
                    calcLocalFactor(lp, lox, loy, lp.xc, l_x, lp.yc, l_y, factorx, factory, lp.trans);

                    if(original->L[y][x] < 32768.f) {
                        float prov, prov100;

                        if(original->L[y][x] / 327.68f > localtype) {
                            prov = original->L[y][x] / 327.68f;
                            //  prov=shmap->map[y][x]/327.68f;
                            //  prov100=shmap->map[y][x]/32768.f;
                            prov100 = prov / 100.f;

                            if(prov >= localtype && prov < 50.f + localtype / 2.f) {
                                core = (lco.alsup2 * prov + lco.blsup2) ;
                                core *= factorx;
                                transformed->L[y][x] = 327.68f * (prov + pm * (prov - localtype) * (core));
                            } else {
                                //core=(lco.alsup3*prov+lco.blsup3) ;
                                core = lco.aDY * (lco.aaa * prov100 * prov100 + lco.bbb * prov100 + lco.ccc);

                                core *= factorx;
                                transformed->L[y][x] = 327.68f * (prov + pm * (prov - localtype) * (core));
                            }
                        } else  { //inferior
                            prov = original->L[y][x] / 327.68f;
                            //  prov=shmap->map[y][x]/327.68f;

                            prov100 = prov / 100.f;
                            //  prov100=shmap->map[y][x]/32768.f;

                            if(prov > localtype / 2.f && prov < localtype)  {
                                core = (lco.alsup * prov + lco.blsup) ;
                                core *= factorx;
                                transformed->L[y][x] = 327.68f * (prov - pm * (localtype - prov) * core);
                            } else if(prov <= localtype / 2.f) {
                                //  core=lco.alinf*prov;
                                core = prov * lco.alinf * (lco.aa * prov100 * prov100 + lco.bb * prov100);
                                //  core=lco.alinf*prov;

                                core *= factorx;
                                transformed->L[y][x] = 327.68f * (prov - pm * (localtype - prov) * core);
                            }
                        }
                    }
                }
            }

        }
    }

}


void ImProcFunctions::InverseContrast_Local(float ave, float pm, bool locL, struct local_contra& lco, float hueref, float dhue, float lumaref, float av, struct local_params& lp, LabImage* original, LabImage* transformed, LabImage* tmp1, int sx, int sy, int cx, int cy, int oW, int oH,  int fw, int fh,  LUTf & localcurve, bool locutili, int sk)
{
    #pragma omp parallel for schedule(dynamic,16) if (multiThread)

    for (int y = 0; y < transformed->H; y++) {
        for (int x = 0; x < transformed->W; x++) {
            int lox = cx + x;
            double factorx = 1.0;
            double factory = 1.0;
            int loy = cy + y;
            float conl;
            float fac;
            factorx = 1.f;
            factory = 1.f;
            float diflc;
            float localtype = lumaref;

            float ach = (float)lp.trans / 100.f;

            bool qu1, qu1nt, qu1wt, qu2, qu2nt, qu2wt, qu3, qu3nt, qu3wt, qu4, qu4nt, qu4wt;
            calc_fourarea (lox, loy, ach, lp, qu1, qu1nt, qu1wt, qu2, qu2nt, qu2wt, qu3, qu3nt, qu3wt, qu4, qu4nt, qu4wt);

            float l_x, l_y;

            if(qu1) {
                l_x = lp.lx;
                l_y = lp.ly;
            }

            if(qu2) {
                l_x = lp.lx;
                l_y = lp.lyT;
            }

            if(qu3) {
                l_x = lp.lxL;
                l_y = lp.lyT;
            }

            if(qu4) {
                l_x = lp.lxL;
                l_y = lp.ly;
            }

            if(qu1 || qu2  || qu3  || qu4) {
                factorx = 1.f;
                factory = 1.f;

                if(qu1nt || qu2nt  || qu3nt  || qu4nt) {//interior ellipse with max transit no application
                    transformed->L[y][x] = original->L[y][x];
                }
                //intermed transition
                else if (qu1wt || qu2wt  || qu3wt  || qu4wt) {
                    calcLocalFactor(lp, lox, loy, lp.xc, l_x, lp.yc, l_y, factorx, factory, lp.trans);
                    factorx = 1.f - factorx;

                    if(original->L[y][x] < 32768.f) {
                        float kh = lco.ah * (original->L[y][x] / 327.68f) + lco.bh;
                        float kl = lco.al * (original->L[y][x] / 327.68f) + 1.f;
                        float prov;

                        if(original->L[y][x] > ave) {
                            prov = original->L[y][x];
                            original->L[y][x] = ave + kh * (original->L[y][x] - ave);
                        } else  {
                            prov = original->L[y][x];
                            original->L[y][x] = ave - kl * (ave - original->L[y][x]);
                        }

                        float diflc = original->L[y][x] - prov;
                        diflc *= factorx;
                        transformed->L[y][x] =  prov + diflc;
                    }
                } else { //exter ellipse full borders no transition
                    factorx = 1.f;
                    factory = 1.f;

                    if(original->L[y][x] < 32768.f) {
                        float kh = lco.ah * (original->L[y][x] / 327.68f) + lco.bh;
                        float kl = lco.al * (original->L[y][x] / 327.68f) + 1.f;
                        float prov;

                        if(original->L[y][x] > ave) {
                            prov = original->L[y][x];
                            original->L[y][x] = ave + kh * (original->L[y][x] - ave);
                        } else  {
                            prov = original->L[y][x];
                            original->L[y][x] = ave - kl * (ave - original->L[y][x]);
                        }

                        float diflc = original->L[y][x] - prov;
                        diflc *= factorx;
                        transformed->L[y][x] =  prov + diflc;
                    }
                }
            } else {
                factorx = 1.f;
                factory = 1.f;

                if(original->L[y][x] < 32768.f) {
                    float kh = lco.ah * (original->L[y][x] / 327.68f) + lco.bh;
                    float kl = lco.al * (original->L[y][x] / 327.68f) + 1.f;
                    float prov;

                    if(original->L[y][x] > ave) {
                        prov = original->L[y][x];
                        original->L[y][x] = ave + kh * (original->L[y][x] - ave);
                    } else  {
                        prov = original->L[y][x];
                        original->L[y][x] = ave - kl * (ave - original->L[y][x]);
                    }

                    float diflc = original->L[y][x] - prov;
                    diflc *= factorx;
                    transformed->L[y][x] =  prov + diflc;
                }
            }

        }
    }
}

void ImProcFunctions::ColorLight_Local(float hueplus, float huemoins, float hueplusred, float huemoinsred, float hueref, float dhue, float chromaref, float lumaref, struct local_params& lp, LabImage* original, LabImage* transformed, LabImage* tmp1, int sx, int sy, int cx, int cy, int oW, int oH,  int fw, int fh,  LUTf & localcurve, bool locutili, int sk)
{
// chroma and lightness
    #pragma omp parallel for schedule(dynamic,16) if (multiThread)
    for (int y = 0; y < transformed->H; y++) {
        for (int x = 0; x < transformed->W; x++) {
            double factor = 1.;
            double factorx = 1.0;
            double factory = 1.0;
            double fac = 1.f;
            int lox = cx + x;
            int loy = cy + y;
            float ach = (float)lp.trans / 100.f;
            float rhue = xatan2f(original->b[y][x], original->a[y][x]);
            float rchro = sqrt(SQR(original->b[y][x]) + SQR(original->a[y][x])) / 327.68f;
            float rL = (original->L[y][x]) / 327.68f;
            float ra = original->a[y][x] / 327.68f;
            float rb = original->b[y][x] / 327.68f;
            float dU = 0.1f; //for transition action
            float realchro = 1.f;
            float dUx = 1.f;

            float achsens;
            float bchsens;
            //chroma
            float amplchsens = 2.5f;
            achsens = (amplchsens - 1.f) / (100.f - 20.f); //20. default locallab.sensi
            bchsens = 1.f - 20.f * achsens;
            float multchro = lp.sens * achsens + bchsens;

            //luma
            float alumsens;
            float blumsens;
            float ampllumsens = 2.f;
            alumsens = (ampllumsens - 1.f) / (100.f - 20.f); //20. default locallab.sensi
            blumsens = 1.f - 20.f * alumsens;
            float multlum = lp.sens * alumsens + blumsens;

            //skin
            float achsensskin;
            float bchsensskin;
            float amplchsensskin = 1.6f;
            achsensskin = (amplchsensskin - 1.f) / (100.f - 20.f); //20. default locallab.sensi
            bchsensskin = 1.f - 20.f * achsensskin;
            float multchroskin = lp.sens * achsensskin + bchsensskin;
            float kdiff = 0.f;
            int kzon = 0;
            //transition = difficult to avoid artifact with scope on flat area (sky...)
            float delhu = 0.1f; //between 0.05 and 0.2
            float aplus = (1.f - lp.chro) / delhu;
            float bplus = 1.f - aplus * hueplus;
            float amoins = (lp.chro - 1.f) / delhu;
            float bmoins = 1.f - amoins * huemoins;

            if((hueref + dhue) < M_PI && rhue < hueplus && rhue > huemoins) {//transition are good
                if(rhue >= hueplus - delhu  && rhue < hueplus)  {
                    realchro = aplus * rhue + bplus;
                } else if(rhue >= huemoins && rhue < huemoins + delhu)  {
                    realchro = amoins * rhue + bmoins;
                } else {
                    realchro = lp.chro;
                }

                kzon = 1;
            } else if((hueref + dhue) >= M_PI && (rhue > huemoins  || rhue < hueplus )) {
                if(rhue >= hueplus - delhu  && rhue < hueplus)  {
                    realchro = aplus * rhue + bplus;
                } else if(rhue >= huemoins && rhue < huemoins + delhu)  {
                    realchro = amoins * rhue + bmoins;
                } else {
                    realchro = lp.chro;
                }

                kzon = 2;
            }

            if((hueref - dhue) > -M_PI && rhue < hueplus && rhue > huemoins) {
                if(rhue >= hueplus - delhu  && rhue < hueplus)  {
                    realchro = aplus * rhue + bplus;
                } else if(rhue >= huemoins && rhue < huemoins + delhu)  {
                    realchro = amoins * rhue + bmoins;
                } else {
                    realchro = lp.chro;
                }

                kzon = 3;
            } else if((hueref - dhue) <= -M_PI && (rhue > huemoins  || rhue < hueplus )) {
                if(rhue >= hueplus - delhu  && rhue < hueplus)  {
                    realchro = aplus * rhue + bplus;
                } else if(rhue >= huemoins && rhue < huemoins + delhu)  {
                    realchro = amoins * rhue + bmoins;
                } else {
                    realchro = lp.chro;
                }

                kzon = 4;
            }

            if(kzon > 0) {
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

                }
            }

            if(kzon > 0) {
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

            //luma
            float lumdelta = 11.f; //11
            float modlum = lumdelta * multlum;

            if(lumaref + modlum >= 100.f) {
                modlum = (100.f - lumaref) / 2.f;
            }

            if(lumaref - modlum <= 0.f) {
                modlum = (lumaref) / 2.f;
            }

            float alu = 1.f / (lumaref + modlum - 100.f); //linear
            float blu = -100.f * alu; //linear
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

            float kLinf = rL / (100.f);
            float kLsup = rL / (100.f);

            if(kzon > 0/*rhue < hueplus && rhue > huemoins*/) {

                if( (rL > (lumaref - modlum) && rL < (lumaref + modlum))) {
                    kdiff = 1.f;
                } else if (rL > 0.f && rL <= (lumaref - modlum)) {
                    kdiff = aa * kLinf * kLinf + bb * kLinf;    //parabolic
                } else if (rL <= 100.f && rL >= (lumaref + modlum)) {
                    kdiff = aaa * kLsup * kLsup + bbb * kLsup + ccc;    //parabolic
                } else {
                    kdiff = 0.f;
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

            bool qu1, qu1nt, qu1wt, qu2, qu2nt, qu2wt, qu3, qu3nt, qu3wt, qu4, qu4nt, qu4wt;
            calc_fourarea (lox, loy, ach, lp, qu1, qu1nt, qu1wt, qu2, qu2nt, qu2wt, qu3, qu3nt, qu3wt, qu4, qu4nt, qu4wt);


            float l_x, l_y;

            if(qu1) {
                l_x = lp.lx;
                l_y = lp.ly;
            }

            if(qu2) {
                l_x = lp.lx;
                l_y = lp.lyT;
            }

            if(qu3) {
                l_x = lp.lxL;
                l_y = lp.lyT;
            }

            if(qu4) {
                l_x = lp.lxL;
                l_y = lp.ly;
            }

            if(qu1 || qu2  || qu3  || qu4) {
                factorx = 1.f;
                factory = 1.f;
                float lightcont = localcurve[original->L[y][x]]; //apply lightness
                fac = (100.f + factorx * realchro) / 100.f; //chroma factor transition
                float diflc = lightcont - original->L[y][x];
                diflc *= kdiff;
                diflc *= factorx; //transition lightess

                if(qu1nt || qu2nt  || qu3nt  || qu4nt) {//interior ellipse with max transit no application
                    transformed->L[y][x] = original->L[y][x] + diflc;
                    transformed->a[y][x] = original->a[y][x] * fac;
                    transformed->b[y][x] = original->b[y][x] * fac;
                }

                else if (qu1wt || qu2wt  || qu3wt  || qu4wt) {
                    float lightcont = localcurve[original->L[y][x]]; //apply lightness
                    calcLocalFactor(lp, lox, loy, lp.xc, l_x, lp.yc, l_y, factorx, factory, lp.trans);
                    fac = (100.f + factorx * realchro) / 100.f; //chroma factor transition
                    float diflc = lightcont - original->L[y][x];
                    diflc *= kdiff;

                    diflc *= factorx; //transition lightess

                    transformed->L[y][x] = original->L[y][x] + diflc;
                    transformed->a[y][x] = original->a[y][x] * fac;
                    transformed->b[y][x] = original->b[y][x] * fac;
                } else { //exter ellipse full borders no transition
                    transformed->L[y][x] = original->L[y][x];
                    transformed->a[y][x] = original->a[y][x];
                    transformed->b[y][x] = original->b[y][x];
                }
            }
        }
    }
}

void ImProcFunctions::InverseColorLight_Local(float hueplus, float huemoins, float hueplusred, float huemoinsred, float hueref, float dhue, float chromaref, float lumaref, struct local_params& lp, LabImage* original, LabImage* transformed, LabImage* tmp1, int sx, int sy, int cx, int cy, int oW, int oH,  int fw, int fh,  LUTf & localcurve, bool locutili, int sk)
{
    #pragma omp parallel for schedule(dynamic,16) if (multiThread)

    for (int y = 0; y < transformed->H; y++) {
        for (int x = 0; x < transformed->W; x++) {
            double factor = 1.;
            double factorx = 1.0;
            double factory = 1.0;
            double fac = 1.f;
            int lox = cx + x;
            int loy = cy + y;
            factorx = 1.f;
            factory = 1.f;
            float lightcont = localcurve[original->L[y][x]]; //apply lightness
            //  calcLocalFactor(lp, lox, loy, lp.xc, lp.lx, lp.yc, lp.ly, factorx, factory, lp.trans);
            fac = (100.f + factorx * lp.chro) / 100.f; //chroma factor transition
            float diflc = lightcont - original->L[y][x];
            diflc *= factorx; //transition lightess
            float ach = (float)lp.trans / 100.f;
            bool qu1, qu1nt, qu1wt, qu2, qu2nt, qu2wt, qu3, qu3nt, qu3wt, qu4, qu4nt, qu4wt;
            calc_fourarea (lox, loy, ach, lp, qu1, qu1nt, qu1wt, qu2, qu2nt, qu2wt, qu3, qu3nt, qu3wt, qu4, qu4nt, qu4wt);


            float l_x, l_y;

            if(qu1) {
                l_x = lp.lx;
                l_y = lp.ly;
            }

            if(qu2) {
                l_x = lp.lx;
                l_y = lp.lyT;
            }

            if(qu3) {
                l_x = lp.lxL;
                l_y = lp.lyT;
            }

            if(qu4) {
                l_x = lp.lxL;
                l_y = lp.ly;
            }

            if(qu1 || qu2  || qu3  || qu4) {
                factorx = 1.f;
                factory = 1.f;
                float lightcont = localcurve[original->L[y][x]]; //apply lightness
                fac = (100.f + factorx * lp.chro) / 100.f; //chroma factor transition
                float diflc = lightcont - original->L[y][x];
                diflc *= factorx; //transition lightness

                if(qu1nt || qu2nt  || qu3nt  || qu4nt) {//interior ellipse with max transit no application
                    transformed->L[y][x] = original->L[y][x];
                    transformed->a[y][x] = original->a[y][x];
                    transformed->b[y][x] = original->b[y][x];
                }
                //intermed transition
                else if (qu1wt || qu2wt  || qu3wt  || qu4wt) {
                    calcLocalFactor(lp, lox, loy, lp.xc, l_x, lp.yc, l_y, factorx, factory, lp.trans);
                    factorx = 1.f - factorx;
                    fac = (100.f + factorx * lp.chro) / 100.f; //chroma factor transition
                    diflc = lightcont - original->L[y][x];
                    diflc *= factorx;
                    transformed->L[y][x] = original->L[y][x] + diflc;
                    transformed->a[y][x] = original->a[y][x] * fac;
                    transformed->b[y][x] = original->b[y][x] * fac;
                } else { //exter ellipse full borders no transition
                    factorx = 1.f;
                    factory = 1.f;
                    float lightcont = localcurve[original->L[y][x]];
                    fac = (100.f + factorx * lp.chro) / 100.f;
                    float diflc = lightcont - original->L[y][x];
                    diflc *= factorx;
                    transformed->L[y][x] = original->L[y][x] + diflc;
                    transformed->a[y][x] = original->a[y][x] * fac;
                    transformed->b[y][x] = original->b[y][x] * fac;
                }
            } else {
                transformed->L[y][x] = original->L[y][x] + diflc;
                transformed->a[y][x] = original->a[y][x] * fac;
                transformed->b[y][x] = original->b[y][x] * fac;
            }
        }
    }

}


void ImProcFunctions::Lab_Local(LabImage* original, LabImage* transformed, int sx, int sy, int cx, int cy, int oW, int oH,  int fw, int fh,  LUTf & localcurve, bool locutili, int sk)
{
    if(params->locallab.enabled) {
#ifdef _DEBUG
        MyTime t1e, t2e;
        t1e.set();
        // init variables to display Munsell corrections
        MunsellDebugInfo* MunsDebugInfo = new MunsellDebugInfo();
#endif

        struct local_params lp;
        calcLocalParams(oW, oH, params->locallab, lp);

        TMatrix wiprof = iccStore->workingSpaceInverseMatrix (params->icm.working);
        double wip[3][3] = {
            {wiprof[0][0], wiprof[0][1], wiprof[0][2]},
            {wiprof[1][0], wiprof[1][1], wiprof[1][2]},
            {wiprof[2][0], wiprof[2][1], wiprof[2][2]}
        };

        LabImage * tmp1;
        int GW = transformed->W;
        int GH = transformed->H;
        tmp1 = new LabImage(transformed->W, transformed->H);
        float radius = lp.rad / (sk * 2.f); //0 to 50 ==> see skip
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            if(lp.rad > 0.1 || lp.stren > 0.1) {
                AlignedBufferMP<double> buffer(max(GW, GH));
                gaussHorizontal<float> (original->L, tmp1->L, buffer, GW, GH, radius);
                gaussHorizontal<float> (original->a, tmp1->a, buffer, GW, GH, radius);
                gaussHorizontal<float> (original->b, tmp1->b, buffer, GW, GH, radius);
                gaussVertical<float>   (tmp1->L, tmp1->L, buffer, GW, GH, radius);
                gaussVertical<float>   (tmp1->a, tmp1->a, buffer, GW, GH, radius);
                gaussVertical<float>   (tmp1->b, tmp1->b, buffer, GW, GH, radius);
            }
        }

        if(lp.stren > 0.1f) {
            double mean = 0.;//0 best result
            double variance = lp.stren ; //(double) SQR(lp.stren)/sk;
            int chroma = 3;//>=3 disabled chroma
            addGaNoise (tmp1, tmp1, mean, variance, chroma, sk) ;
            printf("vari=%f\n", variance);
        }

        if(lp.rad > 0.1 || lp.stren > 0.1) {
            if(!lp.invrad) { //blurr and noise (center)
                BlurNoise_Local(lp, original, transformed, tmp1, sx, sy,  cx, cy,  oW,  oH,   fw, fh, localcurve, locutili, sk);
            }

            else {
                InverseBlurNoise_Local(lp, original, transformed, tmp1, sx, sy,  cx, cy,  oW,  oH,   fw, fh, localcurve, locutili, sk);
            }
        }

        if(lp.rad > 0.1 || lp.stren > 0.1) {
            delete tmp1;
        }

        // begin map
        SHMap* shmap;

        /*
            double radius = 80.;
            int gW =transformed->W;
            int gH =transformed->H;
            if(!shmap)
            shmap = new SHMap (gW, gH, true);
            shmap->updateLab (original, radius, true, sk);
        */


//  if(shmap)
//     delete shmap; shmap = NULL;
        //


        //begin contrast and evalue hue
        float ave = 0.;
        int n = 0;
        int nab = 0;
        float av, avA, avB, avL;
        float aveA = 0.f;
        float aveB = 0.f;
        float aveL = 0.f;
        float aveSH = 0.f;
        //evauate mean luminance for contrast : actually one area
        // evaluate also hue
        bool inter = true;

        if(inter) {//interior
            //evaluate hue, chroma, luma
            int xred, yred, xredL, yredT;
            xred = 0.02f * oW;
            xredL = xred;
            yred = xred;
            yredT = xred;

            for (int y = 0; y < transformed->H; y++) {
                for (int x = 0; x < transformed->W; x++) {
                    int lox = cx + x;
                    int loy = cy + y;

                    if((lox >= lp.xc && lox < (lp.xc + xred)) && (loy >= lp.yc && loy < (lp.yc + yred))) {
                        aveL += original->L[y][x];
                        aveA += original->a[y][x];
                        aveB += original->b[y][x];
                        //  aveSH+=shmap->map[y][x];

                        nab++;
                    }

                    if((lox >= lp.xc && lox < (lp.xc + xred)) && (loy < lp.yc && loy > (lp.yc - yredT))) {
                        aveL += original->L[y][x];
                        aveA += original->a[y][x];
                        aveB += original->b[y][x];
                        //  aveSH+=shmap->map[y][x];

                        nab++;
                    }

                    if((lox < lp.xc && lox > (lp.xc - xredL)) && (loy <= lp.yc && loy > (lp.yc - yredT))) {
                        aveL += original->L[y][x];
                        aveA += original->a[y][x];
                        aveB += original->b[y][x];
                        //  aveSH+=shmap->map[y][x];

                        nab++;
                    }

                    if((lox < lp.xc && lox > (lp.xc - xredL)) && (loy > lp.yc && loy < (lp.yc + yred))) {
                        aveL += original->L[y][x];
                        aveA += original->a[y][x];
                        aveB += original->b[y][x];
                        //  aveSH+=shmap->map[y][x];
                        nab++;
                    }
                }
            }


            for (int y = 0; y < transformed->H; y++) {
                for (int x = 0; x < transformed->W; x++) {
                    int lox = cx + x;
                    int loy = cy + y;

                    if((lox >= lp.xc && lox < (lp.xc + lp.lx)) && (loy >= lp.yc && loy < (lp.yc + lp.ly))) {
                        ave += original->L[y][x];
                        n++;
                    }

                    if((lox >= lp.xc && lox < (lp.xc + lp.lx)) && (loy < lp.yc && loy > (lp.yc - lp.lyT))) {
                        ave += original->L[y][x];
                        n++;
                    }

                    if((lox < lp.xc && lox > (lp.xc - lp.lxL)) && (loy <= lp.yc && loy > (lp.yc - lp.lyT))) {
                        ave += original->L[y][x];
                        n++;
                    }

                    if((lox < lp.xc && lox > (lp.xc - lp.lxL)) && (loy > lp.yc && loy < (lp.yc + lp.ly))) {
                        ave += original->L[y][x];
                        n++;
                    }
                }
            }
        }

        if (lp.inv) {//exterior
            ave = 0.f;
            n = 0;

            for (int y = 0; y < transformed->H; y++) {
                for (int x = 0; x < transformed->W; x++) {
                    int lox = cx + x;
                    int loy = cy + y;

                    if((lox >= lp.xc && lox < (lp.xc + lp.lx)) && (loy >= lp.yc && loy < (lp.yc + lp.ly))) {
                        //  ave1+=original->L[y][x];
                        //  ave+=original->L[y][x];
                        //  n1++;n++;
                    } else if((lox >= lp.xc && lox < (lp.xc + lp.lx)) && (loy < lp.yc && loy > (lp.yc - lp.lyT))) {
                        //ave2+=original->L[y][x];
                        //  ave+=original->L[y][x];
                        //  n2++;n++;
                    } else if((lox < lp.xc && lox > (lp.xc - lp.lxL)) && (loy <= lp.yc && loy > (lp.yc - lp.lyT))) {
                        //  ave3+=original->L[y][x];
                        //  ave+=original->L[y][x];
                        //  n3++;n++;
                    } else if((lox < lp.xc && lox > (lp.xc - lp.lxL)) && (loy > lp.yc && loy < (lp.yc + lp.ly))) {
                        //  ave4+=original->L[y][x];
                        //  ave+=original->L[y][x];
                        //  n4++;n++;
                    } else {
                        ave += original->L[y][x];
                        //  aveL+=original->L[y][x];
                        //  aveA+=original->a[y][x];
                        //  aveB+=original->b[y][x];

                        n++;
                    }
                }
            }

            //ave=13000.f;
            //n=1;
        }

        if(n == 0) {
            ave = 15000.f;
            n = 1;
        }

        ave = ave / n;

        if(!lp.inv) {//interior
            aveL = aveL / nab;
            aveA = aveA / nab;
            aveB = aveB / nab;
            //  aveSH=aveSH/nab;
        } else {
            aveL = aveL / nab;
            aveA = aveA / nab;
            aveB = aveB / nab;
        }

        av = ave / 327.68f;
        avA = aveA / 327.68f;
        avB = aveB / 327.68f;
        avL = aveL / 327.68f;
//  printf("avL=%f avSH=%f\n",avL,aveSH/327.68f);
        float hueref = xatan2f(avB, avA); //mean hue
        float chromaref = sqrt(SQR(avB) + SQR(avA));
        float lumaref = avL;
//  float lumaref=aveSH/327.68f;
        struct local_contra lco;

        // we must here detect : genearl case, skin, sky,...folliages ???
        // delta dhue, luminance and chroma
        float ared = (M_PI - 0.05f) / 100.f;
        float bred = 0.05f;
        float dhue = ared * lp.sens + bred; //delta hue
        //printf("huref=%f dhue=%f \n",hueref, dhue);

        float maxh = 4.f; //amplification contrast above mean
        float maxl = 3.f; //reductio contrast under mean
        float multh = (float) fabs(lp.cont) * (maxh - 1.f) / 100.f + 1.f;
        float mult = (float)fabs(lp.cont) * (maxl - 1.f) / 100.f + 1.f;
        lco.dx = 1.f - 1.f / mult;
        lco.dy = 1.f - 1.f / multh;
        float  pm = 1.f;

        if(lp.cont < 0.f) {
            pm = -1.f;
        }

        lco.DY = 1.f - 1.f / multh;
        float reducac = 0.2f;

        float prov;
        float multL = (float)lp.cont * (maxl - 1.f) / 100.f + 1.f;
        float multH = (float) lp.cont * (maxh - 1.f) / 100.f + 1.f;

        lco.ah = (multH - 1.f) / (av - 100.f); //av ==> lumaref
        lco.bh = 1.f - 100.f * lco.ah;
        lco.al = (multL - 1.f) / av;
        lco.bl = 1.f;

        bool locL = true;
//  SHMap* shmap;

        if(!lp.inv) {   //contrast interior ellipse
            Contrast_Local(pm, locL, lco, hueref, dhue, lumaref, av, lp, original, transformed, shmap, sx, sy, cx, cy, oW, oH,  fw, fh, localcurve, locutili, sk);
        } else {
            InverseContrast_Local(ave, pm, locL, lco, hueref, dhue, lumaref, av, lp, original, transformed, tmp1, sx, sy, cx, cy, oW, oH,  fw, fh, localcurve, locutili, sk);
        }

// end contrast interior and exterior
        if(shmap) {
            delete shmap;
        }

        shmap = NULL;

        float hueplus, huemoins, huemoinsred, hueplusred;
        hueplus = hueref + dhue;
        huemoins = hueref - dhue;
        hueplusred = hueref + 0.9f * dhue;
        huemoinsred = hueref - 0.9f * dhue;

        if(hueplusred > M_PI) {
            hueplusred = hueref + 0.9f * dhue - 2.f * M_PI;
        }

        if(huemoinsred < -M_PI) {
            huemoinsred = hueref - 0.9f * dhue + 2.f * M_PI;
        }

        if(hueplus > M_PI) {
            hueplus = hueref + dhue - 2.f * M_PI;
        }

        if(huemoins < -M_PI) {
            huemoins = hueref - dhue + 2.f * M_PI;
        }



        if(!lp.inv) { //interior ellipse renforced lightness and chroma
            ColorLight_Local(hueplus, huemoins, hueplusred, huemoinsred, hueref, dhue, chromaref, lumaref, lp, original, transformed, tmp1, sx, sy, cx, cy, oW, oH,  fw, fh, localcurve, locutili, sk);
        }
        //inverse
        else {
            InverseColorLight_Local(hueplus, huemoins, hueplusred, huemoinsred, hueref, dhue, chromaref, lumaref, lp, original, transformed, tmp1, sx, sy, cx, cy, oW, oH,  fw, fh, localcurve, locutili, sk);
        }


// Gamut and Munsell control
        if(params->locallab.avoid) {
#ifdef _DEBUG
            #pragma omp parallel for schedule(dynamic,16) firstprivate(MunsDebugInfo) if (multiThread)
#else
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif

            for (int y = 0; y < transformed->H; y++) {
                for (int x = 0; x < transformed->W; x++) {
                    float HH = xatan2f(transformed->b[y][x], transformed->a[y][x]);
                    float Chprov1 = sqrt(SQR(transformed->a[y][x] / 327.68f) + SQR(transformed->b[y][x] / 327.68f));
                    float Lprov1 = transformed->L[y][x] / 327.68f;
                    float Lprov2 = original->L[y][x] / 327.68f;
                    float memChprov = sqrt(SQR(original->a[y][x] / 327.68f) + SQR(original->b[y][x] / 327.68f));
                    bool highlight = params->toneCurve.hrenabled;
                    float R, G, B;
#ifdef _DEBUG
                    bool neg = false;
                    bool more_rgb = false;
                    Color::gamutLchonly(HH, Lprov1, Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
                    Color::gamutLchonly(HH, Lprov1, Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f);
#endif
                    transformed->L[y][x] = Lprov1 * 327.68f;
                    transformed->a[y][x] = 327.68f * Chprov1 * cos(HH);
                    transformed->b[y][x] = 327.68f * Chprov1 * sin(HH);

                    if (lp.chro != 0.f) {
                        float correctionHue = 0.0f; // Munsell's correction
                        float correctlum = 0.0f;
                        Lprov1 = transformed->L[y][x] / 327.68f;
                        float Chprov = sqrt(SQR(transformed->a[y][x] / 327.68f) + SQR(transformed->b[y][x] / 327.68f));
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
void ImProcFunctions::Lab_Tile(LabImage* lab, LabImage* dst, int ska )
{
    const short int imheight = lab->H, imwidth = lab->W;

    int tilesize;
    int overlap;
//  if(settings->leveldnti ==0) {
//      tilesize = 1024;
//      overlap = 128;
//  }
//  if(settings->leveldnti ==1) {
    //  tilesize = 768;
    //  overlap = 96;
//  }
    tilesize = 512 / ska;
    overlap = 64 / ska;

    LabImage * dsttmp = new LabImage(imwidth, imheight);

    for (int n = 0; n < 3 * imwidth * imheight; n++) {
        dsttmp->data[n] = 0;
    }

    int numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip;
    int kall = 2;
    Tile_calc2 (tilesize, overlap, kall, imwidth, imheight, numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip);
    //now we have tile dimensions, overlaps
#ifdef _OPENMP
    // Calculate number of tiles. If less than omp_get_max_threads(), then limit num_threads to number of tiles
    int numtiles = numtiles_W * numtiles_H;
    int numthreads = MIN(numtiles, omp_get_max_threads());

    if(options.rgbDenoiseThreadLimit > 0) {
        numthreads = MIN(numthreads, options.rgbDenoiseThreadLimit);
    }

    if(numthreads == 1 && omp_get_max_threads() > 1) {
        numthreads = 2;
    }

    #pragma omp parallel num_threads(numthreads)
#endif
    {
        int pos;

#ifdef _OPENMP
        #pragma omp for schedule(dynamic) collapse(2)
#endif

        for (int tiletop = 0; tiletop < imheight; tiletop += tileHskip) {
            for (int tileleft = 0; tileleft < imwidth ; tileleft += tileWskip) {
                pos = (tiletop / tileHskip) * numtiles_W + tileleft / tileWskip ;
                //  printf("pos=%d ", pos);
                int tileright = MIN(imwidth, tileleft + tilewidth);
                int tilebottom = MIN(imheight, tiletop + tileheight);
                int width  = tileright - tileleft;
                int height = tilebottom - tiletop;
                LabImage * labdn = new LabImage(width, height);

                for (int i = tiletop/*, i1=0*/; i < tilebottom; i++/*, i1++*/) {
                    int i1 = i - tiletop;

                    for (int j = tileleft/*, j1=0*/; j < tileright; j++/*, j1++*/) {
                        int j1 = j - tileleft;
                        float L = lab->L[i][j];
                        float a = lab->a[i][j];
                        float b = lab->b[i][j];
                        labdn->L[i1][j1] = L;
                        labdn->a[i1][j1] = a;
                        labdn->b[i1][j1] = b;

                        //here local correction: saturation, lightness, hue, contrast...
                        if(pos > 1 && pos < 8 ) {
                            labdn->a[i1][j1] *= 1.8;
                            labdn->b[i1][j1] *= 1.8;
                        }
                    }
                }

                //calculate mask for feathering output tile overlaps
                float * Vmask = new float [height + 1];
                float * Hmask = new float [width + 1];

                for (int i = 0; i < height; i++) {
                    Vmask[i] = 1;
                }

                for (int j = 0; j < width; j++) {
                    Hmask[j] = 1;
                }

                for (int i = 0; i < overlap; i++) {
                    float mask = SQR(sin((M_PI * i) / (2 * overlap)));

                    if (tiletop > 0) {
                        Vmask[i] = mask;
                    }

                    if (tilebottom < imheight) {
                        Vmask[height - i] = mask;
                    }

                    if (tileleft > 0) {
                        Hmask[i] = mask;
                    }

                    if (tileright < imwidth) {
                        Hmask[width - i] = mask;
                    }
                }

                for (int i = tiletop; i < tilebottom; i++) {
                    int i1 = i - tiletop;
                    float X, Y, Z, L, a, b;

                    for (int j = tileleft; j < tileright; j++) {
                        int j1 = j - tileleft;
                        L = labdn->L[i1][j1];
                        a = labdn->a[i1][j1];
                        b = labdn->b[i1][j1];
                        float factor = Vmask[i1] * Hmask[j1];
                        dsttmp->L[i][j] += factor * L;
                        dsttmp->a[i][j] += factor * a;
                        dsttmp->b[i][j] += factor * b;
                    }
                }

                delete labdn;
                delete[] Vmask;
                delete[] Hmask;

            }
        }
    }
    memcpy (dst->data, dsttmp->data, 3 * dst->W * dst->H * sizeof(float));
    delete dsttmp;

}
void ImProcFunctions::Tile_calc2 (int tilesize, int overlap, int kall, int imwidth, int imheight, int &numtiles_W, int &numtiles_H, int &tilewidth, int &tileheight, int &tileWskip, int &tileHskip)
{

    if(kall == 2) {
        if (imwidth < tilesize) {
            numtiles_W = 1;
            tileWskip = imwidth;
            tilewidth = imwidth;
        } else {
            numtiles_W = ceil(((float)(imwidth)) / (tilesize - overlap));
            tilewidth  = ceil(((float)(imwidth)) / (numtiles_W)) + overlap;
            tilewidth += (tilewidth & 1);
            tileWskip = tilewidth - overlap;
        }

        if (imheight < tilesize) {
            numtiles_H = 1;
            tileHskip = imheight;
            tileheight = imheight;
        } else {
            numtiles_H = ceil(((float)(imheight)) / (tilesize - overlap));
            tileheight = ceil(((float)(imheight)) / (numtiles_H)) + overlap;
            tileheight += (tileheight & 1);
            tileHskip = tileheight - overlap;
        }
    }

    if(kall == 0) {
        numtiles_W = 1;
        tileWskip = imwidth;
        tilewidth = imwidth;
        numtiles_H = 1;
        tileHskip = imheight;
        tileheight = imheight;
    }

    //  printf("Nw=%d NH=%d tileW=%d tileH=%d\n",numtiles_W,numtiles_H,tileWskip,tileHskip);
}





}
