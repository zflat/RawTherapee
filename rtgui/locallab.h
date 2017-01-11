/*
 *  This file is part of RawTherapee.
 */
#ifndef _LOCALLAB_H_
#define _LOCALLAB_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "edit.h"
#include "guiutils.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "toolpanel.h"
#include "../rtengine/imagedata.h"
#include <memory>
#include "options.h"
#include <string>
#include "../rtengine/improcfun.h"


class Locallab : public ToolParamBlock, public AdjusterListener, public FoldableToolPanel, public rtengine::localListener, public CurveListener, public EditSubscriber
{
private:
    int lastObject;

protected:
//   Gtk::CheckButton* enabled;
    Gtk::HBox *editHBox;
    Gtk::ToggleButton* edit;
    Adjuster* degree;
    Adjuster* locX;
    Adjuster* locY;
    Adjuster* locXL;
    Adjuster* locYT;
    Adjuster* centerX;
    Adjuster* centerY;
    Adjuster* circrad;
    Adjuster* lightness;
    Adjuster* contrast;
    Adjuster* chroma;
    Adjuster* sensi;
    Adjuster* sensih;
    Adjuster* radius;
    Adjuster* strength;
    Adjuster* transit;
    Adjuster* str;
    Adjuster* neigh;
    Adjuster* vart;
    Adjuster* chrrt;
    Adjuster* nbspot;
    Adjuster* anbspot;
    Adjuster* maxn;
    Adjuster* sharradius;
    Adjuster* sharamount;
    Adjuster* shardamping;
    Adjuster* shariter;
    Adjuster* sensisha;
    Adjuster* thres;
    Adjuster* proxi;
    Adjuster* noiselumf;
    Adjuster* noiselumc;
    Adjuster* noisechrof;
    Adjuster* noisechroc;
    Adjuster* multiplier[5];
    Adjuster* threshold;
    Adjuster* sensicb;
    Adjuster* sensibn;
    Adjuster* stren;
    Adjuster* gamma;
    Adjuster* estop;
    Adjuster* scaltm;
    Adjuster* rewei;
    Adjuster* sensitm;
    Adjuster* retrab;



    sigc::connection lumaneutralPressedConn;
    sigc::connection lumacontrastPlusPressedConn;
    sigc::connection lumacontrastMinusPressedConn;

    Gtk::CheckButton* avoid;
    MyComboBoxText*   Smethod;
    sigc::connection  Smethodconn;
    Gtk::HBox* ctboxS;
    Gtk::CheckButton* invers;
    Gtk::CheckButton* inversrad;
    Gtk::CheckButton* inversret;
    Gtk::CheckButton* activlum;
    Gtk::CheckButton* inverssha;

    Gtk::Button* neutral;
    Gtk::HBox* neutrHBox;

    Gtk::Button* neutral1;
    Gtk::HBox* neutrHBox1;

    MyComboBoxText*   retinexMethod;
    MyComboBoxText*   qualityMethod;
    Gtk::Label* labmdh;
    Gtk::HBox* dhbox;
    CurveEditorGroup* CCWcurveEditorgainT;
    FlatCurveEditor* cTgainshape;
    CurveEditorGroup* CCWcurveEditorgainTrab;
    FlatCurveEditor* cTgainshaperab;
    CurveEditorGroup* llCurveEditorG;
    DiagonalCurveEditor* llshape;
    Gtk::Image* irg;

    int nextdatasp[59];
    int nextlength;
    std::string nextstr;
    std::string nextstr2;
    std::string nextll_str;
    std::string nextll_str2;

    double draggedPointOldAngle;
    double draggedPointAdjusterAngle;
    double draggedFeatherOffset;
    double draggedlocYOffset;
    double draggedlocXOffset;
    double draggedlocYTOffset;
    double draggedlocXLOffset;
    rtengine::Coord draggedCenter;
    bool lastavoid, lastinvers, lastinversrad, lastinversret, lastactivlum, lastinverssha;
    int lastanbspot;
    sigc::connection  editConn, avoidConn, inversConn, activlumConn, inversradConn, inversretConn, inversshaConn, retinexMethodConn, qualityMethodConn, neutralconn, neutralconn1;

    void editToggled ();

public:

    Locallab ();
    ~Locallab ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);

    void setBatchMode   (bool batchMode);

    void updateGeometry (const int centerX_, const int centerY_, const int circrad_, const int locY_, const double degree_, const int locX_, const int locYT_, const int locXL_, const int fullWidth = -1, const int fullHeight = -1);
    void SmethodChanged      ();

    void adjusterChanged (Adjuster* a, double newval);
    void enabledChanged ();
    void setAdjusterBehavior (bool degreeadd, bool locYadd, bool locXadd, bool locYTadd, bool locXLadd, bool centeradd, bool lightnessadd, bool contrastadd, bool chromaadd, bool sensiadd, bool transitadd, bool radiusadd, bool strengthadd);
    void trimValues          (rtengine::procparams::ProcParams* pp);
    void avoidChanged ();
    void activlumChanged ();
    void inversChanged ();
    void inversradChanged ();
    void inversretChanged ();
    void inversshaChanged ();
    void curveChanged (CurveEditor* ce);
    void autoOpenCurve ();
    void localChanged           (int **datasp, std::string datastr, std::string ll_str, int sp, int maxdat);
    void localretChanged           (int **datasp, std::string datastr, std::string ll_str, int sp, int maxdat);
    bool localComputed_         ();
    bool localretComputed_         ();
    void setEditProvider (EditDataProvider* provider);
    void retinexMethodChanged();
    void qualityMethodChanged();
    void lumaneutralPressed ();
    void lumacontrastPlusPressed ();
    void lumacontrastMinusPressed ();
    void neutral_pressed       ();

    // EditSubscriber interface
    CursorShape getCursor(int objectID);
    bool mouseOver(int modifierKey);
    bool button1Pressed(int modifierKey);
    bool button1Released();
    bool drag1(int modifierKey);
    void switchOffEditMode ();
};

#endif
