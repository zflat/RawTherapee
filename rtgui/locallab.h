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

class Locallab : public ToolParamBlock, public AdjusterListener, public FoldableToolPanel, public CurveListener, public EditSubscriber
{

private:
    int lastObject;
    Gtk::HBox* enaBox;

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

    Gtk::CheckButton* avoid;
    MyComboBoxText*   Smethod;
    sigc::connection  Smethodconn;
    Gtk::HBox* ctboxS;
    Gtk::CheckButton* invers;
    Gtk::CheckButton* inversrad;
    Gtk::CheckButton* inversret;

    MyComboBoxText*   retinexMethod;
    Gtk::Label* labmdh;
    Gtk::HBox* dhbox;
    CurveEditorGroup* CCWcurveEditorgainT;
    FlatCurveEditor* cTgainshape;

    double draggedPointOldAngle;
    double draggedPointAdjusterAngle;
    double draggedFeatherOffset;
    double draggedlocYOffset;
    double draggedlocXOffset;
    double draggedlocYTOffset;
    double draggedlocXLOffset;
    rtengine::Coord draggedCenter;
    bool lastavoid, lastinvers, lastinversrad, lastinversret;
    sigc::connection  editConn, avoidConn, inversConn, inversradConn, inversretConn, retinexMethodConn;

    void editToggled ();

public:

    Locallab ();
    ~Locallab ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = NULL);
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = NULL);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = NULL);
    void setBatchMode   (bool batchMode);

    void updateGeometry (const int centerX_, const int centerY_, const int locY_, const double degree_, const int locX_, const int locYT_, const int locXL_, const int fullWidth = -1, const int fullHeight = -1);
    void SmethodChanged      ();

    void adjusterChanged (Adjuster* a, double newval);
    void enabledChanged ();
    void setAdjusterBehavior (bool degreeadd, bool locYadd, bool locXadd, bool locYTadd, bool locXLadd, bool centeradd, bool lightnessadd, bool contrastadd, bool chromaadd, bool sensiadd, bool transitadd, bool radiusadd, bool strengthadd);
    void trimValues          (rtengine::procparams::ProcParams* pp);
    void avoidChanged ();
    void inversChanged ();
    void inversradChanged ();
    void inversretChanged ();
    void curveChanged (CurveEditor* ce);
    void autoOpenCurve ();

    void setEditProvider (EditDataProvider* provider);
    void retinexMethodChanged();

    // EditSubscriber interface
    CursorShape getCursor(int objectID);
    bool mouseOver(int modifierKey);
    bool button1Pressed(int modifierKey);
    bool button1Released();
    bool drag1(int modifierKey);
    void switchOffEditMode ();
};

#endif
