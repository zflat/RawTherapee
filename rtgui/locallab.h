/*
 *  This file is part of RawTherapee.
 */
#ifndef _LOCALLAB_H_
#define _LOCALLAB_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "edit.h"
#include "../rtengine/coord.h"

class Locallab : public ToolParamBlock, public AdjusterListener, public FoldableToolPanel, public EditSubscriber
{

private:
    int lastObject;
    Gtk::HBox* enaBox;

protected:
//   Gtk::CheckButton* enabled;
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
    Adjuster* radius;
    Adjuster* strength;
    Adjuster* transit;
    Gtk::CheckButton* avoid;
    MyComboBoxText*   Smethod;
    sigc::connection  Smethodconn;
    Gtk::HBox* ctboxS;
    Gtk::CheckButton* invers;
    Gtk::CheckButton* inversrad;

    double draggedPointOldAngle;
    double draggedPointAdjusterAngle;
    double draggedFeatherOffset;
    double draggedlocYOffset;
    double draggedlocXOffset;
    double draggedlocYTOffset;
    double draggedlocXLOffset;
    rtengine::Coord draggedCenter;
    bool lastavoid, lastinvers, lastinversrad;
    sigc::connection  editConn, avoidConn, inversConn, inversradConn;

    void editToggled ();

public:

    Locallab ();
    ~Locallab ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = NULL);
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = NULL);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = NULL);
    void setBatchMode   (bool batchMode);

    void updateGeometry (int centerX_, int centerY_, int locY_, double degree_, int locX_, int locYT_, int locXL_);
    void SmethodChanged      ();

    void adjusterChanged (Adjuster* a, double newval);
    void enabledChanged ();
    void setAdjusterBehavior (bool degreeadd, bool locYadd, bool locXadd, bool locYTadd, bool locXLadd, bool centeradd, bool lightnessadd, bool contrastadd, bool chromaadd, bool sensiadd, bool transitadd, bool radiusadd, bool strengthadd);
    void trimValues          (rtengine::procparams::ProcParams* pp);
    void avoidChanged ();
    void inversChanged ();
    void inversradChanged ();

    void setEditProvider (EditDataProvider* provider);

    // EditSubscriber interface
    CursorShape getCursor(int objectID);
    bool mouseOver(int modifierKey);
    bool button1Pressed(int modifierKey);
    bool button1Released();
    bool drag1(int modifierKey);
    void switchOffEditMode ();
};

#endif
