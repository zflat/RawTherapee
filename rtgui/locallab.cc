/*
 *  This file is part of RawTherapee.
 */
#include "locallab.h"
#include "rtimage.h"
#include <iomanip>
#include "../rtengine/rt_math.h"

using namespace rtengine;
using namespace rtengine::procparams;

//Locallab::Locallab () : FoldableToolPanel(this), EditSubscriber(ET_OBJECTS), lastObject(-1), draggedPointOldAngle(-1000.)
Locallab::Locallab (): FoldableToolPanel(this, "gradient", M("TP_LOCALLAB_LABEL"), false, true), EditSubscriber(ET_OBJECTS), lastObject(-1), draggedPointOldAngle(-1000.)
{
    editHBox = Gtk::manage (new Gtk::HBox());
    edit = Gtk::manage (new Gtk::ToggleButton());
    edit->add (*Gtk::manage (new RTImage ("editmodehand.png")));
    edit->set_tooltip_text(M("EDIT_OBJECT_TOOLTIP"));
    editConn = edit->signal_toggled().connect( sigc::mem_fun(*this, &Locallab::editToggled) );
    editHBox->pack_start(*edit, Gtk::PACK_SHRINK, 0);
    pack_start (*editHBox, Gtk::PACK_SHRINK, 0);

    Gtk::Frame* shapeFrame = Gtk::manage (new Gtk::Frame (M("TP_LOCALLAB_SHFR")) );
    shapeFrame->set_border_width(0);
    shapeFrame->set_label_align(0.025, 0.5);

    Gtk::Frame* colorFrame = Gtk::manage (new Gtk::Frame (M("TP_LOCALLAB_COFR")) );
    colorFrame->set_border_width(0);
    colorFrame->set_label_align(0.025, 0.5);

    Gtk::Frame* blurrFrame = Gtk::manage (new Gtk::Frame (M("TP_LOCALLAB_BLUFR")) );
    blurrFrame->set_border_width(0);
    blurrFrame->set_label_align(0.025, 0.5);

    Gtk::VBox *shapeVBox = Gtk::manage ( new Gtk::VBox());
    shapeVBox->set_spacing(2);
    shapeVBox->set_border_width(4);

    ctboxS = Gtk::manage (new Gtk::HBox ());
    Gtk::Label* labmS = Gtk::manage (new Gtk::Label (M("TP_LOCALLAB_STYPE") + ":"));
    ctboxS->pack_start (*labmS, Gtk::PACK_SHRINK, 4);
    ctboxS->set_tooltip_markup (M("TP_LOCALLAB_STYPE_TOOLTIP"));

    Smethod = Gtk::manage (new MyComboBoxText ());
//    Smethod->append_text (M("TP_LOCALLAB_IND"));
    Smethod->append_text (M("TP_LOCALLAB_SYM"));
//    Smethod->append_text (M("TP_LOCALLAB_INDSL"));
    Smethod->append_text (M("TP_LOCALLAB_SYMSL"));
    Smethod->set_active(0);
    Smethodconn = Smethod->signal_changed().connect ( sigc::mem_fun(*this, &Locallab::SmethodChanged) );

    locX = Gtk::manage (new Adjuster (M("TP_LOCAL_WIDTH"), 0, 150, 1, 25));
    //locX->set_tooltip_text (M("TP_LOCAL_WIDTH_TOOLTIP"));
    locX->setAdjusterListener (this);

    locXL = Gtk::manage (new Adjuster (M("TP_LOCAL_WIDTH_L"), 0, 150, 1, 25));
    //locX->set_tooltip_text (M("TP_LOCAL_WIDTH_TOOLTIP"));
    locXL->setAdjusterListener (this);

    degree = Gtk::manage (new Adjuster (M("TP_LOCAL_DEGREE"), -180, 180, 1, 0));
    degree->set_tooltip_text (M("TP_LOCAL_DEGREE_TOOLTIP"));
    degree->setAdjusterListener (this);

    locY = Gtk::manage (new Adjuster (M("TP_LOCAL_HEIGHT"), 0, 150, 1, 25));
    //locY->set_tooltip_text (M("TP_LOCAL_HEIGHT_TOOLTIP"));
    locY->setAdjusterListener (this);

    locYT = Gtk::manage (new Adjuster (M("TP_LOCAL_HEIGHT_T"), 0, 150, 1, 25));
    //locY->set_tooltip_text (M("TP_LOCAL_HEIGHT_TOOLTIP"));
    locYT->setAdjusterListener (this);

    centerX = Gtk::manage (new Adjuster (M("TP_LOCALLAB_CENTER_X"), -100, 100, 1, 0));
    //centerX->set_tooltip_text (M("TP_LOCALLAB_CENTER_X_TOOLTIP"));
    centerX->setAdjusterListener (this);

    centerY = Gtk::manage (new Adjuster (M("TP_LOCALLAB_CENTER_Y"), -100, 100, 1, 0));
    //centerY->set_tooltip_text (M("TP_LOCALLAB_CENTER_Y_TOOLTIP"));
    centerY->setAdjusterListener (this);

    lightness = Gtk::manage (new Adjuster (M("TP_LOCALLAB_LIGHTNESS"), -100, 100, 1, 0));
    //lightness->set_tooltip_text (M("TP_LOCALLAB_LIGHTNESS_TOOLTIP"));
    lightness->setAdjusterListener (this);

    contrast = Gtk::manage (new Adjuster (M("TP_LOCALLAB_CONTRAST"), -100, 100, 1, 0));
    //contrast->set_tooltip_text (M("TP_LOCALLAB_CONTRAST_TOOLTIP"));
    contrast->setAdjusterListener (this);

    chroma = Gtk::manage (new Adjuster (M("TP_LOCALLAB_CHROMA"), -100, 200, 1, 0));
    //chroma->set_tooltip_text (M("TP_LOCALLAB_CHROMA_TOOLTIP"));
    chroma->setAdjusterListener (this);

    sensi = Gtk::manage (new Adjuster (M("TP_LOCALLAB_SENSI"), 0, 100, 1, 20));
    sensi->set_tooltip_text (M("TP_LOCALLAB_SENSI_TOOLTIP"));
    sensi->setAdjusterListener (this);

    radius = Gtk::manage ( new Adjuster (M("TP_LOCALLAB_RADIUS"), 0., 100., 0.1, 0.) );
    //radius->set_tooltip_text (M("TP_LOCALLAB_RADIUS_TOOLTIP"));
    radius->setAdjusterListener (this);
    strength = Gtk::manage ( new Adjuster (M("TP_LOCALLAB_STRENGTH"), 0., 100., 0.1, 0.) );
    //radius->set_tooltip_text (M("TP_LOCALLAB_RADIUS_TOOLTIP"));
    strength->setAdjusterListener (this);

    transit = Gtk::manage (new Adjuster (M("TP_LOCALLAB_TRANSIT"), 5, 95, 1, 60));
    transit->set_tooltip_text (M("TP_LOCALLAB_TRANSIT_TOOLTIP"));
    transit->setAdjusterListener (this);

    invers = Gtk::manage (new Gtk::CheckButton (M("TP_LOCALLAB_INVERS")));
    invers->set_active (false);
    inversConn  = invers->signal_toggled().connect( sigc::mem_fun(*this, &Locallab::inversChanged) );

    inversrad = Gtk::manage (new Gtk::CheckButton (M("TP_LOCALLAB_INVERS")));
    inversrad->set_active (false);
    inversradConn  = inversrad->signal_toggled().connect( sigc::mem_fun(*this, &Locallab::inversradChanged) );

    avoid = Gtk::manage (new Gtk::CheckButton (M("TP_LOCALLAB_AVOID")));
    avoid->set_active (false);
    avoidConn  = avoid->signal_toggled().connect( sigc::mem_fun(*this, &Locallab::avoidChanged) );

    ctboxS->pack_start (*Smethod);
    shapeVBox->pack_start (*editHBox);
    shapeVBox->pack_start (*ctboxS);

    shapeVBox->pack_start (*locX);
    shapeVBox->pack_start (*locXL);
    //pack_start (*degree);
    shapeVBox->pack_start (*locY);
    shapeVBox->pack_start (*locYT);
    shapeVBox->pack_start (*centerX);
    shapeVBox->pack_start (*centerY);

    shapeFrame->add(*shapeVBox);
    pack_start (*shapeFrame);

    Gtk::VBox *colorVBox = Gtk::manage ( new Gtk::VBox());
    colorVBox->set_spacing(2);
    colorVBox->set_border_width(4);

    Gtk::VBox *blurrVBox = Gtk::manage ( new Gtk::VBox());
    blurrVBox->set_spacing(2);
    blurrVBox->set_border_width(4);


    colorVBox->pack_start (*lightness);
    colorVBox->pack_start (*contrast);
    colorVBox->pack_start (*chroma);
    colorVBox->pack_start (*sensi);
    colorVBox->pack_start (*invers);

    colorFrame->add(*colorVBox);
    pack_start (*colorFrame);

    blurrVBox->pack_start (*radius);
    blurrVBox->pack_start (*strength);
    blurrVBox->pack_start (*inversrad);
    blurrFrame->add(*blurrVBox);
    pack_start (*blurrFrame);

    pack_start (*transit);
    pack_start (*avoid);

    // Instantiating the Editing geometry; positions will be initialized later
    Line  *hLine, *vLine, *locYLine[2], *locXLine[2];
    Circle *centerCircle;

    // Visible geometry
    locXLine[0] = new Line();
    locXLine[0]->innerLineWidth = 2;
    locXLine[1] = new Line();
    locXLine[1]->innerLineWidth = 2;
    locXLine[0]->datum  = locXLine[1]->datum = Geometry::IMAGE;

    locYLine[0] = new Line();
    locYLine[0]->innerLineWidth = 2;
    locYLine[1] = new Line();
    locYLine[1]->innerLineWidth = 2;
    locYLine[0]->datum = locYLine[1]->datum = Geometry::IMAGE;

    centerCircle = new Circle();
    centerCircle->datum = Geometry::IMAGE;
    centerCircle->radiusInImageSpace = false;
    centerCircle->radius = 8;
    centerCircle->filled = false;

    EditSubscriber::visibleGeometry.push_back( locXLine[0] );
    EditSubscriber::visibleGeometry.push_back( locXLine[1] );
    EditSubscriber::visibleGeometry.push_back( locYLine[0] );
    EditSubscriber::visibleGeometry.push_back( locYLine[1] );
    EditSubscriber::visibleGeometry.push_back( centerCircle );

    // MouseOver geometry
    locXLine[0] = new Line();
    locXLine[0]->innerLineWidth = 2;
    locXLine[1] = new Line();
    locXLine[1]->innerLineWidth = 2;
    locXLine[0]->datum  = locXLine[1]->datum = Geometry::IMAGE;

    locYLine[0] = new Line();
    locYLine[0]->innerLineWidth = 2;
    locYLine[1] = new Line();
    locYLine[1]->innerLineWidth = 2;
    locYLine[0]->datum = locYLine[1]->datum = Geometry::IMAGE;

    centerCircle = new Circle();
    centerCircle->datum = Geometry::IMAGE;
    centerCircle->radiusInImageSpace = false;
    centerCircle->radius = 8;
    centerCircle->filled = false;

    EditSubscriber::mouseOverGeometry.push_back( locXLine[0] );
    EditSubscriber::mouseOverGeometry.push_back( locXLine[1] );

    EditSubscriber::mouseOverGeometry.push_back( locYLine[0] );
    EditSubscriber::mouseOverGeometry.push_back( locYLine[1] );

    EditSubscriber::mouseOverGeometry.push_back( centerCircle );

    show_all();
}

Locallab::~Locallab()
{
    for (std::vector<Geometry*>::const_iterator i = visibleGeometry.begin(); i != visibleGeometry.end(); ++i) {
        delete *i;
    }

    for (std::vector<Geometry*>::const_iterator i = mouseOverGeometry.begin(); i != mouseOverGeometry.end(); ++i) {
        delete *i;
    }
}

void Locallab::read (const ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();

    if (pedited) {
        degree->setEditedState (pedited->locallab.degree ? Edited : UnEdited);
        locY->setEditedState (pedited->locallab.locY ? Edited : UnEdited);
        locX->setEditedState (pedited->locallab.locX ? Edited : UnEdited);
        locYT->setEditedState (pedited->locallab.locYT ? Edited : UnEdited);
        locXL->setEditedState (pedited->locallab.locXL ? Edited : UnEdited);
        centerX->setEditedState (pedited->locallab.centerX ? Edited : UnEdited);
        centerY->setEditedState (pedited->locallab.centerY ? Edited : UnEdited);
        lightness->setEditedState (pedited->locallab.lightness ? Edited : UnEdited);
        contrast->setEditedState (pedited->locallab.contrast ? Edited : UnEdited);
        chroma->setEditedState (pedited->locallab.chroma ? Edited : UnEdited);
        sensi->setEditedState (pedited->locallab.sensi ? Edited : UnEdited);
        radius->setEditedState (pedited->locallab.radius ? Edited : UnEdited);
        strength->setEditedState (pedited->locallab.strength ? Edited : UnEdited);
        transit->setEditedState (pedited->locallab.transit ? Edited : UnEdited);
//      enabled->set_inconsistent (multiImage && !pedited->locallab.enabled);
        set_inconsistent (multiImage && !pedited->locallab.enabled);
        avoid->set_inconsistent (multiImage && !pedited->locallab.avoid);
        invers->set_inconsistent (multiImage && !pedited->locallab.invers);
        inversrad->set_inconsistent (multiImage && !pedited->locallab.inversrad);

        if (!pedited->locallab.Smethod) {
            Smethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

    }

    Smethodconn.block(true);

    //enaConn.block (true);
    //enabled->set_active (pp->locallab.enabled);
    //enaConn.block (false);
    setEnabled(pp->locallab.enabled);
    avoidConn.block (true);
    avoid->set_active (pp->locallab.avoid);
    avoidConn.block (false);
    inversConn.block (true);
    invers->set_active (pp->locallab.invers);
    inversConn.block (false);
    inversradConn.block (true);
    inversrad->set_active (pp->locallab.inversrad);
    inversradConn.block (false);

    degree->setValue (pp->locallab.degree);
    locY->setValue (pp->locallab.locY);
    locX->setValue (pp->locallab.locX);
    locYT->setValue (pp->locallab.locYT);
    locXL->setValue (pp->locallab.locXL);
    centerX->setValue (pp->locallab.centerX);
    centerY->setValue (pp->locallab.centerY);
    lightness->setValue (pp->locallab.lightness);
    contrast->setValue (pp->locallab.contrast);
    chroma->setValue (pp->locallab.chroma);
    sensi->setValue (pp->locallab.sensi);
    transit->setValue (pp->locallab.transit);
    radius->setValue (pp->locallab.radius);
    strength->setValue (pp->locallab.strength);

//  lastEnabled = pp->locallab.enabled;
    lastavoid = pp->locallab.avoid;
    lastinvers = pp->locallab.invers;
    lastinversrad = pp->locallab.inversrad;
    inversChanged();
    updateGeometry (pp->locallab.centerX, pp->locallab.centerY, pp->locallab.locY, pp->locallab.degree,  pp->locallab.locX, pp->locallab.locYT, pp->locallab.locXL);

//   if (pp->locallab.Smethod == "IND") {
//       Smethod->set_active (0);
//   } else
    if (pp->locallab.Smethod == "SYM") {
        Smethod->set_active (0);
//   } else if (pp->locallab.Smethod == "INDSL") {
//       Smethod->set_active (2);
    } else if (pp->locallab.Smethod == "SYMSL") {
        Smethod->set_active (1);
    }

    SmethodChanged();
    Smethodconn.block(false);

    if (pp->locallab.Smethod == "SYM" || pp->locallab.Smethod == "SYMSL") {
        locXL->setValue (locX->getValue());
        locYT->setValue (locY->getValue());
    } else if (pp->locallab.Smethod == "LOC") {
        locXL->setValue (locX->getValue());
        locYT->setValue (locX->getValue());
        locY->setValue (locX->getValue());
    } else if (pp->locallab.Smethod == "INDSL" || pp->locallab.Smethod == "IND") {
        locX->setValue (pp->locallab.locX);
        locY->setValue (pp->locallab.locY);
        locXL->setValue (pp->locallab.locXL);
        locYT->setValue (pp->locallab.locYT);

    }

    enableListener ();
}

void Locallab::updateGeometry(const int centerX_, const int centerY_, const int locY_, const double degree_, const int locX_, const int locYT_, const int locXL_, const int fullWidth, const int fullHeight)
{
    EditDataProvider* dataProvider = getEditProvider();


    if (!dataProvider) {
        return;
    }

    int imW = 0;
    int imH = 0;

    if (fullWidth != -1 && fullHeight != -1) {
        imW = fullWidth;
        imH = fullHeight;
    } else {
        dataProvider->getImageSize(imW, imH);

        if (!imW || !imH) {
            return;
        }
    }

    PolarCoord polCoord1, polCoord2, polCoord0;
    // dataProvider->getImageSize(imW, imH);
    double decayY = (locY_) * double(imH) / 200.;
    double decayYT = (locYT_) * double(imH) / 200.;
    double decayX = (locX_) * (double(imW)) / 200.;
    double decayXL = (locXL_) * (double(imW)) / 200.;
    rtengine::Coord origin(imW / 2 + centerX_ * imW / 200.f, imH / 2 + centerY_ * imH / 200.f);
    printf("deX=%f dexL=%f deY=%f deyT=%f\n", decayX, decayXL, decayY, decayYT);

    //  if (Smethod->get_active_row_number()==2) decayY=decayYT=decayXL=decayX;
//   if (Smethod->get_active_row_number() == 1 || Smethod->get_active_row_number() == 3) {
    if (Smethod->get_active_row_number() == 0 || Smethod->get_active_row_number() == 1) {
        decayYT = decayY;
        decayXL = decayX;
    }

    Line *currLine;
    Circle *currCircle;
    double decay;
    const auto updateLine = [&](Geometry * geometry, const float radius, const float begin, const float end) {
        const auto line = static_cast<Line*>(geometry);
        line->begin = PolarCoord(radius, -degree_ + begin);
        line->begin += origin;
        line->end = PolarCoord(radius, -degree_ + end);
        line->end += origin;
    };

    const auto updateLineWithDecay = [&](Geometry * geometry, const float radius, const float decal, const float offSetAngle) {
        const auto line = static_cast<Line*>(geometry);//180
        line->begin = PolarCoord (radius, -degree_ + decal) + PolarCoord (decay, -degree_ + offSetAngle);
        line->begin += origin;//0
        line->end = PolarCoord (radius, -degree_ + (decal - 180)) + PolarCoord (decay, -degree_ + offSetAngle);
        line->end += origin;
    };

    const auto updateCircle = [&](Geometry * geometry) {
        const auto circle = static_cast<Circle*>(geometry);
        circle->center = origin;
    };

    decay = decayX;
    updateLineWithDecay (visibleGeometry.at(0), 100., 90., 0.);
    updateLineWithDecay (mouseOverGeometry.at(0), 100., 90., 0.);

    decay = decayXL;

    updateLineWithDecay (visibleGeometry.at(1), 100., 90., 180.);
    updateLineWithDecay (mouseOverGeometry.at(1), 100., 90., 180.);

    decay = decayYT;
    updateLineWithDecay (visibleGeometry.at(2), 100., 180., 270.);
    updateLineWithDecay (mouseOverGeometry.at(2), 100., 180., 270.);

    decay = decayY;

    updateLineWithDecay (visibleGeometry.at(3), 100., 180, 90.);
    updateLineWithDecay (mouseOverGeometry.at(3), 100., 180., 90.);


    updateCircle (visibleGeometry.at(4));
    updateCircle (mouseOverGeometry.at(4));


}

void Locallab::write (ProcParams* pp, ParamsEdited* pedited)
{
    pp->locallab.degree = degree->getValue ();
    pp->locallab.locY = locY->getIntValue ();
    pp->locallab.locX = locX->getValue ();
    pp->locallab.locYT = locYT->getIntValue ();
    pp->locallab.locXL = locXL->getValue ();
    pp->locallab.centerX = centerX->getIntValue ();
    pp->locallab.centerY = centerY->getIntValue ();
    pp->locallab.lightness = lightness->getIntValue ();
    pp->locallab.contrast = contrast->getIntValue ();
    pp->locallab.chroma = chroma->getIntValue ();
    pp->locallab.sensi = sensi->getIntValue ();
    pp->locallab.radius = radius->getValue ();
    pp->locallab.strength = strength->getValue ();
    //pp->locallab.enabled = enabled->get_active();
    pp->locallab.enabled = getEnabled();
    pp->locallab.transit = transit->getIntValue ();
    pp->locallab.avoid = avoid->get_active();
    pp->locallab.invers = invers->get_active();
    pp->locallab.inversrad = inversrad->get_active();

    if (pedited) {
        pedited->locallab.degree = degree->getEditedState ();
        pedited->locallab.Smethod  = Smethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->locallab.locY = locY->getEditedState ();
        pedited->locallab.locX = locX->getEditedState ();
        pedited->locallab.locYT = locYT->getEditedState ();
        pedited->locallab.locXL = locXL->getEditedState ();
        pedited->locallab.centerX = centerX->getEditedState ();
        pedited->locallab.centerY = centerY->getEditedState ();
        pedited->locallab.lightness = lightness->getEditedState ();
        pedited->locallab.contrast = contrast->getEditedState ();
        pedited->locallab.chroma = chroma->getEditedState ();
        pedited->locallab.sensi = sensi->getEditedState ();
        pedited->locallab.radius = radius->getEditedState ();
        pedited->locallab.strength = strength->getEditedState ();
        pedited->locallab.transit = transit->getEditedState ();
//      pedited->locallab.enabled = !enabled->get_inconsistent();
        pedited->locallab.enabled = !get_inconsistent();
        pedited->locallab.avoid = !avoid->get_inconsistent();
        pedited->locallab.invers = !invers->get_inconsistent();
        pedited->locallab.inversrad = !inversrad->get_inconsistent();
    }

//    if (Smethod->get_active_row_number() == 0) {
//       pp->locallab.Smethod = "IND";
//   } else
    if (Smethod->get_active_row_number() == 0) {
        pp->locallab.Smethod = "SYM";
        //  } else if (Smethod->get_active_row_number() == 2) {
//      pp->locallab.Smethod = "INDSL";
    } else if (Smethod->get_active_row_number() == 1) {
        pp->locallab.Smethod = "SYMSL";
    }

//   if(Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
    if(Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 1) {
        pp->locallab.locX = locX->getValue();
        pp->locallab.locY = locY->getValue();

        pp->locallab.locXL = pp->locallab.locX;
        pp->locallab.locYT = pp->locallab.locY;
    }
    /*  else if(Smethod->get_active_row_number()==2){
            pp->locallab.locXL=pp->locallab.locX;
            pp->locallab.locYT=pp->locallab.locX;
            pp->locallab.locY=pp->locallab.locX;
        }
        */
    else {
        pp->locallab.locXL = locXL->getValue();
        pp->locallab.locX = locX->getValue();
        pp->locallab.locY = locY->getValue();
        pp->locallab.locYT = locYT->getValue();
    }
}

void Locallab::SmethodChanged ()
{
    if (!batchMode) {
        if(Smethod->get_active_row_number() == 5) { //IND 0
            locX->hide();
            locXL->hide();
            locY->hide();
            locYT->hide();
            centerX->hide();
            centerY->hide();
        } else if(Smethod->get_active_row_number() == 0) {          // 1 SYM
            locX->hide();
            locXL->hide();
            locY->hide();
            locYT->hide();
            centerX->hide();
            centerY->hide();

        } else if(Smethod->get_active_row_number() == 2) {          //2 SYM
            locX->show();
            locXL->show();
            locY->show();
            locYT->show();
            centerX->show();
            centerY->show();

        } else if(Smethod->get_active_row_number() == 1) {          // 3 SYM
            locX->show();
            locXL->hide();
            locY->show();
            locYT->hide();
            centerX->show();
            centerY->show();

        }

        /*      else if(Smethod->get_active_row_number()==2) {              // LOC
                    locX->show();
                    locXL->hide();
                    locY->hide();
                    locYT->hide();
                }   */
    }

    //if (listener && ()) ) {
    if (listener && getEnabled()) {
        //  if(Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
        if(Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 1) {
            listener->panelChanged (EvlocallabSmet, Smethod->get_active_text ());
            locXL->setValue (locX->getValue());
            locYT->setValue (locY->getValue());
        }
        //   else if(Smethod->get_active_row_number()==2) {
        //          listener->panelChanged (EvlocallabSmet, Smethod->get_active_text ());
        //           locXL->setValue (locX->getValue());
        //           locYT->setValue (locX->getValue());
        //          locY->setValue (locX->getValue());
        //     }
        else

        {
            listener->panelChanged (EvlocallabSmet, Smethod->get_active_text ());

        }
    }
}
void Locallab::inversChanged ()
{

    if (batchMode) {
        if (invers->get_inconsistent()) {
            invers->set_inconsistent (false);
            inversConn.block (true);
            invers->set_active (false);
            inversConn.block (false);
        } else if (lastinvers) {
            invers->set_inconsistent (true);
        }

        lastinvers = invers->get_active ();
    }

    if(invers->get_active ()) {
        sensi->hide();
    } else {
        sensi->show();
    }

    if (listener) {
        if (getEnabled()) {
            listener->panelChanged (Evlocallabinvers, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (Evlocallabinvers, M("GENERAL_DISABLED"));
        }
    }
}

void Locallab::inversradChanged ()
{

    if (batchMode) {
        if (inversrad->get_inconsistent()) {
            inversrad->set_inconsistent (false);
            inversradConn.block (true);
            inversrad->set_active (false);
            inversradConn.block (false);
        } else if (lastinversrad) {
            inversrad->set_inconsistent (true);
        }

        lastinversrad = inversrad->get_active ();
    }

    if (listener) {
        if (getEnabled()) {
            listener->panelChanged (Evlocallabinversrad, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (Evlocallabinversrad, M("GENERAL_DISABLED"));
        }
    }
}


void Locallab::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{
    degree->setDefault (defParams->locallab.degree);
    locY->setDefault (defParams->locallab.locY);
    locX->setDefault (defParams->locallab.locX);
    locYT->setDefault (defParams->locallab.locYT);
    locXL->setDefault (defParams->locallab.locXL);
    centerX->setDefault (defParams->locallab.centerX);
    centerY->setDefault (defParams->locallab.centerY);
    lightness->setDefault (defParams->locallab.lightness);
    contrast->setDefault (defParams->locallab.contrast);
    chroma->setDefault (defParams->locallab.chroma);
    sensi->setDefault (defParams->locallab.sensi);
    transit->setDefault (defParams->locallab.transit);
    radius->setDefault (defParams->locallab.radius);
    strength->setDefault (defParams->locallab.strength);


    if (pedited) {
        degree->setDefaultEditedState (pedited->locallab.degree ? Edited : UnEdited);
        locY->setDefaultEditedState (pedited->locallab.locY ? Edited : UnEdited);
        locX->setDefaultEditedState (pedited->locallab.locX ? Edited : UnEdited);
        locYT->setDefaultEditedState (pedited->locallab.locYT ? Edited : UnEdited);
        locXL->setDefaultEditedState (pedited->locallab.locXL ? Edited : UnEdited);
        centerX->setDefaultEditedState (pedited->locallab.centerX ? Edited : UnEdited);
        centerY->setDefaultEditedState (pedited->locallab.centerY ? Edited : UnEdited);
        lightness->setDefaultEditedState (pedited->locallab.lightness ? Edited : UnEdited);
        contrast->setDefaultEditedState (pedited->locallab.contrast ? Edited : UnEdited);
        chroma->setDefaultEditedState (pedited->locallab.chroma ? Edited : UnEdited);
        sensi->setDefaultEditedState (pedited->locallab.sensi ? Edited : UnEdited);
        radius->setDefaultEditedState (pedited->locallab.radius ? Edited : UnEdited);
        strength->setDefaultEditedState (pedited->locallab.strength ? Edited : UnEdited);
        transit->setDefaultEditedState (pedited->locallab.transit ? Edited : UnEdited);
    } else {
        degree->setDefaultEditedState (Irrelevant);
        locY->setDefaultEditedState (Irrelevant);
        locX->setDefaultEditedState (Irrelevant);
        locYT->setDefaultEditedState (Irrelevant);
        locXL->setDefaultEditedState (Irrelevant);
        centerX->setDefaultEditedState (Irrelevant);
        centerY->setDefaultEditedState (Irrelevant);
        lightness->setDefaultEditedState (Irrelevant);
        contrast->setDefaultEditedState (Irrelevant);
        chroma->setDefaultEditedState (Irrelevant);
        sensi->setDefaultEditedState (Irrelevant);
        radius->setDefaultEditedState (Irrelevant);
        strength->setDefaultEditedState (Irrelevant);
        transit->setDefaultEditedState (Irrelevant);
    }
}

void Locallab::adjusterChanged (Adjuster* a, double newval)
{

    updateGeometry (int(centerX->getValue()), int(centerY->getValue()), (int)locY->getValue(), degree->getValue(), (int)locX->getValue(), (int)locYT->getValue(), (int)locXL->getValue() );

    if (listener && getEnabled()) {
        if (a == degree) {
            listener->panelChanged (EvlocallabDegree, degree->getTextValue());
        } else if (a == locY) {
            if(Smethod->get_active_row_number() == 4  || Smethod->get_active_row_number() == 2) {  // 0 2
                // if(Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 2) {  // 0 2
                listener->panelChanged (EvlocallablocY, locY->getTextValue());
            }
            /*  else if(Smethod->get_active_row_number()==2) {
                    listener->panelChanged (EvlocallablocY, locY->getTextValue());
                    locXL->setValue (locX->getValue());
                    locY->setValue (locX->getValue());
                    locYT->setValue (locX->getValue());
                    }*/
            //else if(Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
            else if(Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 1) {
                listener->panelChanged (EvlocallablocY, locY->getTextValue());
                locYT->setValue (locY->getValue());
            }
        } else if (a == locX) {
            //listener->panelChanged (EvlocallablocX, locX->getTextValue());
            // if(Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 2) {
            if(Smethod->get_active_row_number() == 4  || Smethod->get_active_row_number() == 2) {
                listener->panelChanged (EvlocallablocX, locX->getTextValue());
            }
            /*  else if(Smethod->get_active_row_number()==2) {
                    listener->panelChanged (EvlocallablocX, locX->getTextValue());
                    locXL->setValue (locX->getValue());
                    locY->setValue (locX->getValue());
                    locYT->setValue (locX->getValue());
                    }*/
            //   else if(Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
            else if(Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 1) {
                listener->panelChanged (EvlocallablocX, locX->getTextValue());
                locXL->setValue (locX->getValue());
            }
        } else if (a == locYT) {
            // if(Smethod->get_active_row_number() == 0 || Smethod->get_active_row_number() == 2) {
            if(Smethod->get_active_row_number() == 4 || Smethod->get_active_row_number() == 2) {
                listener->panelChanged (EvlocallablocYT, locYT->getTextValue());
            }
            /*  else if(Smethod->get_active_row_number()==2) {
                    listener->panelChanged (EvlocallablocYT, locYT->getTextValue());
                    locXL->setValue (locX->getValue());
                    locY->setValue (locX->getValue());
                    locYT->setValue (locX->getValue());
                    }*/
            //  else if(Smethod->get_active_row_number() == 1 || Smethod->get_active_row_number() == 3) {
            else if(Smethod->get_active_row_number() == 0 || Smethod->get_active_row_number() == 1) {
                listener->panelChanged (EvlocallablocYT, locYT->getTextValue());
                locYT->setValue (locY->getValue());
            }
        } else if (a == locXL) {
            //  if(Smethod->get_active_row_number() == 0 || Smethod->get_active_row_number() == 2) {
            if(Smethod->get_active_row_number() == 4 || Smethod->get_active_row_number() == 2) {
                listener->panelChanged (EvlocallablocXL, locXL->getTextValue());
            }
            /*  else if(Smethod->get_active_row_number()==2) {
                    listener->panelChanged (EvlocallablocXL, locXL->getTextValue());
                    locXL->setValue (locX->getValue());
                    locY->setValue (locX->getValue());
                    locYT->setValue (locX->getValue());
                    }*/
            //   else if(Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
            else if(Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 1) {
                listener->panelChanged (EvlocallablocXL, locXL->getTextValue());
                locXL->setValue (locX->getValue());
            }
        } else if (a == lightness) {
            listener->panelChanged (Evlocallablightness, lightness->getTextValue());
        } else if (a == contrast) {
            listener->panelChanged (Evlocallabcontrast, contrast->getTextValue());
        } else if (a == chroma) {
            listener->panelChanged (Evlocallabchroma, chroma->getTextValue());
        } else if (a == sensi) {
            listener->panelChanged (Evlocallabsensi, sensi->getTextValue());
        } else if (a == radius) {
            listener->panelChanged (Evlocallabradius, radius->getTextValue());
        } else if (a == strength) {
            listener->panelChanged (Evlocallabstrength, strength->getTextValue());
        } else if (a == transit) {
            listener->panelChanged (Evlocallabtransit, transit->getTextValue());
        }

        else if (a == centerX || a == centerY) {
            listener->panelChanged (EvlocallabCenter, Glib::ustring::compose ("X=%1\nY=%2", centerX->getTextValue(), centerY->getTextValue()));
        }
    }
}

void Locallab::enabledChanged ()
{

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvlocallabEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvlocallabEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvlocallabEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void Locallab::avoidChanged ()
{

    if (batchMode) {
        if (avoid->get_inconsistent()) {
            avoid->set_inconsistent (false);
            avoidConn.block (true);
            avoid->set_active (false);
            avoidConn.block (false);
        } else if (lastavoid) {
            avoid->set_inconsistent (true);
        }

        lastavoid = avoid->get_active ();
    }

    if (listener) {
        if (getEnabled()) {
            listener->panelChanged (Evlocallabavoid, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (Evlocallabavoid, M("GENERAL_DISABLED"));
        }
    }
}

void Locallab::setAdjusterBehavior (bool degreeadd, bool locYadd, bool locXadd, bool locYTadd, bool locXLadd, bool centeradd, bool lightnessadd, bool contrastadd, bool chromaadd, bool sensiadd, bool transitadd, bool radiusadd,  bool strengthadd)
{
    degree->setAddMode(degreeadd);
    locY->setAddMode(locYadd);
    locX->setAddMode(locXadd);
    locYT->setAddMode(locYTadd);
    locXL->setAddMode(locXLadd);
    centerX->setAddMode(centeradd);
    centerY->setAddMode(centeradd);
    lightness->setAddMode(lightnessadd);
    contrast->setAddMode(contrastadd);
    chroma->setAddMode(chromaadd);
    sensi->setAddMode(sensiadd);
    transit->setAddMode(transitadd);
    radius->setAddMode(radiusadd);
    strength->setAddMode(strengthadd);

}

void Locallab::trimValues (rtengine::procparams::ProcParams* pp)
{
    degree->trimValue(pp->locallab.degree);
    locY->trimValue(pp->locallab.locY);
    locX->trimValue(pp->locallab.locX);
    locYT->trimValue(pp->locallab.locYT);
    locXL->trimValue(pp->locallab.locXL);
    centerX->trimValue(pp->locallab.centerX);
    centerY->trimValue(pp->locallab.centerY);
    lightness->trimValue(pp->locallab.lightness);
    contrast->trimValue(pp->locallab.contrast);
    chroma->trimValue(pp->locallab.chroma);
    sensi->trimValue(pp->locallab.sensi);
    radius->trimValue(pp->locallab.radius);
    strength->trimValue(pp->locallab.strength);
    transit->trimValue(pp->locallab.transit);
}

void Locallab::setBatchMode (bool batchMode)
{
    //removeIfThere(enaBox, edit, false);
    removeIfThere(this, edit, false);
    ToolPanel::setBatchMode (batchMode);
    degree->showEditedCB ();
    locY->showEditedCB ();
    locX->showEditedCB ();
    locYT->showEditedCB ();
    locXL->showEditedCB ();
    centerX->showEditedCB ();
    centerY->showEditedCB ();
    lightness->showEditedCB ();
    contrast->showEditedCB ();
    chroma->showEditedCB ();
    sensi->showEditedCB ();
    radius->showEditedCB ();
    strength->showEditedCB ();
    transit->showEditedCB ();
    Smethod->append_text (M("GENERAL_UNCHANGED"));

}

void Locallab::setEditProvider (EditDataProvider* provider)
{
    EditSubscriber::setEditProvider(provider);
}

void Locallab::editToggled ()
{
    if (edit->get_active()) {
        subscribe();
    } else {
        unsubscribe();
    }
}

CursorShape Locallab::getCursor(int objectID)
{
    switch (objectID) {
    case (2): {
        int angle = degree->getIntValue();

        if (angle < -135 || (angle >= -45 && angle <= 45) || angle > 135) {
            return CSMove1DV;
        }

        return CSMove1DH;
    }

    case (3): {
        int angle = degree->getIntValue();

        if (angle < -135 || (angle >= -45 && angle <= 45) || angle > 135) {
            return CSMove1DV;
        }

        return CSMove1DH;
    }

    case (0): {
        int angle = degree->getIntValue();

        if (angle < -135 || (angle >= -45 && angle <= 45) || angle > 135) {
            return CSMove1DH;
        }

        return CSMove1DV;
    }

    case (1): {
        int angle = degree->getIntValue();

        if (angle < -135 || (angle >= -45 && angle <= 45) || angle > 135) {
            return CSMove1DH;
        }

        return CSMove1DV;
    }

    case (4):
        return CSMove2D;

    default:
        return CSOpenHand;
    }
}

bool Locallab::mouseOver(int modifierKey)
{
    EditDataProvider* editProvider = getEditProvider();

    if (editProvider && editProvider->object != lastObject) {
        if (lastObject > -1) {
            if (lastObject == 2 || lastObject == 3) {
                EditSubscriber::visibleGeometry.at(2)->state = Geometry::NORMAL;
                EditSubscriber::visibleGeometry.at(3)->state = Geometry::NORMAL;
            } else if (lastObject == 0 || lastObject == 1) {
                EditSubscriber::visibleGeometry.at(0)->state = Geometry::NORMAL;
                EditSubscriber::visibleGeometry.at(1)->state = Geometry::NORMAL;
            }

            else {
                EditSubscriber::visibleGeometry.at(lastObject)->state = Geometry::NORMAL;
            }
        }

        if (editProvider->object > -1) {
            if (editProvider->object == 2 || editProvider->object == 3) {
                EditSubscriber::visibleGeometry.at(2)->state = Geometry::PRELIGHT;
                EditSubscriber::visibleGeometry.at(3)->state = Geometry::PRELIGHT;
            } else if (editProvider->object == 0 || editProvider->object == 1) {
                EditSubscriber::visibleGeometry.at(0)->state = Geometry::PRELIGHT;
                EditSubscriber::visibleGeometry.at(1)->state = Geometry::PRELIGHT;
            }

            else {
                EditSubscriber::visibleGeometry.at(editProvider->object)->state = Geometry::PRELIGHT;
            }
        }

        lastObject = editProvider->object;
        return true;
    }

    return false;
}

bool Locallab::button1Pressed(int modifierKey)
{
    if (lastObject < 0) {
        return false;
    }

    EditDataProvider *provider = getEditProvider();

    if (!(modifierKey & GDK_CONTROL_MASK)) {
        // button press is valid (no modifier key)
        PolarCoord pCoord;
        //  EditDataProvider *provider = getEditProvider();
        int imW, imH;
        provider->getImageSize(imW, imH);
        double halfSizeW = imW / 2.;
        double halfSizeH = imH / 2.;
        draggedCenter.set(int(halfSizeW + halfSizeW * (centerX->getValue() / 100.)), int(halfSizeH + halfSizeH * (centerY->getValue() / 100.)));

        // trick to get the correct angle (clockwise/counter-clockwise)
        rtengine::Coord p1 = draggedCenter;
        rtengine::Coord p2 = provider->posImage;
        int p = p1.y;
        p1.y = p2.y;
        p2.y = p;
        pCoord = p2 - p1;
        draggedPointOldAngle = pCoord.angle;
        draggedPointAdjusterAngle = degree->getValue();

        // if(Smethod->get_active_row_number() == 0 || Smethod->get_active_row_number() == 2) {
        if(Smethod->get_active_row_number() == 4 || Smethod->get_active_row_number() == 2) {
            if (lastObject == 2) {
                PolarCoord draggedPoint;
                rtengine::Coord currPos;
                currPos = provider->posImage;
                rtengine::Coord centerPos = draggedCenter;
                double verti = double(imH);
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint = currPos - centerPos;
                // compute the projected value of the dragged point
                draggedlocYOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue()) / 180.*M_PI);

                if (lastObject == 2) {
                    //draggedlocYOffset = -draggedlocYOffset;
                    draggedlocYOffset -= (locYT->getValue() / 200. * verti);
                }
            } else if (lastObject == 3) {
                // Dragging a line to change the angle
                PolarCoord draggedPoint;
                rtengine::Coord currPos;
                currPos = provider->posImage;
                rtengine::Coord centerPos = draggedCenter;

                double verti = double(imH);

                // trick to get the correct angle (clockwise/counter-clockwise)
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint = currPos - centerPos;

                // draggedPoint.setFromCartesian(centerPos, currPos);
                // compute the projected value of the dragged point
                draggedlocYOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue()) / 180.*M_PI);

                if (lastObject == 3) {
                    draggedlocYOffset = -draggedlocYOffset;
                    draggedlocYOffset -= (locY->getValue() / 200. * verti);
                }

            }

            //   } else if(Smethod->get_active_row_number() == 1 || Smethod->get_active_row_number() == 3) {
        } else if(Smethod->get_active_row_number() == 0 || Smethod->get_active_row_number() == 1) {
            if (lastObject == 2 || lastObject == 3) {
                // Dragging a line to change the angle
                PolarCoord draggedPoint;
                rtengine::Coord currPos;
                currPos = provider->posImage;
                rtengine::Coord centerPos = draggedCenter;
                double verti = double(imH);
                // trick to get the correct angle (clockwise/counter-clockwise)
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint = currPos - centerPos;

                //    draggedPoint.setFromCartesian(centerPos, currPos);
                // compute the projected value of the dragged point
                draggedlocYOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue()) / 180.*M_PI);

                if (lastObject == 3) {
                    draggedlocYOffset = -draggedlocYOffset;
                }

                draggedlocYOffset -= (locY->getValue() / 200. * verti);
            }
        }

        // if(Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 2) {
        if(Smethod->get_active_row_number() == 4  || Smethod->get_active_row_number() == 2) {
            //else if (lastObject==0) {
            if (lastObject == 0) {
                // Dragging a line to change the angle
                PolarCoord draggedPoint;
                rtengine::Coord currPos;
                currPos = provider->posImage;
                rtengine::Coord centerPos = draggedCenter;

                double horiz = double(imW);

                // trick to get the correct angle (clockwise/counter-clockwise)
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint = currPos - centerPos;

                //     draggedPoint.setFromCartesian(centerPos, currPos);
                // compute the projected value of the dragged point
                printf("rad=%f ang=%f\n", draggedPoint.radius, draggedPoint.angle - degree->getValue());
                draggedlocXOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue() + 90.) / 180.*M_PI);
                //  if (lastObject==1)
                //      draggedlocXOffset = -draggedlocXOffset;//-
                draggedlocXOffset -= (locX->getValue() / 200. * horiz);
            } else if (lastObject == 1) {
                // Dragging a line to change the angle
                PolarCoord draggedPoint;
                rtengine::Coord currPos;
                currPos = provider->posImage;
                rtengine::Coord centerPos = draggedCenter;
                double horiz = double(imW);
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint = currPos - centerPos;

                //     draggedPoint.setFromCartesian(centerPos, currPos);
                printf("rad=%f ang=%f\n", draggedPoint.radius, draggedPoint.angle - degree->getValue());
                draggedlocXOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue() + 90.) / 180.*M_PI);

                if (lastObject == 1) {
                    draggedlocXOffset = -draggedlocXOffset;    //-
                }

                draggedlocXOffset -= (locXL->getValue() / 200. * horiz);
            }

            //     } else if(Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
        } else if(Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 1) {

            if (lastObject == 0 || lastObject == 1) {
                PolarCoord draggedPoint;
                rtengine::Coord currPos;
                currPos = provider->posImage;
                rtengine::Coord centerPos = draggedCenter;
                double horiz = double(imW);
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint = currPos - centerPos;

                //    draggedPoint.setFromCartesian(centerPos, currPos);
                printf("rad=%f ang=%f\n", draggedPoint.radius, draggedPoint.angle - degree->getValue());
                draggedlocXOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue() + 90.) / 180.*M_PI);

                if (lastObject == 1) {
                    draggedlocXOffset = -draggedlocXOffset;    //-
                }

                draggedlocXOffset -= (locX->getValue() / 200. * horiz);
            }
        }

        /*  else if(Smethod->get_active_row_number()==2) {
                if (lastObject==0 || lastObject==1 || lastObject==2 || lastObject==3) {
                if (lastObject==2 || lastObject==3) {
                    // Dragging a line to change the angle
                    PolarCoord draggedPoint;
                    Coord currPos;
                    currPos = provider->posImage;
                    Coord centerPos = draggedCenter;
                    double verti = double(imH);
                    // trick to get the correct angle (clockwise/counter-clockwise)
                    int p = centerPos.y;
                    centerPos.y = currPos.y;
                    currPos.y = p;

                    draggedPoint.setFromCartesian(centerPos, currPos);
                    // compute the projected value of the dragged point
                    draggedlocYOffset = draggedPoint.radius * sin((draggedPoint.angle-degree->getValue())/180.*M_PI);
                    if (lastObject==3)
                        draggedlocYOffset = -draggedlocYOffset;
                    draggedlocYOffset -= (locY->getValue() / 200. * verti);
                }


                if (lastObject==0 || lastObject==1) {
                    PolarCoord draggedPoint;
                    Coord currPos;
                    currPos = provider->posImage;
                    Coord centerPos = draggedCenter;
                    double horiz = double(imW);
                    int p = centerPos.y;
                    centerPos.y = currPos.y;
                    currPos.y = p;
                    draggedPoint.setFromCartesian(centerPos, currPos);
                    printf("rad=%f ang=%f\n",draggedPoint.radius,draggedPoint.angle-degree->getValue());
                    draggedlocXOffset = draggedPoint.radius * sin((draggedPoint.angle-degree->getValue()+90.)/180.*M_PI);
                    if (lastObject==1)
                        draggedlocXOffset = -draggedlocXOffset;//-
                    draggedlocXOffset -= (locX->getValue() / 200. * horiz);
                }

                }
            }
            */
        //    EditSubscriber::dragging = true;
        EditSubscriber::action = ES_ACTION_DRAGGING;
        return false;
    } else {
        // this will let this class ignore further drag events
        if (lastObject > -1) { // should theoretically always be true
            if (lastObject == 2 || lastObject == 3) {
                EditSubscriber::visibleGeometry.at(2)->state = Geometry::NORMAL;
                EditSubscriber::visibleGeometry.at(3)->state = Geometry::NORMAL;
            }

            if (lastObject == 0 || lastObject == 1) {
                EditSubscriber::visibleGeometry.at(0)->state = Geometry::NORMAL;
                EditSubscriber::visibleGeometry.at(1)->state = Geometry::NORMAL;
            } else {
                EditSubscriber::visibleGeometry.at(lastObject)->state = Geometry::NORMAL;
            }
        }

        lastObject = -1;
        return true;
    }
}

bool Locallab::button1Released()
{
    draggedPointOldAngle = -1000.;
    EditSubscriber::action = ES_ACTION_NONE;

    return true;
}

bool Locallab::drag1(int modifierKey)
{
    // compute the polar coordinate of the mouse position
    EditDataProvider *provider = getEditProvider();
    int imW, imH;
    provider->getImageSize(imW, imH);
    double halfSizeW = imW / 2.;
    double halfSizeH = imH / 2.;

    //if(Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 2) {
    if(Smethod->get_active_row_number() == 4  || Smethod->get_active_row_number() == 2) {
        if (lastObject == 2) {
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            rtengine::Coord currPos;
            currPos = provider->posImage + provider->deltaImage;
            rtengine::Coord centerPos = draggedCenter;
            double verti = double(imH);
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint = currPos - centerPos;

            //  draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedlocYOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue()) / 180.*M_PI);

            if (lastObject == 2) {
                currDraggedlocYOffset -= draggedlocYOffset;
            }

            //else if (lastObject==3)
            // Dragging the lower locY bar
            //  currDraggedlocYOffset = -currDraggedlocYOffset + draggedlocYOffset;
            currDraggedlocYOffset = currDraggedlocYOffset * 200. / verti;

            if (int(currDraggedlocYOffset) != locYT->getIntValue()) {
                locYT->setValue((int(currDraggedlocYOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();
                updateGeometry (centX, centY, locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

                if (listener) {
                    listener->panelChanged (EvlocallablocY, locYT->getTextValue());
                }

                return true;
            }
        } else if (lastObject == 3) {
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            rtengine::Coord currPos;
            currPos = provider->posImage + provider->deltaImage;
            rtengine::Coord centerPos = draggedCenter;
            double verti = double(imH);
            // trick to get the correct angle (clockwise/counter-clockwise)
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint = currPos - centerPos;

            //  draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedlocYOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue()) / 180.*M_PI);

            //  if (lastObject==2)
            // Dragging the upper locY bar
            //      currDraggedlocYOffset -= draggedlocYOffset;
            //  else
            if (lastObject == 3)
                // Dragging the lower locY bar
            {
                currDraggedlocYOffset = -currDraggedlocYOffset + draggedlocYOffset;
            }

            currDraggedlocYOffset = currDraggedlocYOffset * 200. / verti;

            if (int(currDraggedlocYOffset) != locY->getIntValue()) {

                locY->setValue((int(currDraggedlocYOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();

                updateGeometry (centX, centY, locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

                if (listener) {
                    listener->panelChanged (EvlocallablocY, locY->getTextValue());
                }

                return true;
            }
        }

//   } else if(Smethod->get_active_row_number() == 1 || Smethod->get_active_row_number() == 3) {
    } else if(Smethod->get_active_row_number() == 0 || Smethod->get_active_row_number() == 1) {
        if (lastObject == 2 || lastObject == 3) {
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            rtengine::Coord currPos;
            currPos = provider->posImage + provider->deltaImage;
            rtengine::Coord centerPos = draggedCenter;
            double verti = double(imH);
            // trick to get the correct angle (clockwise/counter-clockwise)
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint = currPos - centerPos;

            //   draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedlocYOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue()) / 180.*M_PI);

            if (lastObject == 2)
                // Dragging the upper locY bar
            {
                currDraggedlocYOffset -= draggedlocYOffset;
            } else if (lastObject == 3)
                // Dragging the lower locY bar
            {
                currDraggedlocYOffset = -currDraggedlocYOffset + draggedlocYOffset;
            }

            currDraggedlocYOffset = currDraggedlocYOffset * 200. / verti;

            if (int(currDraggedlocYOffset) != locY->getIntValue()) {
                locY->setValue((int(currDraggedlocYOffset)));
                //Smethod->get_active_row_number()==2
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();

                updateGeometry (centX, centY, locY->getValue(), degree->getValue(), locX->getValue(),  locYT->getValue(), locXL->getValue());

                if (listener) {
                    //   if(Smethod->get_active_row_number() == 1 || Smethod->get_active_row_number() == 3) {
                    if(Smethod->get_active_row_number() == 0 || Smethod->get_active_row_number() == 1) {
                        listener->panelChanged (EvlocallablocY, locY->getTextValue());
                    }

                    //  else listener->panelChanged (EvlocallablocY, locX->getTextValue());

                }

                return true;
            }
        }

    }

    //  if(Smethod->get_active_row_number() == 0 || Smethod->get_active_row_number() == 2) {
    if(Smethod->get_active_row_number() == 4 || Smethod->get_active_row_number() == 2) {
        //else if (lastObject==0) {
        if (lastObject == 0) {
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            rtengine::Coord currPos;
            currPos = provider->posImage + provider->deltaImage;
            rtengine::Coord centerPos = draggedCenter;
            double horiz = double(imW);
            // trick to get the correct angle (clockwise/counter-clockwise)
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint = currPos - centerPos;

            //    draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedStrOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue() + 90.) / 180.*M_PI);

            if (lastObject == 0)
                // Dragging the upper locY bar
            {
                currDraggedStrOffset -= draggedlocXOffset;
            } else if (lastObject == 1)
                // Dragging the lower locY bar
            {
                currDraggedStrOffset = - currDraggedStrOffset - draggedlocXOffset;    //-
            }

            currDraggedStrOffset = currDraggedStrOffset * 200. / horiz;

            if (int(currDraggedStrOffset) != locX->getIntValue()) {
                locX->setValue((int(currDraggedStrOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();
                updateGeometry (centX, centY, locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

                if (listener) {
                    listener->panelChanged (EvlocallablocX, locX->getTextValue());
                }

                return true;
            }
        } else if (lastObject == 1) {
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            rtengine::Coord currPos;
            currPos = provider->posImage + provider->deltaImage;
            rtengine::Coord centerPos = draggedCenter;
            double horiz = double(imW);
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint = currPos - centerPos;

            //draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedStrOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue() + 90.) / 180.*M_PI);

            if (lastObject == 0)
                // Dragging the upper locY bar
            {
                currDraggedStrOffset -= draggedlocXOffset;
            } else if (lastObject == 1)
                // Dragging the lower locY bar
            {
                currDraggedStrOffset = - currDraggedStrOffset - draggedlocXOffset;    //-
            }

            currDraggedStrOffset = currDraggedStrOffset * 200. / horiz;

            if (int(currDraggedStrOffset) != locXL->getIntValue()) {
                locXL->setValue((int(currDraggedStrOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();
                updateGeometry (centX, centY, locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

                if (listener) {
                    listener->panelChanged (EvlocallablocX, locX->getTextValue());
                }

                return true;
            }
        }

        //  } else if(Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
    } else if(Smethod->get_active_row_number() == 0  || Smethod->get_active_row_number() == 1) {
        if (lastObject == 0 || lastObject == 1) {
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            rtengine::Coord currPos;
            currPos = provider->posImage + provider->deltaImage;
            rtengine::Coord centerPos = draggedCenter;
            double horiz = double(imW);
            // trick to get the correct angle (clockwise/counter-clockwise)
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint = currPos - centerPos;

            // draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedStrOffset = draggedPoint.radius * sin((draggedPoint.angle - degree->getValue() + 90.) / 180.*M_PI);

            if (lastObject == 0)
                // Dragging the upper locY bar
            {
                currDraggedStrOffset -= draggedlocXOffset;
            } else if (lastObject == 1)
                // Dragging the lower locY bar
            {
                currDraggedStrOffset = - currDraggedStrOffset - draggedlocXOffset;    //-
            }

            currDraggedStrOffset = currDraggedStrOffset * 200. / horiz;

            if (int(currDraggedStrOffset) != locX->getIntValue()) {
                locX->setValue((int(currDraggedStrOffset)));
                double centX, centY;
                centX = centerX->getValue();
                centY = centerY->getValue();
                updateGeometry (centX, centY, locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

                if (listener) {
                    listener->panelChanged (EvlocallablocX, locX->getTextValue());
                }

                return true;
            }
        }
    }

    /*  else if(Smethod->get_active_row_number()==2) {
            if (lastObject==0 || lastObject==1 || lastObject==2 || lastObject==3) {
        if (lastObject==2 || lastObject==3) {
            // Dragging the upper or lower locY bar
            PolarCoord draggedPoint;
            Coord currPos;
            currPos = provider->posImage+provider->deltaImage;
            Coord centerPos = draggedCenter;
            double verti = double(imH);
            // trick to get the correct angle (clockwise/counter-clockwise)
            int p = centerPos.y;
            centerPos.y = currPos.y;
            currPos.y = p;
            draggedPoint.setFromCartesian(centerPos, currPos);
            double currDraggedlocYOffset = draggedPoint.radius * sin((draggedPoint.angle-degree->getValue())/180.*M_PI);
            double currDraggedStrOffset = draggedPoint.radius * sin((draggedPoint.angle-degree->getValue() +90.)/180.*M_PI);

            if (lastObject==2)
                currDraggedlocYOffset -= draggedlocYOffset;
            else if (lastObject==3)
                currDraggedlocYOffset = -currDraggedlocYOffset + draggedlocYOffset;
            currDraggedlocYOffset = currDraggedlocYOffset * 200. / verti;
        //  if (int(currDraggedlocYOffset) != locY->getIntValue()) {
        //      locY->setValue((int(currDraggedlocYOffset)));
            if (int(currDraggedlocYOffset) != locX->getIntValue()) {//locX
        //  if (int(currDraggedStrOffset) != locX->getIntValue()) {//locX
                locX->setValue((int(currDraggedlocYOffset)));
                double centX,centY;
                centX=centerX->getValue();
                centY=centerY->getValue();

            //  updateGeometry (centX, centY, locY->getValue(), degree->getValue(), locX->getValue(),  locYT->getValue(), locXL->getValue());
                updateGeometry (centX, centY, locX->getValue(), degree->getValue(), locX->getValue(),  locX->getValue(), locX->getValue());
                if (listener) {
                    if(Smethod->get_active_row_number()==1) listener->panelChanged (EvlocallablocY, locY->getTextValue());

                    }
                return true;
            }
        }
            if (lastObject==0 || lastObject==1) {
                // Dragging the upper or lower locY bar
                PolarCoord draggedPoint;
                Coord currPos;
                currPos = provider->posImage+provider->deltaImage;
                Coord centerPos = draggedCenter;
                double horiz = double(imW);
                int p = centerPos.y;
                centerPos.y = currPos.y;
                currPos.y = p;
                draggedPoint.setFromCartesian(centerPos, currPos);
                double currDraggedStrOffset = draggedPoint.radius * sin((draggedPoint.angle-degree->getValue() +90.)/180.*M_PI);
                if (lastObject==0)
                    currDraggedStrOffset -= draggedlocXOffset;
                else if (lastObject==1)
                    currDraggedStrOffset = - currDraggedStrOffset - draggedlocXOffset;//-
                    currDraggedStrOffset = currDraggedStrOffset * 200. / horiz;

                if (int(currDraggedStrOffset) != locX->getIntValue()) {
                    locX->setValue((int(currDraggedStrOffset)));
                    double centX,centY;
                    centX=centerX->getValue();
                    centY=centerY->getValue();
                    updateGeometry (centX, centY, locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(),locXL->getValue());
                    if (listener)
                        listener->panelChanged (EvlocallablocX, locX->getTextValue());
                    return true;
                }
            }


            }
        }
        */
    //else if (lastObject==4) {
    if (lastObject == 4) {

        // Dragging the circle to change the center
        rtengine::Coord currPos;
        draggedCenter += provider->deltaPrevImage;
        currPos = draggedCenter;
        currPos.clip(imW, imH);
        int newCenterX = int((double(currPos.x) - halfSizeW) / halfSizeW * 100.);
        int newCenterY = int((double(currPos.y) - halfSizeH) / halfSizeH * 100.);

        if (newCenterX != centerX->getIntValue() || newCenterY != centerY->getIntValue()) {
            centerX->setValue(newCenterX);
            centerY->setValue(newCenterY);
            updateGeometry (newCenterX, newCenterY, locY->getValue(), degree->getValue(), locX->getValue(), locYT->getValue(), locXL->getValue());

            if (listener) {
                listener->panelChanged (EvlocallabCenter, Glib::ustring::compose ("X=%1\nY=%2", centerX->getTextValue(), centerY->getTextValue()));
            }

            return true;
        }
    }

    return false;
}

void Locallab::switchOffEditMode ()
{
    if (edit->get_active()) {
        // switching off the toggle button
        bool wasBlocked = editConn.block(true);
        edit->set_active(false);

        if (!wasBlocked) {
            editConn.block(false);
        }
    }

    EditSubscriber::switchOffEditMode();  // disconnect
}

