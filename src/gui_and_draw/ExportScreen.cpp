//
// This file is released under the terms of the NASA Open Source Agreement (NOSA)
// version 1.3 as detailed in the LICENSE file which accompanies this software.
//
//
//////////////////////////////////////////////////////////////////////

#include "ExportScreen.h"
#include "ScreenMgr.h"
#include "EventMgr.h"
#include "Vehicle.h"
#include "StlHelper.h"

#include <assert.h>

//==== Constructor ====//
ExportScreen::ExportScreen( ScreenMgr* mgr ) : VspScreen( mgr )
{
    ExportFileUI* ui = m_ExportFileUI = new ExportFileUI();
    m_FLTK_Window = ui->UIWindow;
    m_SelectedSetIndex = 0;

    ui->setChoice->callback( staticScreenCB, this );

    ui->sterolithButton->callback( staticScreenCB, this );
    ui->xsecButton->callback( staticScreenCB, this );
    ui->rhinoButton->callback( staticScreenCB, this );
    ui->nascartButton->callback( staticScreenCB, this );
    ui->cart3dButton->callback( staticScreenCB, this );
    ui->povrayButton->callback( staticScreenCB, this );
    ui->gmshButton->callback( staticScreenCB, this );
    ui->x3dButton->callback( staticScreenCB, this );
    ui->stepButton->callback( staticScreenCB, this );
    ui->bezButton->callback( staticScreenCB, this );

}

//==== Update Screen ====//
bool ExportScreen::Update()
{
    LoadSetChoice();

    return true;
}

//==== Show Screen ====//
void ExportScreen::Show()
{
    Update();
    m_FLTK_Window->show();
}


//==== Hide Screen ====//
void ExportScreen::Hide()
{
    m_FLTK_Window->hide();
}

//==== Load Type Choice ====//
void ExportScreen::LoadSetChoice()
{
    m_ExportFileUI->setChoice->clear();

    Vehicle* veh = m_ScreenMgr->GetVehiclePtr();
    vector< string > set_name_vec = veh->GetSetNameVec();

    for ( int i = 0 ; i < ( int )set_name_vec.size() ; i++ )
    {
        m_ExportFileUI->setChoice->add( set_name_vec[i].c_str() );
    }

    m_ExportFileUI->setChoice->value( m_SelectedSetIndex );

}

void ExportScreen::ExportFile( string &newfile, int write_set, int type )
{
    Vehicle* veh = m_ScreenMgr->GetVehiclePtr();

    if ( type == ExportScreen::XSEC )
    {
        newfile = m_ScreenMgr->GetSelectFileScreen()->FileChooser( "Write XSec File?", "*.hrm" );
    }
    else if ( type == ExportScreen::STEREOLITH )
    {
        newfile = m_ScreenMgr->GetSelectFileScreen()->FileChooser( "Write STL File?", "*.stl" );
    }
    else if ( type == ExportScreen::RHINO3D )
    {
        newfile = m_ScreenMgr->GetSelectFileScreen()->FileChooser( "Write Rhino3D File?", "*.3dm" );
    }
    else if ( type == ExportScreen::NASCART )
    {
        newfile = m_ScreenMgr->GetSelectFileScreen()->FileChooser( "Write NASCART Files?", "*.dat" );
    }
    else if ( type == ExportScreen::CART3D )
    {
        newfile = m_ScreenMgr->GetSelectFileScreen()->FileChooser( "Write Cart3D Files?", "*.tri" );
    }
    else if ( type == ExportScreen::GMSH )
    {
        newfile = m_ScreenMgr->GetSelectFileScreen()->FileChooser( "Write GMsh Files?", "*.msh" );
    }
    else if ( type == ExportScreen::POVRAY )
    {
        newfile = m_ScreenMgr->GetSelectFileScreen()->FileChooser( "Write POVRAY File?", "*.pov" );
    }
    else if ( type == ExportScreen::X3D )
    {
        newfile = m_ScreenMgr->GetSelectFileScreen()->FileChooser( "Write X3D File?", "*.x3d" );
    }
    else if ( type == ExportScreen::STEP )
    {
        newfile = m_ScreenMgr->GetSelectFileScreen()->FileChooser( "Write STEP File?", "*.stp" );
    }
    else if ( type == ExportScreen::BEZ )
    {
        newfile = m_ScreenMgr->GetSelectFileScreen()->FileChooser( "Write Bezier File?", "*.bez" );
    }
    else if ( type == -1 )
    {
        m_ExportFileUI->UIWindow->show();
        return;
    }

    if ( newfile.size() != 0 && newfile[ newfile.size() - 1] != '/' )
//jrg back() only in c++11 and remove strcmp
//  if ( newfile.compare("") != 0 && strcmp( &newfile.back(), "/") != 0 )
    {
        if ( type == ExportScreen::XSEC )
        {
            veh->WriteXSecFile( newfile, write_set );
        }
        else if ( type == ExportScreen::STEREOLITH )
        {
            veh->WriteSTLFile( newfile, write_set );
        }
        else if ( type == ExportScreen::CART3D )
        {
            veh->WriteTRIFile( newfile, write_set );
        }
        else if ( type == ExportScreen::NASCART )
        {
            veh->WriteNascartFiles( newfile, write_set );
        }
        else if ( type == ExportScreen::GMSH )
        {
            veh->WriteGmshFile( newfile, write_set );
        }
        else if ( type == ExportScreen::POVRAY )
        {
            veh->WritePovRayFile( newfile, write_set );
        }
        else if ( type == ExportScreen::X3D )
        {
            veh->WriteX3DFile( newfile, write_set );
        }
        else if ( type == ExportScreen::STEP )
        {
            veh->WriteSTEPFile( newfile, write_set );
        }
        else if ( type == ExportScreen::BEZ )
        {
            veh->WriteBezFile( newfile, write_set );
        }
    }
    m_ExportFileUI->UIWindow->hide();

}


//==== Callbacks ====//
void ExportScreen::CallBack( Fl_Widget *w )
{
    Vehicle* veh = m_ScreenMgr->GetVehiclePtr();
    string newfile;

    if ( w ==   m_ExportFileUI->xsecButton )        // Export CrossSection File
    {
        ExportFile( newfile, m_SelectedSetIndex, ExportScreen::XSEC );
    }
    else if ( w == m_ExportFileUI->sterolithButton )
    {
        ExportFile( newfile, m_SelectedSetIndex, ExportScreen::STEREOLITH );
    }
    else if ( w == m_ExportFileUI->nascartButton )
    {
        ExportFile( newfile, m_SelectedSetIndex, ExportScreen::NASCART );
    }
    else if ( w == m_ExportFileUI->cart3dButton )   // Export Tri File
    {
        ExportFile( newfile, m_SelectedSetIndex, ExportScreen::CART3D );
    }
    else if ( w == m_ExportFileUI->gmshButton )
    {
        ExportFile( newfile, m_SelectedSetIndex, ExportScreen::GMSH );
    }
    else if ( w == m_ExportFileUI->povrayButton )
    {
        ExportFile( newfile, m_SelectedSetIndex, ExportScreen::POVRAY );
    }
    else if ( w == m_ExportFileUI->x3dButton )
    {
        ExportFile( newfile, m_SelectedSetIndex, ExportScreen::X3D );
    }
    else if ( w == m_ExportFileUI->stepButton )
    {
        ExportFile( newfile, m_SelectedSetIndex, ExportScreen::STEP );
    }
    else if ( w == m_ExportFileUI->bezButton )
    {
        ExportFile( newfile, m_SelectedSetIndex, ExportScreen::BEZ );
    }
    else if ( w == m_ExportFileUI->setChoice )
    {
        m_SelectedSetIndex = m_ExportFileUI->setChoice->value();
    }

    m_ScreenMgr->SetUpdateFlag( true );
//  m_ScreenMgr->UpdateAllScreens();
}


