# Python API Test
#
#
from vsp_wrapper import Geometry, GeometryType
import vsp_g
vsp=vsp_g
#import vsp

stdout = vsp.cvar.cstdout
errorMgr = vsp.ErrorMgrSingleton_getInstance()

#==== Use Case 1 ==== #
vsp.VSPCheckSetup();
errorMgr.PopErrorAndPrint( stdout )

types = vsp.GetGeomTypes()
errorMgr.PopErrorAndPrint( stdout )

print types

# Add Fuse
fuse = Geometry(GeometryType.FUSELAGE)
errorMgr.PopErrorAndPrint( stdout )

# Add Pod
pod = fuse.create_child(GeometryType.POD)
errorMgr.PopErrorAndPrint( stdout )

# Set Name
pod.name = "Pod"
errorMgr.PopErrorAndPrint( stdout )

# Change Length
pod.length = 7.0

# Change Finess Ratio
pod.fineratio = 10.0

# Change Y Location
pod.y_location = 1.0

# Change X Location
pod.x_location = 3.0

# Change Symmetry
pod.sym_planar_flag = vsp.SYM_XZ

# Copy Pod Geom
vsp.CopyGeomToClipboard( pod.id )
vsp.PasteGeomClipboard( fuse.id ) # make fuse parent

# Set Name
pod.name = "Original_Pod"
second_pod_id = vsp.FindGeom( "Pod", 0 )

# Change Location and Symmetry
vsp.SetParmVal( second_pod_id, "Sym_Planar_Flag", "Sym", 0 )
vsp.SetParmVal( second_pod_id, "Y_Location", "XForm", 0.0 )
vsp.SetParmVal( second_pod_id, "Z_Location", "XForm", 1.0 )

fname = "apitest.vsp3"

vsp.WriteVSPFile( fname )

geoms = vsp.FindGeoms()

print "All geoms in Vehicle."
print geoms

errorMgr.PopErrorAndPrint( stdout )

#==== Use Case 2 ====#

vsp.VSPRenew();
errorMgr.PopErrorAndPrint( stdout )

geoms = vsp.FindGeoms()

print "All geoms in Vehicle."
print geoms

# Add Fuse
fuse_id = vsp.AddGeom( "FUSELAGE" )

# Get XSec Surf ID
xsurf_id = vsp.GetXSecSurf( fuse_id, 0 )

# Change Type of First XSec
vsp.ChangeXSecType( xsurf_id, 0, vsp.SUPER_ELLIPSE )
errorMgr.PopErrorAndPrint( stdout )

# Change Type First XSec Properties
xsec_id = vsp.GetXSec( xsurf_id, 0 )
width_id = vsp.GetXSecParm( xsec_id, "Super_Width" )
height_id = vsp.GetXSecParm( xsec_id, "Super_Height" )
vsp.SetParmVal( width_id, 4.0 )
vsp.SetParmVal( height_id, 2.0 )

# Copy Cross-Section to Clipboard
vsp.CopyXSec( xsurf_id, 0)

# Paste Cross-Section
vsp.PasteXSec( xsurf_id, 1 )
vsp.PasteXSec( xsurf_id, 2 )
vsp.PasteXSec( xsurf_id, 3 )

# Change Type to File XSec

vsp.ChangeXSecType( xsurf_id, 0, vsp.FILE_FUSE )
file_xsec_id = vsp.GetXSec( xsurf_id, 0 )

# Build Point Vec

pnt_vec = vsp.Vec3dVec();
pnt_vec.push_back( vsp.vec3d( 0.0, 0.0, 2.0 ) )
pnt_vec.push_back( vsp.vec3d( 0.0, 1.0, 0.0 ) )
pnt_vec.push_back( vsp.vec3d( 0.0, 0.0,-2.0 ) )
pnt_vec.push_back( vsp.vec3d( 0.0,-1.0, 0.0 ) )
pnt_vec.push_back( vsp.vec3d( 0.0, 0.0, 2.0 ) )

# Load Points Into XSec
vsp.SetXSecPnts( file_xsec_id, pnt_vec )

geoms = vsp.FindGeoms()

print "End of second use case, all geoms in Vehicle."
print geoms

#==== Use Case 3 ====#

print "Start of third use case, read in first-case file."

#==== Read Geometry From File ====#
vsp.VSPRenew();
errorMgr.PopErrorAndPrint( stdout )

vsp.ReadVSPFile( fname )

geoms = vsp.FindGeoms()

print "All geoms in Vehicle."
print geoms

# Check for errors

num_err = errorMgr.GetNumTotalErrors()
for i in range(0,num_err):
	err = errorMgr.PopLastError()
	print "error = ", err.m_ErrorString

vsp.StartGui()

