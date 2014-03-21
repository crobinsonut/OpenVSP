# Python API Test
#
#
from vsp_wrapper import Geometry, GeometryType
import vsp_g as vsp

stdout = vsp.cvar.cstdout
errorMgr = vsp.ErrorMgrSingleton_getInstance()

#==== Use Case 1 ==== #
vsp.VSPCheckSetup();
errorMgr.PopErrorAndPrint( stdout )

types = vsp.GetGeomTypes()
errorMgr.PopErrorAndPrint( stdout )

print types

# Add Fuse
fuse = Geometry.add_geometry(GeometryType.FUSELAGE)
errorMgr.PopErrorAndPrint( stdout )

# Add Pod
pod = Geometry.add_geometry(GeometryType.POD, fuse)
errorMgr.PopErrorAndPrint( stdout )

# Set Name
pod.name = "Pod"
errorMgr.PopErrorAndPrint( stdout )

# Change Length
pod.Length = 7.0

# Change Finess Ratio
pod.FineRatio = 10.0

# Change Y Location
pod.Y_Location = 1.0

# Change X Location
pod.X_Location = 3.0

# Change Symmetry
pod.Sym_Planar_Flag = vsp.SYM_XZ

# Copy Pod Geom
vsp.CopyGeomToClipboard(pod.id)
vsp.PasteGeomClipboard(fuse.id) #make fuse parent

pod.name = "Original_Pod"
second_pod =  Geometry(vsp.FindGeom("Pod", 0)) 

# Set Name
pod.name = "Original_Pod"

# Change Location and Symmetry
second_pod.Sym_Planar_Flag = 0
second_pod.Y_Location = 0.0
second_pod.Z_Location = 1.0

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
fuse = Geometry.add_geometry(GeometryType.FUSELAGE)

# Get XSec Surf ID
xsurf = fuse.xsec_surfs[0]

# Change Type of First XSec
vsp.ChangeXSecType( xsurf.id, 0, vsp.SUPER_ELLIPSE )
errorMgr.PopErrorAndPrint( stdout )

# Change Type First XSec Properties
xsec = xsurf.xsecs[0]
xsec.Super_Width = 4.0
xsec.Super_Height = 2.0

# Copy Cross-Section to Clipboard
vsp.CopyXSec( xsurf.id, 0)

# Paste Cross-Section
#vsp.PasteXSec( xsurf.id, 1 )
#vsp.PasteXSec( xsurf.id, 2 )
#vsp.PasteXSec( xsurf.id, 3 )

# Change Type to File XSec

vsp.ChangeXSecType( xsurf.id, 0, vsp.FILE_FUSE )
file_xsec = xsurf.xsecs[0]

# Build Point Vec
pnt_vec = vsp.Vec3dVec();
pnt_vec.push_back( vsp.vec3d( 0.0, 0.0, 2.0 ) )
pnt_vec.push_back( vsp.vec3d( 0.0, 1.0, 0.0 ) )
pnt_vec.push_back( vsp.vec3d( 0.0, 0.0,-2.0 ) )
pnt_vec.push_back( vsp.vec3d( 0.0,-1.0, 0.0 ) )
pnt_vec.push_back( vsp.vec3d( 0.0, 0.0, 2.0 ) )

# Load Points Into XSec
file_xsec.set_points(pnt_vec)

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

