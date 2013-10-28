
#include "xsec_surf.h"
#include "geom.h"
#include <cmath>

//==== Get area normal vector for PLANAR cross section ====//
vec3d Xsec_surf::get_area_normal( int ixs )
{
	vec3d areaNormal, refVec1, refVec2;
	int halfPnts = (num_pnts + 1) / 2;

	refVec1 = pnts_xsecs(ixs, halfPnts - 1) - pnts_xsecs(ixs, 0);
	refVec2 = pnts_xsecs(ixs, 1) - pnts_xsecs(ixs, num_pnts - 2);

	areaNormal = cross(refVec2, refVec1);
	areaNormal.normalize();

	return areaNormal;
}

//==== Get reflected area normal vector for PLANAR cross section ====//
vec3d Xsec_surf::get_refl_area_normal( int ixs )
{
	vec3d areaNormal, refVec1, refVec2;
	int halfPnts = (num_pnts + 1) / 2;

	refVec1 = refl_pnts_xsecs(ixs, halfPnts - 1) - refl_pnts_xsecs(ixs, 0);
	refVec2 = refl_pnts_xsecs(ixs, 1) - refl_pnts_xsecs(ixs, num_pnts - 2);

	areaNormal = cross(refVec2, refVec1);
	areaNormal.normalize();

	return areaNormal;
}

//==== Get Cross Section Area ====//
// JBB: This only works for PLANAR cross sections.
double Xsec_surf::get_xsec_area( int ixs )
{
	double area = 0;
	int j;
	vec3d xAxis(1,0,0), yAxis(0,1,0), zAxis(0,0,1), areaNormal = get_area_normal(ixs);
	vector<vec3d> pnts;

	// Load cross section points
	for( int i = 0; i < num_pnts; i++ )
	{
		pnts.push_back( pnts_xsecs(ixs,i) );
	}

	// Get rotation in xy plane to align areaNormal with yaxis
	vec2d an2(areaNormal.x(), areaNormal.y()), ax1(1,0), ax2(0,1);
	double theta = angle(an2, ax2);
	if (areaNormal.x() < 0) theta = -theta;
	// Rotate areaNormal around z to y axis
	areaNormal   = RotateArbAxis( areaNormal, theta, zAxis );
	// Rotate points as well
	for ( int i = 0; i < num_pnts; i++ )
		pnts[i] = RotateArbAxis( pnts[i], theta, zAxis );

	// Get rotation in yz plane to align areaNormal with yaxis
	an2.set_xy(areaNormal.y(), areaNormal.z());
	theta = angle(an2, ax1);
	if(areaNormal.z() >= 0) theta = -theta;
	// Rotate points
	for ( int i = 0; i < num_pnts; i++ )
		pnts[i] = RotateArbAxis( pnts[i], theta, xAxis );

	// Cross section should now be in x-z plane. Calculate area.
	for( j = 1; j < num_pnts; j++ )
		area += pnts[j-1].x()*pnts[j].z() - pnts[j].x()*pnts[j-1].z();
	area += pnts[j-1].x()*pnts[0].z() - pnts[0].x()*pnts[j-1].z();

	return (0.5*area);
}

//==== Cross sectional area in a given plane. ====//
double Xsec_surf::get_xsec_plane_area( int ixs, int plane, float mat[4][4] )
{
	double area = 0;
	int j;
	vector<vec3d> pnts;

	for ( int i = 0; i < num_pnts; i++ )
	{
		pnts.push_back( pnts_xsecs(ixs,i).transform(mat) );
	}

	switch (plane)
	{
	case XY_PLANE:
		for( j = 1; j < num_pnts; j++ )
			area += pnts[j-1].x()*pnts[j].y() - pnts[j].x()*pnts[j-1].y();
		area += pnts[j-1].x()*pnts[0].y() - pnts[0].x()*pnts[j-1].y();
		break;
	case XZ_PLANE:
		for( j = 1; j < num_pnts; j++ )
			area += pnts[j-1].x()*pnts[j].z() - pnts[j].x()*pnts[j-1].z();
		area += pnts[j-1].x()*pnts[0].z() - pnts[0].x()*pnts[j-1].z();
		break;
	case YZ_PLANE:
		for( j = 1; j < num_pnts; j++ )
			area += pnts[j-1].y()*pnts[j].z() - pnts[j].y()*pnts[j-1].z();
		area += pnts[j-1].y()*pnts[0].z() - pnts[0].y()*pnts[j-1].z();
		break;
	}

	return (0.5 * area);
}

//==== Cross sectional area in a given plane. ====//
double Xsec_surf::get_refl_xsec_plane_area( int ixs, int plane, float refl_mat[4][4] )
{
	double area = 0;
	int j;
	vector<vec3d> pnts;

	for ( int i = 0; i < num_pnts; i++ )
	{
		pnts.push_back( refl_pnts_xsecs(ixs,i).transform(refl_mat) );
	}

	switch (plane)
	{
	case XY_PLANE:
		for( j = 1; j < num_pnts; j++ )
			area += pnts[j-1].x()*pnts[j].y() - pnts[j].x()*pnts[j-1].y();
		area += pnts[j-1].x()*pnts[0].y() - pnts[0].x()*pnts[j-1].y();
		break;
	case XZ_PLANE:
		for( j = 1; j < num_pnts; j++ )
			area += pnts[j-1].x()*pnts[j].z() - pnts[j].x()*pnts[j-1].z();
		area += pnts[j-1].x()*pnts[0].z() - pnts[0].x()*pnts[j-1].z();
		break;
	case YZ_PLANE:
		for( j = 1; j < num_pnts; j++ )
			area += pnts[j-1].y()*pnts[j].z() - pnts[j].y()*pnts[j-1].z();
		area += pnts[j-1].y()*pnts[0].z() - pnts[0].y()*pnts[j-1].z();
		break;
	}

	return (0.5 * area);
}

//==== Get cross section centroid location. Only works for PLANAR cross sections! ====//
vec3d Xsec_surf::get_xsec_centroid( int ixs )
{
	double xcg = 0, zcg = 0, area = get_xsec_area(ixs);
	int j;
	vec3d xAxis(1,0,0), zAxis(0,0,1), areaNormal = get_area_normal(ixs);

	vector<vec3d> pnts;
	// Load cross section points
	for( int i = 0; i < num_pnts; i++ )
	{
		pnts.push_back( pnts_xsecs(ixs,i) );
	}

	// Get rotation in xy plane to align areaNormal with yaxis
	vec2d an2(areaNormal.x(), areaNormal.y()), ax1(1,0), ax2(0,1);
	double theta2, theta1 = angle(an2, ax2);
	if (areaNormal.x() < 0) theta1 = -theta1;
	// Rotate areaNormal around z to y axis
	areaNormal   = RotateArbAxis( areaNormal, theta1, zAxis );
	// Rotate points as well
	for ( int i = 0; i < num_pnts; i++ )
		pnts[i] = RotateArbAxis( pnts[i], theta1, zAxis );

	// Get rotation in yz plane to align areaNormal with yaxis
	an2.set_xy(areaNormal.y(), areaNormal.z());
	theta2 = angle(an2, ax1);
	if(areaNormal.z() >= 0) theta2 = -theta2;
	// Rotate points about x axis
	for ( int i = 0; i < num_pnts; i++ )
		pnts[i] = RotateArbAxis( pnts[i], theta2, xAxis );

	// Cross section should now be in x-z plane. Calculate centroid.
	for ( j = 1; j < num_pnts; j++ )
	{
		xcg += (pnts[j-1].x()*pnts[j].z() - pnts[j].x()*pnts[j-1].z()) \
				* ( pnts[j-1].x() + pnts[j].x() );
		zcg += (pnts[j-1].x()*pnts[j].z() - pnts[j].x()*pnts[j-1].z()) \
				* ( pnts[j-1].z() + pnts[j].z() );
	}
	xcg += (pnts[j-1].x()*pnts[0].z() - pnts[0].x()*pnts[j-1].z()) \
			* ( pnts[j-1].x() + pnts[0].x() );
	zcg += (pnts[j-1].x()*pnts[0].z() - pnts[0].x()*pnts[j-1].z()) \
			* ( pnts[j-1].z() + pnts[0].z() );
	xcg /= (6*area);
	zcg /= (6*area);

	// Xcg vector including y position (same for all points since xsec rotated into xz plane)
	vec3d cgLoc(xcg, pnts[0].y(), zcg);
	// Rotate back into original coordinate frame
	cgLoc = RotateArbAxis( cgLoc, -theta2, xAxis );
	cgLoc = RotateArbAxis( cgLoc, -theta1, zAxis );

	return cgLoc;
}

//==== Get cross section centroid location. Only works for PLANAR cross sections! ====//
vec3d Xsec_surf::get_refl_xsec_centroid( int ixs )
{
	double xcg = 0, zcg = 0, area = get_xsec_area(ixs);
	int j;
	vec3d xAxis(1,0,0), zAxis(0,0,1), areaNormal = get_refl_area_normal(ixs);

	vector<vec3d> pnts;
	// Load cross section points
	for( int i = 0; i < num_pnts; i++ )
	{
		pnts.push_back( refl_pnts_xsecs(ixs,i) );
	}

	// Get rotation in xy plane to align areaNormal with yaxis
	vec2d an2(areaNormal.x(), areaNormal.y()), ax1(1,0), ax2(0,1);
	double theta2, theta1 = angle(an2, ax2);
	if (areaNormal.x() < 0) theta1 = -theta1;
	// Rotate areaNormal around z to y axis
	areaNormal   = RotateArbAxis( areaNormal, theta1, zAxis );
	// Rotate points as well
	for ( int i = 0; i < num_pnts; i++ )
		pnts[i] = RotateArbAxis( pnts[i], theta1, zAxis );

	// Get rotation in yz plane to align areaNormal with yaxis
	an2.set_xy(areaNormal.y(), areaNormal.z());
	theta2 = angle(an2, ax1);
	if(areaNormal.z() >= 0) theta2 = -theta2;
	// Rotate points about x axis
	for ( int i = 0; i < num_pnts; i++ )
		pnts[i] = RotateArbAxis( pnts[i], theta2, xAxis );

	// Cross section should now be in x-z plane. Calculate centroid.
	for ( j = 1; j < num_pnts; j++ )
	{
		xcg += (pnts[j-1].x()*pnts[j].z() - pnts[j].x()*pnts[j-1].z()) \
				* ( pnts[j-1].x() + pnts[j].x() );
		zcg += (pnts[j-1].x()*pnts[j].z() - pnts[j].x()*pnts[j-1].z()) \
				* ( pnts[j-1].z() + pnts[j].z() );
	}
	xcg += (pnts[j-1].x()*pnts[0].z() - pnts[0].x()*pnts[j-1].z()) \
			* ( pnts[j-1].x() + pnts[0].x() );
	zcg += (pnts[j-1].x()*pnts[0].z() - pnts[0].x()*pnts[j-1].z()) \
			* ( pnts[j-1].z() + pnts[0].z() );
	xcg /= (6*area);
	zcg /= (6*area);

	// Xcg vector including y position (same for all points since xsec rotated into xz plane)
	vec3d cgLoc(xcg, pnts[0].y(), zcg);
	// Rotate back into original coordinate frame
	cgLoc = RotateArbAxis( cgLoc, -theta2, xAxis );
	cgLoc = RotateArbAxis( cgLoc, -theta1, zAxis );

	return cgLoc;
}

vec2d Xsec_surf::get_xsec_centroid_in_plane(int ixs, int plane, float mat[4][4])
{
	double xcg = 0, ycg = 0, area = get_xsec_plane_area(ixs, plane, mat);
	int j;
	vector<vec2d> pnts;
	vec2d tmpPnts;

	switch(plane)
	{
	case XY_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( pnts_xsecs(ixs,i).transform(mat).x(), pnts_xsecs(ixs,i).transform(mat).y() );
			pnts.push_back(tmpPnts);
		}
		break;
	case XZ_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( pnts_xsecs(ixs,i).transform(mat).x(), pnts_xsecs(ixs,i).transform(mat).z() );
			pnts.push_back(tmpPnts);
		}
		break;
	case YZ_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( pnts_xsecs(ixs,i).transform(mat).y(), pnts_xsecs(ixs,i).transform(mat).z() );
			pnts.push_back(tmpPnts);
		}
		break;
	}

	for ( j = 1; j < num_pnts; j++ )
	{
		xcg += (pnts[j-1].x()*pnts[j].y() - pnts[j].x()*pnts[j-1].y()) \
				* ( pnts[j-1].x() + pnts[j].x() );
		ycg += (pnts[j-1].x()*pnts[j].y() - pnts[j].x()*pnts[j-1].y()) \
				* ( pnts[j-1].y() + pnts[j].y() );
	}
	xcg += (pnts[j-1].x()*pnts[0].y() - pnts[0].x()*pnts[j-1].y()) \
			* ( pnts[j-1].x() + pnts[0].x() );
	ycg += (pnts[j-1].x()*pnts[0].y() - pnts[0].x()*pnts[j-1].y()) \
			* ( pnts[j-1].y() + pnts[0].y() );
	xcg /= (6*area);
	ycg /= (6*area);

	vec2d cgLoc( xcg, ycg );
	return cgLoc;
}

vec2d Xsec_surf::get_refl_xsec_centroid_in_plane(int ixs, int plane, float refl_mat[4][4])
{
	double xcg = 0, ycg = 0, area = get_refl_xsec_plane_area(ixs, plane, refl_mat);
	int j;
	vector<vec2d> pnts;
	vec2d tmpPnts;

	switch(plane)
	{
	case XY_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( refl_pnts_xsecs(ixs,i).transform(refl_mat).x(), refl_pnts_xsecs(ixs,i).transform(refl_mat).y() );
			pnts.push_back(tmpPnts);
		}
		break;
	case XZ_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( refl_pnts_xsecs(ixs,i).transform(refl_mat).x(), refl_pnts_xsecs(ixs,i).transform(refl_mat).z() );
			pnts.push_back(tmpPnts);
		}
		break;
	case YZ_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( refl_pnts_xsecs(ixs,i).transform(refl_mat).y(), refl_pnts_xsecs(ixs,i).transform(refl_mat).z() );
			pnts.push_back(tmpPnts);
		}
		break;
	}

	for ( j = 1; j < num_pnts; j++ )
	{
		xcg += (pnts[j-1].x()*pnts[j].y() - pnts[j].x()*pnts[j-1].y()) \
				* ( pnts[j-1].x() + pnts[j].x() );
		ycg += (pnts[j-1].x()*pnts[j].y() - pnts[j].x()*pnts[j-1].y()) \
				* ( pnts[j-1].y() + pnts[j].y() );
	}
	xcg += (pnts[j-1].x()*pnts[0].y() - pnts[0].x()*pnts[j-1].y()) \
			* ( pnts[j-1].x() + pnts[0].x() );
	ycg += (pnts[j-1].x()*pnts[0].y() - pnts[0].x()*pnts[j-1].y()) \
			* ( pnts[j-1].y() + pnts[0].y() );
	xcg /= (6*area);
	ycg /= (6*area);

	vec2d cgLoc( xcg, ycg );
	return cgLoc;
}

DegenGeom* Xsec_surf::createSurfDegenGeom(Geom* parentGeom, int sym_code_in, float mat[4][4], float refl_mat[4][4])
{
	DegenGeom*	degenGeom = new DegenGeom();
	degenGeom->setParentGeom(parentGeom);
	degenGeom->setType(DegenGeom::SURFACE_TYPE);

	degenGeom->setNumXSecs( num_xsecs );
	degenGeom->setNumPnts( num_pnts );

	createDegenSurface(degenGeom, sym_code_in, mat);
	if(sym_code_in != NO_SYM)
		createDegenSurface_refl(degenGeom, sym_code_in, refl_mat);

	createSurfDegenPlate(degenGeom, sym_code_in, mat);
	if(sym_code_in != NO_SYM)
		createSurfDegenPlate_refl(degenGeom, sym_code_in, refl_mat);

	createSurfDegenStick(degenGeom, sym_code_in, mat, refl_mat);

	return degenGeom;
}

DegenGeom* Xsec_surf::createBodyDegenGeom(Geom* parentGeom, int sym_code_in, float mat[4][4], float refl_mat[4][4])
{
	DegenGeom*	degenGeom = new DegenGeom();
	degenGeom->setParentGeom(parentGeom);
	degenGeom->setType(DegenGeom::BODY_TYPE);

	degenGeom->setNumXSecs( num_xsecs );
	degenGeom->setNumPnts( num_pnts );

	createDegenSurface(degenGeom, sym_code_in, mat);
	if(sym_code_in != NO_SYM)
		createDegenSurface_refl(degenGeom, sym_code_in, refl_mat);

	createBodyDegenPlate(degenGeom, sym_code_in, mat, refl_mat);
	createBodyDegenStick(degenGeom, sym_code_in, mat, refl_mat);

	return degenGeom;
}

void Xsec_surf::createDegenSurface(DegenGeom* degenGeom, int sym_code_in, float mat[4][4])
{
	DegenSurface	degenSurface = degenGeom->getDegenSurface();

	int nLow = 0, nHigh = num_xsecs;

	vector< vector<vec3d> > xSurfMat;
	vector< vector<vec3d> > nSurfMat;

	vector<vec3d> xVec( num_pnts );
	vector<vec3d> nVec( num_pnts );

	if ( degenGeom->getType() == DegenGeom::SURFACE_TYPE )
	{
		if ( degenGeom->getParentGeom()->getTypeStr() == "wing" || degenGeom->getParentGeom()->getTypeStr() == "prop" )
		{
			nLow  = 1;
			nHigh = num_xsecs - 1;
		}
	}

	for ( int i = nLow; i < nHigh-1; i++ )
	{
		nVec.assign(num_pnts,vec3d(NAN,NAN,NAN));

		for ( int j = 0; j < num_pnts-1; j++ )
		{
			vec3d sVec1 = pnts_xsecs(i+1,j).transform(mat) - pnts_xsecs(i,j).transform(mat);
			vec3d sVec2 = pnts_xsecs(i,j+1).transform(mat) - pnts_xsecs(i,j).transform(mat);
			nVec[j]     = cross(sVec1,sVec2);
			nVec[j].normalize();
		}

		nSurfMat.push_back(nVec);
	}
	nVec.assign(num_pnts,vec3d(NAN,NAN,NAN));
	nSurfMat.push_back(nVec);



	for ( int i = nLow; i < nHigh; i++ )
	{
		for ( int j = 0; j < num_pnts; j++ )
		{
			xVec[j] = pnts_xsecs(i,j).transform(mat);
		}

		degenSurface.u.push_back( uArray[i] );

		xSurfMat.push_back(xVec);
	}


	for ( int j = 0; j < num_pnts; j++ )
	{
		degenSurface.w.push_back( wArray[j] );
	}

	degenSurface.x    = xSurfMat;
	degenSurface.nvec = nSurfMat;
	degenGeom->setDegenSurface(degenSurface);

	return;
}

void Xsec_surf::createDegenSurface_refl(DegenGeom* degenGeom, int sym_code_in, float refl_mat[4][4])
{
	DegenSurface	degenSurface = degenGeom->getDegenSurface();

	int nLow = 0, nHigh = num_xsecs;

	vector< vector<vec3d> > xSurfMat = degenSurface.x;
	vector< vector<vec3d> > nSurfMat = degenSurface.nvec;

	vector<vec3d> xVec( num_pnts );
	vector<vec3d> nVec( num_pnts );

	if ( degenGeom->getType() == DegenGeom::SURFACE_TYPE )
	{
		if ( degenGeom->getParentGeom()->getTypeStr() == "wing" || degenGeom->getParentGeom()->getTypeStr() == "prop" )
		{
			nLow  = 1;
			nHigh = num_xsecs - 1;
		}
	}

	if ( sym_code_in != refl_pnts_xsecs_code )
	{
		refl_pnts_xsecs_code = sym_code_in;
		load_refl_pnts_xsecs();
	}

	for ( int i = nLow; i < nHigh-1; i++ )
	{
		nVec.assign(num_pnts,vec3d(NAN,NAN,NAN));

		for ( int j = 0; j < num_pnts-1; j++ )
		{
			vec3d sVec2 = refl_pnts_xsecs(i+1,j).transform(refl_mat) - refl_pnts_xsecs(i,j).transform(refl_mat);
			vec3d sVec1 = refl_pnts_xsecs(i,j+1).transform(refl_mat) - refl_pnts_xsecs(i,j).transform(refl_mat);
			nVec[j]     = cross(sVec1,sVec2);
			nVec[j].normalize();
		}

		nSurfMat.push_back(nVec);
	}
	nVec.assign(num_pnts,vec3d(NAN,NAN,NAN));
	nSurfMat.push_back(nVec);

	for ( int i = nLow; i < nHigh; i++ )
	{
		for ( int j = 0; j < num_pnts; j++ )
		{
			xVec[j] = refl_pnts_xsecs(i,j).transform(refl_mat);
		}

		degenSurface.u.push_back( uArray[i] );

		xSurfMat.push_back(xVec);
	}


	for ( int j = 0; j < num_pnts; j++ )
	{
		degenSurface.w.push_back( wArray[j] );
	}

	degenSurface.x    = xSurfMat;
	degenSurface.nvec = nSurfMat;
	degenGeom->setDegenSurface(degenSurface);
}


void Xsec_surf::createSurfDegenPlate(DegenGeom* degenGeom, int sym_code_in, float mat[4][4])
{
	DegenPlate	degenPlate = degenGeom->getDegenPlate();

	int nLow = 0, nHigh = num_xsecs;
	int platePnts = (num_pnts + 1) / 2;

	if ( degenGeom->getParentGeom()->getTypeStr() == "wing" || degenGeom->getParentGeom()->getTypeStr() == "prop" )
	{
		// Keep only airfoil sections, discard endcap close-out lines
		nLow  = 1;
		nHigh = num_xsecs - 1;
	}

	vector< vector<vec3d>  >	xMat;
	vector< vector<vec3d>  >	nCamberMat;
	vector< vector<double> >	tMat;
	vector< vector<double> >	zMat;

	vector<vec3d>	xVec(  platePnts );
	vector<vec3d>	nCVec( platePnts );
	vector<vec3d>	nPVec( platePnts );
	vector<double>	tVec(  platePnts );
	vector<double>	zVec(  platePnts );

	vec3d  topPnt, botPnt, chordVec, camberPnt, chordPnt, nPlate;

	for ( int i = nLow; i < nHigh; i++ )
	{
		// Set first point (trailing edge)
		xVec[0]  = pnts_xsecs(i,0).transform(mat);
		nCVec[0] = vec3d(0,0,0); // on camber line
		tVec[0]  = 0;
		zVec[0]  = 0;

		// Set last point (leading edge)
		xVec[platePnts-1]  = pnts_xsecs( i, platePnts-1 ).transform(mat);
		nCVec[platePnts-1] = vec3d(0,0,0); // on camber line
		tVec[platePnts-1]  = 0;
		zVec[platePnts-1]  = 0;

		//== Compute plate normal ==//
		// rotated chord vector
		vec3d rcv = pnts_xsecs(i, platePnts-1).transform(mat) - pnts_xsecs(i,0).transform(mat);
		// rotated area normal vector
		vec3d anv = get_area_normal(i).transform(mat) - vec3d(0,0,0).transform(mat);
		// plate normal vector
		nPlate = cross(rcv,anv);
		nPlate.normalize();
		degenPlate.nPlate.push_back( nPlate );

		// normalized, unrotated chord vector (te->le)
		chordVec = pnts_xsecs(i, platePnts-1) - pnts_xsecs(i,0);
		chordVec.normalize();

		// Set points along af section
		for ( int j = 1, k = num_pnts-2; j < platePnts-1; j++, k-- )
		{
			topPnt = pnts_xsecs(i,j);
			botPnt = pnts_xsecs(i,k);

			camberPnt = (topPnt + botPnt) / 2;
			chordPnt  = pnts_xsecs(i,0) + chordVec * dot(camberPnt-pnts_xsecs(i,0),chordVec);

			xVec[j]  = chordPnt.transform(mat);
			nCVec[j] = topPnt.transform(mat)-botPnt.transform(mat);
			nCVec[j].normalize();
			tVec[j]  = dist(topPnt,botPnt);
			zVec[j]  = dist(camberPnt,chordPnt);
			if ( angle(camberPnt - chordPnt, nPlate) != 0) zVec[j] *= -1;
		}

		xMat.push_back(xVec);
		nCamberMat.push_back(nCVec);
		tMat.push_back(tVec);
		zMat.push_back(zVec);

		degenPlate.u.push_back( uArray[i] );
	}

	for ( int j = 0, k = num_pnts-1; j < platePnts-1; j++, k-- )
	{
		degenPlate.wTop.push_back( wArray[j] );
		degenPlate.wBot.push_back( wArray[k] );
	}
	degenPlate.wTop.push_back( wArray[platePnts-1] );
	degenPlate.wBot.push_back( wArray[platePnts-1] );

	degenPlate.x    	= xMat;
	degenPlate.nCamber 	= nCamberMat;
	degenPlate.t		= tMat;
	degenPlate.zcamber	= zMat;
	degenGeom->setDegenPlate(degenPlate);

	return;
}


void Xsec_surf::createSurfDegenPlate_refl(DegenGeom* degenGeom, int sym_code_in, float refl_mat[4][4])
{
	DegenPlate	degenPlate = degenGeom->getDegenPlate();

	int nLow = 0, nHigh = num_xsecs;
	int platePnts = (num_pnts + 1) / 2;

	if ( degenGeom->getParentGeom()->getTypeStr() == "wing" || degenGeom->getParentGeom()->getTypeStr() == "prop" )
	{
		// Keep only airfoil sections, discard endcap close-out lines
		nLow  = 1;
		nHigh = num_xsecs - 1;
	}

	vector< vector<vec3d>  >	xMat = degenPlate.x;
	vector< vector<vec3d>  >	nCamberMat = degenPlate.nCamber;
	vector< vector<double> >	tMat = degenPlate.t;
	vector< vector<double> >	zMat = degenPlate.zcamber;

	vector<vec3d>	xVec(  platePnts );
	vector<vec3d>	nCVec( platePnts );
	vector<vec3d>	nPVec( platePnts );
	vector<double>	tVec(  platePnts );
	vector<double>	zVec(  platePnts );

	vec3d  topPnt, botPnt, chordVec, camberPnt, chordPnt, nPlate;


	if ( sym_code_in != refl_pnts_xsecs_code )
	{
		refl_pnts_xsecs_code = sym_code_in;
		load_refl_pnts_xsecs();
	}

	for ( int i = nLow; i < nHigh; i++ )
	{
		// Set first point (trailing edge)
		xVec[0]  = refl_pnts_xsecs(i,0).transform(refl_mat);
		nCVec[0] = vec3d(0,0,0); // on camber line
		tVec[0]  = 0;
		zVec[0]  = 0;

		// Set last point (leading edge)
		xVec[platePnts-1]  = refl_pnts_xsecs( i, platePnts-1 ).transform(refl_mat);
		nCVec[platePnts-1] = vec3d(0,0,0); // on camber line
		tVec[platePnts-1]  = 0;
		zVec[platePnts-1]  = 0;

		//== Compute plate normal ==//
		// rotated chord vector
		vec3d rcv = refl_pnts_xsecs(i, platePnts-1).transform(refl_mat) \
				  - refl_pnts_xsecs(i,0).transform(refl_mat);
		// rotated area normal vector
		vec3d anv = get_refl_area_normal(i).transform(refl_mat) - vec3d(0,0,0).transform(refl_mat);
		// plate normal vector
		nPlate = cross(rcv,anv);
		nPlate.normalize();
		degenPlate.nPlate.push_back( nPlate );

		// normalized, unrotated chord vector (te->le)
		chordVec = refl_pnts_xsecs(i, platePnts-1) - refl_pnts_xsecs(i,0);
		chordVec.normalize();

		// Set points along af section
		for ( int j = 1, k = num_pnts-2; j < platePnts-1; j++, k-- )
		{
			topPnt = refl_pnts_xsecs(i,j);
			botPnt = refl_pnts_xsecs(i,k);

			camberPnt = (topPnt + botPnt) / 2;
			chordPnt  = refl_pnts_xsecs(i,0) + chordVec * dot(camberPnt-refl_pnts_xsecs(i,0),chordVec);

			xVec[j]  = chordPnt.transform(refl_mat);
			nCVec[j] = topPnt.transform(refl_mat)-botPnt.transform(refl_mat);
			nCVec[j].normalize();
			tVec[j]  = dist(topPnt,botPnt);
			zVec[j]  = dist(camberPnt,chordPnt);
			if ( angle(camberPnt - chordPnt, nPlate) != 0) zVec[j] *= -1;
		}

		xMat.push_back(xVec);
		nCamberMat.push_back(nCVec);
		tMat.push_back(tVec);
		zMat.push_back(zVec);

		degenPlate.u.push_back( uArray[i] );
	}

	for ( int j = 0, k = num_pnts-1; j < platePnts-1; j++, k-- )
	{
		degenPlate.wTop.push_back( wArray[j] );
		degenPlate.wBot.push_back( wArray[k] );
	}
	degenPlate.wTop.push_back( wArray[platePnts-1] );
	degenPlate.wBot.push_back( wArray[platePnts-1] );

	degenPlate.x    	= xMat;
	degenPlate.nCamber 	= nCamberMat;
	degenPlate.t		= tMat;
	degenPlate.zcamber	= zMat;
	degenGeom->setDegenPlate(degenPlate);
}

void Xsec_surf::createBodyDegenPlate(DegenGeom* degenGeom, int sym_code_in, float mat[4][4], float refl_mat[4][4])
{
	DegenPlate	degenPlate = degenGeom->getDegenPlate();

	int nLow = 0, nHigh = num_xsecs;
	int platePnts = (num_pnts + 1) / 2;

	vector< vector<vec3d>  >	xMat;
	vector< vector<vec3d>  >	nCamberMat;
	vector< vector<double> >	tMat;
	vector< vector<double> >	zMat;

	vector<vec3d>	xVec(  platePnts );
	vector<vec3d>	nCVec( platePnts );
	vector<vec3d>	nPVec( platePnts );
	vector<double>	tVec(  platePnts );
	vector<double>	zVec(  platePnts );

	vec3d  topPnt, botPnt, chordVec, camberPnt, chordPnt, nPlate;

	// "Vertical" plate
	for ( int i = nLow; i < nHigh; i++ )
	{
		// Set first point (trailing edge)
		xVec[0]  = pnts_xsecs(i,0).transform(mat);
		nCVec[0] = vec3d(0,0,0); // on camber line
		tVec[0]  = 0;
		zVec[0]  = 0;

		// Set last point (leading edge)
		xVec[platePnts-1]  = pnts_xsecs( i, platePnts-1 ).transform(mat);
		nCVec[platePnts-1] = vec3d(0,0,0); // on camber line
		tVec[platePnts-1]  = 0;
		zVec[platePnts-1]  = 0;

		//== Compute plate normal ==//
		// rotated chord vector
		vec3d rcv = pnts_xsecs(i, platePnts-1).transform(mat) - pnts_xsecs(i,0).transform(mat);
		// rotated area normal vector
		vec3d anv = get_area_normal(i).transform(mat) - vec3d(0,0,0).transform(mat);
		// plate normal vector
		nPlate = cross(rcv,anv);
		nPlate.normalize();
		degenPlate.nPlate.push_back( nPlate );

		// normalized, unrotated chord vector (te->le)
		chordVec = pnts_xsecs(i, platePnts-1) - pnts_xsecs(i,0);
		chordVec.normalize();

		// Set points along af section
		for ( int j = 1, k = num_pnts-2; j < platePnts-1; j++, k-- )
		{
			topPnt = pnts_xsecs(i,j);
			botPnt = pnts_xsecs(i,k);

			camberPnt = (topPnt + botPnt) / 2;
			chordPnt  = pnts_xsecs(i,0) + chordVec * dot(camberPnt-pnts_xsecs(i,0),chordVec);

			xVec[j]  = chordPnt.transform(mat);
			nCVec[j] = topPnt.transform(mat)-botPnt.transform(mat);
			nCVec[j].normalize();
			tVec[j]  = dist(topPnt,botPnt);
			zVec[j]  = dist(camberPnt,chordPnt);
			if ( angle(camberPnt - chordPnt, nPlate) != 0) zVec[j] *= -1;
		}

		xMat.push_back(xVec);
		nCamberMat.push_back(nCVec);
		tMat.push_back(tVec);
		zMat.push_back(zVec);

		degenPlate.u.push_back( uArray[i] );
	}

	for ( int j = 0, k = num_pnts-1; j < platePnts-1; j++, k-- )
	{
		degenPlate.wTop.push_back( wArray[j] );
		degenPlate.wBot.push_back( wArray[k] );
	}
	degenPlate.wTop.push_back( wArray[platePnts-1] );
	degenPlate.wBot.push_back( wArray[platePnts-1] );

	// "Horizontal" Plate
	int startPnt = (num_pnts - 1) / 4;
	for ( int i = nLow; i < nHigh; i++ )
	{
		// Set first point (trailing edge)
		xVec[0]  = pnts_xsecs(i,startPnt).transform(mat);
		nCVec[0] = vec3d(0,0,0); // on camber line
		tVec[0]  = 0;
		zVec[0]  = 0;

		// Set last point (leading edge)
		xVec[platePnts-1]  = pnts_xsecs( i, startPnt+platePnts-1 ).transform(mat);
		nCVec[platePnts-1] = vec3d(0,0,0); // on camber line
		tVec[platePnts-1]  = 0;
		zVec[platePnts-1]  = 0;

		//== Compute plate normal ==//
		// rotated chord vector
		vec3d rcv = pnts_xsecs(i, startPnt+platePnts-1).transform(mat) \
				  - pnts_xsecs(i,startPnt).transform(mat);
		// rotated area normal vector
		vec3d anv = get_area_normal(i).transform(mat) - vec3d(0,0,0).transform(mat);
		// plate normal vector
		nPlate = cross(rcv,anv);
		nPlate.normalize();
		degenPlate.nPlate.push_back( nPlate );

		// normalized, unrotated chord vector (te->le)
		chordVec = pnts_xsecs(i, startPnt+platePnts-1) - pnts_xsecs(i,startPnt);
		chordVec.normalize();

		// Set points along af section
		for ( int j = 1, k = startPnt+num_pnts-2; j < platePnts-1; j++, k-- )
		{
			topPnt = pnts_xsecs(i,startPnt+j);
			botPnt = pnts_xsecs(i, k % (num_pnts-1) );

			camberPnt = (topPnt + botPnt) / 2;
			chordPnt  = pnts_xsecs(i,startPnt) + chordVec * dot(camberPnt-pnts_xsecs(i,startPnt),chordVec);

			xVec[j]  = chordPnt.transform(mat);
			nCVec[j] = topPnt.transform(mat)-botPnt.transform(mat);
			nCVec[j].normalize();
			tVec[j]  = dist(topPnt,botPnt);
			zVec[j]  = dist(camberPnt,chordPnt);
			if ( angle(camberPnt - chordPnt, nPlate) != 0) zVec[j] *= -1;
		}

		xMat.push_back(xVec);
		nCamberMat.push_back(nCVec);
		tMat.push_back(tVec);
		zMat.push_back(zVec);

		degenPlate.u.push_back( uArray[i] );
	}


	for ( int j = startPnt, k = startPnt+num_pnts-1; j < startPnt+platePnts-1; j++, k-- )
	{
		degenPlate.wTop.push_back( wArray[j] );
		degenPlate.wBot.push_back( wArray[ k % (num_pnts-1) ] );
	}
	degenPlate.wTop.push_back( wArray[startPnt+platePnts-1] );
	degenPlate.wBot.push_back( wArray[startPnt+platePnts-1] );

	degenPlate.x    	= xMat;
	degenPlate.nCamber 	= nCamberMat;
	degenPlate.t		= tMat;
	degenPlate.zcamber	= zMat;
	degenGeom->setDegenPlate(degenPlate);

	if ( sym_code_in != NO_SYM )
		createBodyDegenPlate_refl(degenGeom, sym_code_in, mat, refl_mat);

	return;
}

void Xsec_surf::createBodyDegenPlate_refl(DegenGeom* degenGeom, int sym_code_in, float mat[4][4], float refl_mat[4][4])
{
	DegenPlate	degenPlate = degenGeom->getDegenPlate();

	int nLow = 0, nHigh = num_xsecs;
	int platePnts = (num_pnts + 1) / 2;

	vector< vector<vec3d>  >	xMat = degenPlate.x;
	vector< vector<vec3d>  >	nCamberMat = degenPlate.nCamber;
	vector< vector<double> >	tMat = degenPlate.t;
	vector< vector<double> >	zMat = degenPlate.zcamber;

	vector<vec3d>	xVec(  platePnts );
	vector<vec3d>	nCVec( platePnts );
	vector<vec3d>	nPVec( platePnts );
	vector<double>	tVec(  platePnts );
	vector<double>	zVec(  platePnts );

	vec3d  topPnt, botPnt, chordVec, camberPnt, chordPnt, nPlate;

	if ( sym_code_in != refl_pnts_xsecs_code )
	{
		refl_pnts_xsecs_code = sym_code_in;
		load_refl_pnts_xsecs();
	}

	// "Vertical" plate
	for ( int i = nLow; i < nHigh; i++ )
	{
		// Set first point (trailing edge)
		xVec[0]  = refl_pnts_xsecs(i,0).transform(refl_mat);
		nCVec[0] = vec3d(0,0,0); // on camber line
		tVec[0]  = 0;
		zVec[0]  = 0;

		// Set last point (leading edge)
		xVec[platePnts-1]  = refl_pnts_xsecs( i, platePnts-1 ).transform(refl_mat);
		nCVec[platePnts-1] = vec3d(0,0,0); // on camber line
		tVec[platePnts-1]  = 0;
		zVec[platePnts-1]  = 0;

		//== Compute plate normal ==//
		// rotated chord vector
		vec3d rcv = refl_pnts_xsecs(i, platePnts-1).transform(refl_mat) \
				  - refl_pnts_xsecs(i,0).transform(refl_mat);
		// rotated area normal vector
		vec3d anv = get_refl_area_normal(i).transform(refl_mat) - vec3d(0,0,0).transform(refl_mat);
		// plate normal vector
		nPlate = cross(rcv,anv);
		nPlate.normalize();
		degenPlate.nPlate.push_back( nPlate );

		// normalized, unrotated chord vector (te->le)
		chordVec = refl_pnts_xsecs(i, platePnts-1) - refl_pnts_xsecs(i,0);
		chordVec.normalize();

		// Set points along af section
		for ( int j = 1, k = num_pnts-2; j < platePnts-1; j++, k-- )
		{
			topPnt = refl_pnts_xsecs(i,j);
			botPnt = refl_pnts_xsecs(i,k);

			camberPnt = (topPnt + botPnt) / 2;
			chordPnt  = refl_pnts_xsecs(i,0) + chordVec * dot(camberPnt-refl_pnts_xsecs(i,0),chordVec);

			xVec[j]  = chordPnt.transform(refl_mat);
			nCVec[j] = topPnt.transform(refl_mat)-botPnt.transform(refl_mat);
			nCVec[j].normalize();
			tVec[j]  = dist(topPnt,botPnt);
			zVec[j]  = dist(camberPnt,chordPnt);
			if ( angle(camberPnt - chordPnt, nPlate) != 0) zVec[j] *= -1;
		}

		xMat.push_back(xVec);
		nCamberMat.push_back(nCVec);
		tMat.push_back(tVec);
		zMat.push_back(zVec);

		degenPlate.u.push_back( uArray[i] );
	}

	for ( int j = 0, k = num_pnts-1; j < platePnts-1; j++, k-- )
	{
		degenPlate.wTop.push_back( wArray[j] );
		degenPlate.wBot.push_back( wArray[k] );
	}
	degenPlate.wTop.push_back( wArray[platePnts-1] );
	degenPlate.wBot.push_back( wArray[platePnts-1] );

	// "Horizontal" Plate
	int startPnt = (num_pnts - 1) / 4;
	for ( int i = nLow; i < nHigh; i++ )
	{
		// Set first point (trailing edge)
		xVec[0]  = refl_pnts_xsecs(i,startPnt).transform(refl_mat);
		nCVec[0] = vec3d(0,0,0); // on camber line
		tVec[0]  = 0;
		zVec[0]  = 0;

		// Set last point (leading edge)
		xVec[platePnts-1]  = refl_pnts_xsecs( i, startPnt+platePnts-1 ).transform(refl_mat);
		nCVec[platePnts-1] = vec3d(0,0,0); // on camber line
		tVec[platePnts-1]  = 0;
		zVec[platePnts-1]  = 0;

		//== Compute plate normal ==//
		// rotated chord vector
		vec3d rcv = refl_pnts_xsecs(i, startPnt+platePnts-1).transform(refl_mat) \
				  - refl_pnts_xsecs(i,startPnt).transform(refl_mat);
		// rotated area normal vector
		vec3d anv = get_refl_area_normal(i).transform(refl_mat) - vec3d(0,0,0).transform(refl_mat);
		// plate normal vector
		nPlate = cross(rcv,anv);
		nPlate.normalize();
		degenPlate.nPlate.push_back( nPlate );

		// normalized, unrotated chord vector (te->le)
		chordVec = refl_pnts_xsecs(i, startPnt+platePnts-1) - refl_pnts_xsecs(i,startPnt);
		chordVec.normalize();

		// Set points along af section
		for ( int j = 1, k = startPnt+num_pnts-2; j < platePnts-1; j++, k-- )
		{
			topPnt = refl_pnts_xsecs(i,startPnt+j);
			botPnt = refl_pnts_xsecs(i, k % (num_pnts-1) );

			camberPnt = (topPnt + botPnt) / 2;
			chordPnt  = refl_pnts_xsecs(i,startPnt) + chordVec * dot(camberPnt-refl_pnts_xsecs(i,startPnt),chordVec);

			xVec[j]  = chordPnt.transform(refl_mat);
			nCVec[j] = topPnt.transform(refl_mat)-botPnt.transform(refl_mat);
			nCVec[j].normalize();
			tVec[j]  = dist(topPnt,botPnt);
			zVec[j]  = dist(camberPnt,chordPnt);
			if ( angle(camberPnt - chordPnt, nPlate) != 0) zVec[j] *= -1;
		}

		xMat.push_back(xVec);
		nCamberMat.push_back(nCVec);
		tMat.push_back(tVec);
		zMat.push_back(zVec);

		degenPlate.u.push_back( uArray[i] );
	}

	for ( int j = startPnt, k = startPnt+num_pnts-1; j < startPnt+platePnts-1; j++, k-- )
	{
		degenPlate.wTop.push_back( wArray[j] );
		degenPlate.wBot.push_back( wArray[ k % (num_pnts-1) ] );
	}
	degenPlate.wTop.push_back( wArray[startPnt+platePnts-1] );
	degenPlate.wBot.push_back( wArray[startPnt+platePnts-1] );

	degenPlate.x    	= xMat;
	degenPlate.nCamber 	= nCamberMat;
	degenPlate.t		= tMat;
	degenPlate.zcamber	= zMat;
	degenGeom->setDegenPlate(degenPlate);
}

void Xsec_surf::createSurfDegenStick(DegenGeom* degenGeom, int sym_code_in, float mat[4][4], float refl_mat[4][4])
{
	DegenStick	degenStick = degenGeom->getDegenStick();

	int nLow = 0, nHigh = num_xsecs;
	int platePnts = (num_pnts + 1) / 2;
	vec3d chordVec, camberPnt, prevCamberPnt;

	if ( degenGeom->getParentGeom()->getTypeStr() == "wing" || degenGeom->getParentGeom()->getTypeStr() == "prop" )
	{
		// Keep only airfoil sections, discard endcap close-out lines
		nLow  = 1;
		nHigh = num_xsecs - 1;
	}

	for ( int i = nLow; i < nHigh; i++ )
	{
		vec3d topPnt, botPnt;
		double tempThickness = 0, perimTop = 0, perimBot = 0;
		int    maxThickIdx[2];

		// normalized, unrotated chord vector (te->le)
		chordVec = pnts_xsecs(i, platePnts-1) - pnts_xsecs(i,0);
		chordVec.normalize();

		degenStick.xle.push_back( pnts_xsecs( i, platePnts-1 ).transform(mat) );
		degenStick.xte.push_back( pnts_xsecs( i, 0 ).transform(mat) );
		degenStick.chord.push_back( dist(pnts_xsecs(i, platePnts-1), pnts_xsecs(i, 0)) );
		degenStick.u.push_back( uArray[i] );
		degenStick.Ishell.push_back(calculate_shell_inertias_in_plane(i,XZ_PLANE, mat));
		degenStick.Isolid.push_back(calculate_solid_inertias_in_plane(i,XZ_PLANE, mat));
		degenStick.xcgSolid.push_back( get_xsec_centroid(i).transform(mat) );
		degenStick.xcgShell.push_back( get_xsec_shellCG(i).transform(mat)  );
		degenStick.area.push_back( get_xsec_area(i) );

		vec3d areaNormal = get_area_normal(i).transform(mat) - vec3d(0,0,0).transform(mat);
		areaNormal.normalize();
		degenStick.areaNormal.push_back( areaNormal );

		for ( int j = 1, k = num_pnts-2; j < platePnts-1; j++, k-- )
		{
			topPnt = pnts_xsecs(i,j);
			botPnt = pnts_xsecs(i,k);

			camberPnt = ( topPnt + botPnt ) / 2;

			if( dist(topPnt, botPnt) > tempThickness)
			{
				tempThickness  = dist(topPnt, botPnt);
				maxThickIdx[0] = j;
				maxThickIdx[1] = k;
			}
			perimTop += dist( pnts_xsecs(i,j), pnts_xsecs(i,j-1) );
			perimBot += dist( pnts_xsecs(i,k), pnts_xsecs(i,k+1) );

			prevCamberPnt = camberPnt;
		}

		camberPnt = ( pnts_xsecs(i,maxThickIdx[0]) + pnts_xsecs(i,maxThickIdx[1]) ) / 2;
		degenStick.tLoc.push_back( 1 - (dot(camberPnt-pnts_xsecs(i,0),chordVec) / degenStick.chord.back()) );
		degenStick.toc.push_back( tempThickness / degenStick.chord.back() );

		perimTop += dist( pnts_xsecs(i, platePnts-1), pnts_xsecs(i, platePnts-2) );
		perimBot += dist( pnts_xsecs(i, platePnts), pnts_xsecs(i, platePnts-1) );
		degenStick.perimTop.push_back( perimTop );
		degenStick.perimBot.push_back( perimBot );
	}

	// Calculate sweep angle
	for( int i = nLow; i < nHigh - 1; i++ )
	{
		vec3d  cvCurrent, cvNext, qcCurrent, qcNext;
		double chordCurrent, chordNext;

		// Get current section chord vector
		cvCurrent = pnts_xsecs(i, platePnts-1) - pnts_xsecs(i,0);
		chordCurrent = cvCurrent.mag();
		cvCurrent.normalize();
		// Get current section quarter chord point
		qcCurrent = (pnts_xsecs(i,0) + cvCurrent*0.75*chordCurrent).transform(mat);

		// Get next section chord vector
		cvNext = pnts_xsecs(i+1, platePnts-1) - pnts_xsecs(i+1,0);
		chordNext = cvNext.mag();
		cvNext.normalize();
		// Get next section quarter chord point
		qcNext = (pnts_xsecs(i+1,0) + cvNext*0.75*chordNext).transform(mat);

		// Get vector from current to next quarter chord
		vec3d qcVec = qcNext - qcCurrent;

		// Zero out z component so angle is only in x-y plane
		qcVec.set_z(0);

		vec3d yAxis(0,1,0), zAxis(0,0,-1);
		if(qcVec.y() < 0)
		{
			yAxis.set_y(-1);
			zAxis.set_z(1);
		}

		// Get signed angle between qc vector and y axis
		double lambda =  ((double)180 / 3.1415927) * signed_angle(yAxis, qcVec, zAxis);
		degenStick.sweep.push_back( lambda );
	}
	degenStick.sweep.push_back( NAN );

	if ( sym_code_in == NO_SYM )
	{
		degenGeom->setDegenStick(degenStick);
		return;
	}

	for ( int i = nLow; i < nHigh; i++ )
	{
		vec3d topPnt, botPnt;
		double tempThickness = 0, perimTop = 0, perimBot = 0;
		int    maxThickIdx[2];

		// normalized, unrotated chord vector (te->le)
		chordVec = refl_pnts_xsecs(i, platePnts-1) - refl_pnts_xsecs(i,0);
		chordVec.normalize();

		degenStick.xle.push_back( refl_pnts_xsecs( i, platePnts-1 ).transform(refl_mat) );
		degenStick.xte.push_back( refl_pnts_xsecs( i, 0 ).transform(refl_mat) );
		degenStick.chord.push_back( dist(refl_pnts_xsecs(i, platePnts-1), refl_pnts_xsecs(i, 0)) );
		degenStick.u.push_back( uArray[i] );
		degenStick.Ishell.push_back(calculate_refl_shell_inertias_in_plane(i,XZ_PLANE, refl_mat));
		degenStick.Isolid.push_back(calculate_refl_solid_inertias_in_plane(i,XZ_PLANE, refl_mat));
		degenStick.xcgSolid.push_back( get_refl_xsec_centroid(i).transform(refl_mat) );
		degenStick.xcgShell.push_back( get_refl_xsec_shellCG(i).transform(refl_mat)  );
		degenStick.area.push_back( get_xsec_area(i) );

		vec3d areaNormal = get_refl_area_normal(i).transform(refl_mat) - vec3d(0,0,0).transform(refl_mat);
		areaNormal.normalize();
		degenStick.areaNormal.push_back( areaNormal );

		for ( int j = 1, k = num_pnts-2; j < platePnts-1; j++, k-- )
		{
			topPnt = refl_pnts_xsecs(i,j);
			botPnt = refl_pnts_xsecs(i,k);

			camberPnt = ( topPnt + botPnt ) / 2;

			if( dist(topPnt, botPnt) > tempThickness)
			{
				tempThickness  = dist(topPnt, botPnt);
				maxThickIdx[0] = j;
				maxThickIdx[1] = k;
			}
			perimTop += dist( refl_pnts_xsecs(i,j), refl_pnts_xsecs(i,j-1) );
			perimBot += dist( refl_pnts_xsecs(i,k), refl_pnts_xsecs(i,k+1) );

			prevCamberPnt = camberPnt;
		}

		camberPnt = ( refl_pnts_xsecs(i,maxThickIdx[0]) + refl_pnts_xsecs(i,maxThickIdx[1]) ) / 2;
		degenStick.tLoc.push_back( 1 - (dot(camberPnt-refl_pnts_xsecs(i,0),chordVec) / degenStick.chord.back() ) );
		degenStick.toc.push_back( tempThickness / degenStick.chord.back() );

		perimTop += dist( refl_pnts_xsecs(i, platePnts-1), refl_pnts_xsecs(i, platePnts-2) );
		perimBot += dist( refl_pnts_xsecs(i, platePnts), refl_pnts_xsecs(i, platePnts-1) );
		degenStick.perimTop.push_back( perimTop );
		degenStick.perimBot.push_back( perimBot );
	}

	// Calculate sweep angle
	for( int i = nLow; i < nHigh - 1; i++ )
	{
		vec3d  cvCurrent, cvNext, qcCurrent, qcNext;
		double chordCurrent, chordNext;

		// Get current section chord vector
		cvCurrent = refl_pnts_xsecs(i, platePnts-1) - refl_pnts_xsecs(i,0);
		chordCurrent = cvCurrent.mag();
		cvCurrent.normalize();
		// Get current section quarter chord point
		qcCurrent = (refl_pnts_xsecs(i,0) + cvCurrent*0.75*chordCurrent).transform(refl_mat);

		// Get next section chord vector
		cvNext = refl_pnts_xsecs(i+1, platePnts-1) - refl_pnts_xsecs(i+1,0);
		chordNext = cvNext.mag();
		cvNext.normalize();
		// Get next section quarter chord point
		qcNext = (refl_pnts_xsecs(i+1,0) + cvNext*0.75*chordNext).transform(refl_mat);

		// Get vector from current to next quarter chord
		vec3d qcVec = qcNext - qcCurrent;

		// Zero out z component so angle is only in x-y plane
		qcVec.set_z(0);

		vec3d yAxis(0,1,0), zAxis(0,0,-1);
		if(qcVec.y() < 0)
		{
			yAxis.set_y(-1);
			zAxis.set_z(1);
		}

		// Get signed angle between qc vector and y axis
		double lambda =  ((double)180 / 3.1415927) * signed_angle(yAxis, qcVec, zAxis);
		degenStick.sweep.push_back( lambda );
	}
	degenStick.sweep.push_back( NAN );

	degenGeom->setDegenStick( degenStick );

}

void Xsec_surf::createBodyDegenStick(DegenGeom* degenGeom, int sym_code_in, float mat[4][4], float refl_mat[4][4])
{
	DegenStick	degenStick = degenGeom->getDegenStick();

	int nLow = 0, nHigh = num_xsecs;
	int platePnts = (num_pnts + 1) / 2;
	vec3d chordVec, camberPnt, prevCamberPnt;

	// "Vertical"
	for ( int i = nLow; i < nHigh; i++ )
	{
		vec3d topPnt, botPnt;
		double tempThickness = 0, perimTop = 0, perimBot = 0;
		int    maxThickIdx[2] = {0,0};

		// normalized, unrotated chord vector (te->le)
		chordVec = pnts_xsecs(i, platePnts-1) - pnts_xsecs(i,0);
		chordVec.normalize();

		degenStick.xle.push_back( pnts_xsecs( i, platePnts-1 ).transform(mat) );
		degenStick.xte.push_back( pnts_xsecs( i, 0 ).transform(mat) );
		degenStick.chord.push_back( dist(pnts_xsecs(i, platePnts-1), pnts_xsecs(i, 0)) );
		degenStick.u.push_back( uArray[i] );
		degenStick.Ishell.push_back(calculate_shell_inertias_in_plane(i,YZ_PLANE, mat));
		degenStick.Isolid.push_back(calculate_solid_inertias_in_plane(i,YZ_PLANE, mat));
		degenStick.xcgSolid.push_back( get_xsec_centroid(i).transform(mat) );
		degenStick.xcgShell.push_back( get_xsec_shellCG(i).transform(mat)  );
		degenStick.area.push_back( get_xsec_area(i) );

		vec3d areaNormal = get_area_normal(i).transform(mat) - vec3d(0,0,0).transform(mat);
		areaNormal.normalize();
		degenStick.areaNormal.push_back( areaNormal );

		for ( int j = 1, k = num_pnts-2; j < platePnts-1; j++, k-- )
		{

			topPnt = pnts_xsecs(i,j);
			botPnt = pnts_xsecs(i,k);

			camberPnt = ( topPnt + botPnt ) / 2;

			if( dist(topPnt, botPnt) > tempThickness)
			{
				tempThickness  = dist(topPnt, botPnt);
				maxThickIdx[0] = j;
				maxThickIdx[1] = k;
			}
			perimTop += dist( pnts_xsecs(i,j), pnts_xsecs(i,j-1) );
			perimBot += dist( pnts_xsecs(i,k), pnts_xsecs(i,k+1) );

			prevCamberPnt = camberPnt;
		}

		camberPnt = ( pnts_xsecs(i,maxThickIdx[0]) + pnts_xsecs(i,maxThickIdx[1]) ) / 2;
		degenStick.tLoc.push_back( 1 - (dot(camberPnt-pnts_xsecs(i,0),chordVec) / degenStick.chord.back()) );
		degenStick.toc.push_back( tempThickness / degenStick.chord.back() );

		perimTop += dist( pnts_xsecs(i, platePnts-1), pnts_xsecs(i, platePnts-2) );
		perimBot += dist( pnts_xsecs(i, platePnts), pnts_xsecs(i, platePnts-1) );
		degenStick.perimTop.push_back( perimTop );
		degenStick.perimBot.push_back( perimBot );
	}

	// Calculate sweep angle
	for( int i = nLow; i < nHigh; i++ )
	{
		degenStick.sweep.push_back( NAN );
	}

	// "Horizontal"
	int startPnt = (num_pnts - 1) / 4;
	for ( int i = nLow; i < nHigh; i++ )
	{
		vec3d topPnt, botPnt;
		double tempThickness = 0, perimTop = 0, perimBot = 0;
		int    maxThickIdx[2] = {0,0};

		// normalized, unrotated chord vector (te->le)
		chordVec = pnts_xsecs(i, startPnt+platePnts-1) - pnts_xsecs(i, startPnt);
		chordVec.normalize();

		degenStick.xle.push_back( pnts_xsecs( i, startPnt+platePnts-1 ).transform(mat) );
		degenStick.xte.push_back( pnts_xsecs( i, startPnt ).transform(mat) );
		degenStick.chord.push_back( dist(pnts_xsecs(i, startPnt+platePnts-1), pnts_xsecs(i, startPnt)) );
		degenStick.u.push_back( uArray[i] );
		degenStick.Ishell.push_back(calculate_shell_inertias_in_plane(i,YZ_PLANE, mat));
		degenStick.Isolid.push_back(calculate_solid_inertias_in_plane(i,YZ_PLANE, mat));
		degenStick.xcgSolid.push_back( get_xsec_centroid(i).transform(mat) );
		degenStick.xcgShell.push_back( get_xsec_shellCG(i).transform(mat)  );
		degenStick.area.push_back( get_xsec_area(i) );

		vec3d areaNormal = get_area_normal(i).transform(mat) - vec3d(0,0,0).transform(mat);
		areaNormal.normalize();
		degenStick.areaNormal.push_back( areaNormal );

		for ( int j = 1, k = startPnt+num_pnts-2; j < platePnts-1; j++, k-- )
		{
			topPnt = pnts_xsecs(i,startPnt+j);
			botPnt = pnts_xsecs(i, k % (num_pnts-1) );

			camberPnt = ( topPnt + botPnt ) / 2;

			if( dist(topPnt, botPnt) > tempThickness)
			{
				tempThickness  = dist(topPnt, botPnt);
				maxThickIdx[0] = j+startPnt;
				maxThickIdx[1] = k % (num_pnts-1);
			}
			perimTop += dist( pnts_xsecs(i,startPnt+j), pnts_xsecs(i,startPnt+j-1) );
			perimBot += dist( pnts_xsecs(i, k % (num_pnts-1) ), pnts_xsecs(i, k % (num_pnts-1)+1) );

			prevCamberPnt = camberPnt;
		}

		camberPnt = ( pnts_xsecs(i,maxThickIdx[0]) + pnts_xsecs(i,maxThickIdx[1]) ) / 2;
		degenStick.tLoc.push_back( 1 - (dot(camberPnt-pnts_xsecs(i,startPnt),chordVec) / degenStick.chord.back()) );
		degenStick.toc.push_back( tempThickness / degenStick.chord.back() );

		perimTop += dist( pnts_xsecs(i, startPnt+platePnts-1), pnts_xsecs(i, startPnt+platePnts-2) );
		perimBot += dist( pnts_xsecs(i, startPnt+platePnts), pnts_xsecs(i, startPnt+platePnts-1) );
		degenStick.perimTop.push_back( perimTop );
		degenStick.perimBot.push_back( perimBot );
	}

	// Calculate sweep angle
	for( int i = nLow; i < nHigh; i++ )
	{
		degenStick.sweep.push_back( NAN );
	}

	if ( sym_code_in == NO_SYM )
	{
		degenGeom->setDegenStick(degenStick);
		return;
	}

	// "Vertical"
	for ( int i = nLow; i < nHigh; i++ )
	{
		vec3d topPnt, botPnt;
		double tempThickness = 0, perimTop = 0, perimBot = 0;
		int    maxThickIdx[2] = {0,0};

		// normalized, unrotated chord vector (te->le)
		chordVec = refl_pnts_xsecs(i, platePnts-1) - refl_pnts_xsecs(i,0);
		chordVec.normalize();

		degenStick.xle.push_back( refl_pnts_xsecs( i, platePnts-1 ).transform(refl_mat) );
		degenStick.xte.push_back( refl_pnts_xsecs( i, 0 ).transform(refl_mat) );
		degenStick.chord.push_back( dist(refl_pnts_xsecs(i, platePnts-1), refl_pnts_xsecs(i, 0)) );
		degenStick.u.push_back( uArray[i] );
		degenStick.Ishell.push_back(calculate_refl_shell_inertias_in_plane(i,YZ_PLANE, refl_mat));
		degenStick.Isolid.push_back(calculate_refl_solid_inertias_in_plane(i,YZ_PLANE, refl_mat));
		degenStick.xcgSolid.push_back( get_refl_xsec_centroid(i).transform(refl_mat) );
		degenStick.xcgShell.push_back( get_refl_xsec_shellCG(i).transform(refl_mat)  );
		degenStick.area.push_back( get_xsec_area(i) );

		vec3d areaNormal = get_refl_area_normal(i).transform(refl_mat) - vec3d(0,0,0).transform(refl_mat);
		areaNormal.normalize();
		degenStick.areaNormal.push_back( areaNormal );

		for ( int j = 1, k = num_pnts-2; j < platePnts-1; j++, k-- )
		{
			topPnt = refl_pnts_xsecs(i,j);
			botPnt = refl_pnts_xsecs(i,k);

			camberPnt = ( topPnt + botPnt ) / 2;

			if( dist(topPnt, botPnt) > tempThickness)
			{
				tempThickness  = dist(topPnt, botPnt);
				maxThickIdx[0] = j;
				maxThickIdx[1] = k;
			}
			perimTop += dist( refl_pnts_xsecs(i,j), refl_pnts_xsecs(i,j-1) );
			perimBot += dist( refl_pnts_xsecs(i,k), refl_pnts_xsecs(i,k+1) );

			prevCamberPnt = camberPnt;
		}

		camberPnt = ( refl_pnts_xsecs(i,maxThickIdx[0]) + refl_pnts_xsecs(i,maxThickIdx[1]) ) / 2;
		degenStick.tLoc.push_back( 1 - (dot(camberPnt-refl_pnts_xsecs(i,0),chordVec) / degenStick.chord.back()) );
		degenStick.toc.push_back( tempThickness / degenStick.chord.back() );

		perimTop += dist( refl_pnts_xsecs(i, platePnts-1), refl_pnts_xsecs(i, platePnts-2) );
		perimBot += dist( refl_pnts_xsecs(i, platePnts), refl_pnts_xsecs(i, platePnts-1) );
		degenStick.perimTop.push_back( perimTop );
		degenStick.perimBot.push_back( perimBot );
	}

	// Calculate sweep angle
	for( int i = nLow; i < nHigh; i++ )
	{
		degenStick.sweep.push_back( NAN );
	}

	// "Horizontal"
	for ( int i = nLow; i < nHigh; i++ )
	{
		vec3d topPnt, botPnt;
		double tempThickness = 0, perimTop = 0, perimBot = 0;
		int    maxThickIdx[2] = {0,0};

		// normalized, unrotated chord vector (te->le)
		chordVec = refl_pnts_xsecs(i, startPnt+platePnts-1) - refl_pnts_xsecs(i, startPnt);
		chordVec.normalize();

		degenStick.xle.push_back( refl_pnts_xsecs( i, startPnt+platePnts-1 ).transform(refl_mat) );
		degenStick.xte.push_back( refl_pnts_xsecs( i, startPnt ).transform(refl_mat) );
		degenStick.chord.push_back( dist(refl_pnts_xsecs(i, startPnt+platePnts-1), refl_pnts_xsecs(i, startPnt)) );
		degenStick.u.push_back( uArray[i] );
		degenStick.Ishell.push_back(calculate_refl_shell_inertias_in_plane(i,YZ_PLANE, refl_mat));
		degenStick.Isolid.push_back(calculate_refl_solid_inertias_in_plane(i,YZ_PLANE, refl_mat));
		degenStick.xcgSolid.push_back( get_refl_xsec_centroid(i).transform(refl_mat) );
		degenStick.xcgShell.push_back( get_refl_xsec_shellCG(i).transform(refl_mat)  );
		degenStick.area.push_back( get_xsec_area(i) );

		vec3d areaNormal = get_refl_area_normal(i).transform(refl_mat) - vec3d(0,0,0).transform(refl_mat);
		areaNormal.normalize();
		degenStick.areaNormal.push_back( areaNormal );

		for ( int j = 1, k = startPnt+num_pnts-2; j < platePnts-1; j++, k-- )
		{
			topPnt = refl_pnts_xsecs(i,startPnt+j);
			botPnt = refl_pnts_xsecs(i, k % (num_pnts-1) );

			camberPnt = ( topPnt + botPnt ) / 2;

			if( dist(topPnt, botPnt) > tempThickness)
			{
				tempThickness  = dist(topPnt, botPnt);
				maxThickIdx[0] = j;
				maxThickIdx[1] = k % (num_pnts-1);
			}
			perimTop += dist( refl_pnts_xsecs(i,startPnt+j), refl_pnts_xsecs(i,startPnt+j-1) );
			perimBot += dist( refl_pnts_xsecs(i, k % (num_pnts-1) ), refl_pnts_xsecs(i, k % (num_pnts-1)+1) );

			prevCamberPnt = camberPnt;
		}

		camberPnt = ( refl_pnts_xsecs(i,maxThickIdx[0]) + refl_pnts_xsecs(i,maxThickIdx[1]) ) / 2;
		degenStick.tLoc.push_back( 1 - (dot(camberPnt-refl_pnts_xsecs(i,startPnt),chordVec) / degenStick.chord.back()) );
		degenStick.toc.push_back( tempThickness / degenStick.chord.back() );

		perimTop += dist( refl_pnts_xsecs(i, startPnt+platePnts-1), refl_pnts_xsecs(i, startPnt+platePnts-2) );
		perimBot += dist( refl_pnts_xsecs(i, startPnt+platePnts), refl_pnts_xsecs(i, startPnt+platePnts-1) );
		degenStick.perimTop.push_back( perimTop );
		degenStick.perimBot.push_back( perimBot );
	}

	// Calculate sweep angle
	for( int i = nLow; i < nHigh; i++ )
	{
		degenStick.sweep.push_back( NAN );
	}

	degenGeom->setDegenStick( degenStick );

}

vec3d Xsec_surf::get_xsec_shellCG( int ixs )
{
	double xcg = 0, zcg = 0, area = get_xsec_area(ixs);
	int j;
	vec3d xAxis(1,0,0), zAxis(0,0,1), areaNormal = get_area_normal(ixs);

	vector<vec3d> pnts;
	// Load cross section points
	for( int i = 0; i < num_pnts; i++ )
	{
		pnts.push_back( pnts_xsecs(ixs,i) );
	}

	// Get rotation in xy plane to align areaNormal with yaxis
	vec2d an2(areaNormal.x(), areaNormal.y()), ax1(1,0), ax2(0,1);
	double theta2, theta1 = angle(an2, ax2);
	if (areaNormal.x() < 0) theta1 = -theta1;
	// Rotate areaNormal around z to y axis
	areaNormal   = RotateArbAxis( areaNormal, theta1, zAxis );
	// Rotate points as well
	for ( int i = 0; i < num_pnts; i++ )
		pnts[i] = RotateArbAxis( pnts[i], theta1, zAxis );

	// Get rotation in yz plane to align areaNormal with yaxis
	an2.set_xy(areaNormal.y(), areaNormal.z());
	theta2 = angle(an2, ax1);
	if(areaNormal.z() >= 0) theta2 = -theta2;
	// Rotate points about x axis
	for ( int i = 0; i < num_pnts; i++ )
		pnts[i] = RotateArbAxis( pnts[i], theta2, xAxis );

	vec3d  cg;
	double segLen, totLen = 0;
	// Cross section should now be in x-z plane. Calculate shell cg.
	for ( j = 1; j < num_pnts; j++ )
	{
		// Line segment center point
		cg     = (pnts[j-1] + pnts[j])/2;
		// Line segment length
		segLen = (pnts[j-1] - pnts[j]).mag();
		// Sum of segment lengths (proxy for area, which is proxy for mass)
		totLen += segLen;

		xcg += cg.x() * segLen;
		zcg += cg.z() * segLen;
	}

	// Line segment center point
	cg     = (pnts[j-1] + pnts[0])/2;
	// Line segment length
	segLen = (pnts[j-1] - pnts[0]).mag();
	// Sum of segment lengths (proxy for area, which is proxy for mass)
	totLen += segLen;

	xcg += cg.x() * segLen;
	zcg += cg.z() * segLen;

	xcg /= totLen;
	zcg /= totLen;

	// Xcg vector including y position (same for all points since xsec rotated into xz plane)
	vec3d cgLoc(xcg, pnts[0].y(), zcg);
	// Rotate back into original coordinate frame
	cgLoc = RotateArbAxis( cgLoc, -theta2, xAxis );
	cgLoc = RotateArbAxis( cgLoc, -theta1, zAxis );

	return cgLoc;
}

vec3d Xsec_surf::get_refl_xsec_shellCG( int ixs )
{
	double xcg = 0, zcg = 0, area = get_xsec_area(ixs);
	int j;
	vec3d xAxis(1,0,0), zAxis(0,0,1), areaNormal = get_refl_area_normal(ixs);

	vector<vec3d> pnts;
	// Load cross section points
	for( int i = 0; i < num_pnts; i++ )
	{
		pnts.push_back( refl_pnts_xsecs(ixs,i) );
	}

	// Get rotation in xy plane to align areaNormal with yaxis
	vec2d an2(areaNormal.x(), areaNormal.y()), ax1(1,0), ax2(0,1);
	double theta2, theta1 = angle(an2, ax2);
	if (areaNormal.x() < 0) theta1 = -theta1;
	// Rotate areaNormal around z to y axis
	areaNormal   = RotateArbAxis( areaNormal, theta1, zAxis );
	// Rotate points as well
	for ( int i = 0; i < num_pnts; i++ )
		pnts[i] = RotateArbAxis( pnts[i], theta1, zAxis );

	// Get rotation in yz plane to align areaNormal with yaxis
	an2.set_xy(areaNormal.y(), areaNormal.z());
	theta2 = angle(an2, ax1);
	if(areaNormal.z() >= 0) theta2 = -theta2;
	// Rotate points about x axis
	for ( int i = 0; i < num_pnts; i++ )
		pnts[i] = RotateArbAxis( pnts[i], theta2, xAxis );

	vec3d  cg;
	double segLen, totLen = 0;
	// Cross section should now be in x-z plane. Calculate shell cg.
	for ( j = 1; j < num_pnts; j++ )
	{
		// Line segment center point
		cg     = (pnts[j-1] + pnts[j])/2;
		// Line segment length
		segLen = (pnts[j-1] - pnts[j]).mag();
		// Sum of segment lengths (proxy for area, which is proxy for mass)
		totLen += segLen;

		xcg += cg.x() * segLen;
		zcg += cg.z() * segLen;
	}

	// Line segment center point
	cg     = (pnts[j-1] + pnts[0])/2;
	// Line segment length
	segLen = (pnts[j-1] - pnts[0]).mag();
	// Sum of segment lengths (proxy for area, which is proxy for mass)
	totLen += segLen;

	xcg += cg.x() * segLen;
	zcg += cg.z() * segLen;

	xcg /= totLen;
	zcg /= totLen;

	// Xcg vector including y position (same for all points since xsec rotated into xz plane)
	vec3d cgLoc(xcg, pnts[0].y(), zcg);
	// Rotate back into original coordinate frame
	cgLoc = RotateArbAxis( cgLoc, -theta2, xAxis );
	cgLoc = RotateArbAxis( cgLoc, -theta1, zAxis );

	return cgLoc;
}

//===== Jxx, Jzz, Jyy coefficients. Use Jxx = inertias[0]*t^3 + inertias[1]*t, =====================//
//===== Jzz = inertias[2]*t^3 + inertias[3]*t, etc with t as shell thickness to get inertias. =====//
vector<double> Xsec_surf::calculate_shell_inertias(int ixs)
{
	int j, platePnts = (num_pnts + 1) / 2;
	vector<double> inertias(6,0);
	vec3d xAxis(1,0,0), zAxis(0,0,1);
	vec3d areaNormal = get_area_normal(ixs), CG = get_xsec_centroid(ixs);

	vector<vec3d> pnts;
	// Load cross section points
	for( int i = 0; i < num_pnts; i++ )
	{
		pnts.push_back( pnts_xsecs(ixs,i) );
	}

	// Get rotation in xy plane to align areaNormal with yaxis
	vec2d an2(areaNormal.x(), areaNormal.y()), ax1(1,0), ax2(0,1);
	double theta = angle(an2, ax2);
	if (areaNormal.x() < 0) theta = -theta;
	// Rotate areaNormal around z to y axis
	areaNormal = RotateArbAxis( areaNormal, theta, zAxis );
	// Rotate centroid
	CG   = RotateArbAxis( CG, theta, zAxis );
	// Rotate points as well
	for ( int i = 0; i < num_pnts; i++ )
		pnts[i] = RotateArbAxis( pnts[i], theta, zAxis );

	// Get rotation in yz plane to align areaNormal with yaxis
	an2.set_xy(areaNormal.y(), areaNormal.z());
	theta = angle(an2, ax1);
	if(areaNormal.z() >= 0) theta = -theta;
	// Rotate points about x axis
	for ( int i = 0; i < num_pnts; i++ )
		pnts[i] = RotateArbAxis( pnts[i], theta, xAxis );
	// Rotate centroid
	CG = RotateArbAxis( CG, theta, xAxis );

	// Cross section is now in x-z plane. Compute shell inertias.
	vec3d segVec, chordVec = pnts[0] - pnts[platePnts-1];
	vec3d distToCG, cg;
	chordVec.normalize();
	double segLen, phi;

	for ( j = 1; j < num_pnts; j++ )
	{
		// Line segment center point
		cg     = (pnts[j-1] + pnts[j])/2;
		// Line segment vector
		segVec = pnts[j-1] - pnts[j];
		// Line segment length
		segLen = segVec.mag();
		segVec.normalize();
		// vector from line segment center point to xsec centroid
		distToCG = cg - CG;
		// angle between line segment and "chord"
		phi = angle(segVec,chordVec);

		inertias[0] += segLen/24 * ( 1+cos(2*phi) );
		inertias[1] += pow(segLen,3)/24 * ( 1-cos(2*phi) ) + segLen*pow(distToCG.z(),2);

		inertias[2] += segLen/24 * ( 1-cos(2*phi) );
		inertias[3] += pow(segLen,3)/24 * ( 1+cos(2*phi) ) + segLen*pow(distToCG.x(),2);
	}

	// Line segment center point
	cg     = (pnts[j-1] + pnts[0])/2;
	// Line segment vector
	segVec = pnts[j-1] - pnts[0];
	// Line segment length
	segLen = segVec.mag();
	segVec.normalize();
	// vector from line segment center point to xsec centroid
	distToCG = cg - CG;
	// angle between line segment and "chord"
	phi = angle(segVec,chordVec);

	inertias[0] += segLen/24 * ( 1+cos(2*phi) );
	inertias[1] += pow(segLen,3)/24 * ( 1-cos(2*phi) ) + segLen*pow(distToCG.z(),2);

	inertias[2] += segLen/24 * ( 1-cos(2*phi) );
	inertias[3] += pow(segLen,3)/24 * ( 1+cos(2*phi) ) + segLen*pow(distToCG.x(),2);

	inertias[4] = inertias[0] + inertias[2];
	inertias[5] = inertias[1] + inertias[3];

	return inertias;
}

vector<double> Xsec_surf::calculate_shell_inertias_in_plane(int ixs, int plane, float mat[4][4])
{
	int j, platePnts = (num_pnts + 1) / 2;
	vector<double> inertias(6,0);
	vector<vec2d> pnts;
	vec2d tmpPnts;

	switch(plane)
	{
	case XY_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( pnts_xsecs(ixs,i).transform(mat).x(), pnts_xsecs(ixs,i).transform(mat).y() );
			pnts.push_back(tmpPnts);
		}
		break;
	case XZ_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( pnts_xsecs(ixs,i).transform(mat).x(), pnts_xsecs(ixs,i).transform(mat).z() );
			pnts.push_back(tmpPnts);
		}
		break;
	case YZ_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( pnts_xsecs(ixs,i).transform(mat).y(), pnts_xsecs(ixs,i).transform(mat).z() );
			pnts.push_back(tmpPnts);
		}
		break;
	}

	vec2d segVec, chordVec = pnts[0] - pnts[platePnts-1];
	vec2d distToCG, cg, CG = get_xsec_centroid_in_plane(ixs, plane, mat);
	chordVec.normalize();
	double segLen, phi;

	for ( j = 1; j < num_pnts; j++ )
	{
		// Line segment center point
		cg     = (pnts[j-1] + pnts[j])/2;
		// Line segment vector
		segVec = pnts[j-1] - pnts[j];
		// Line segment length
		segLen = segVec.mag();
		segVec.normalize();
		// vector from line segment center point to xsec centroid
		distToCG = cg - CG;
		// angle between line segment and x-axis
		phi = angle(segVec, vec2d(1,0) );

		inertias[0] += segLen/24 * ( 1+cos(2*phi) );
		inertias[1] += pow(segLen,3)/24 * ( 1-cos(2*phi) ) + segLen*pow(distToCG.y(),2);

		inertias[2] += segLen/24 * ( 1-cos(2*phi) );
		inertias[3] += pow(segLen,3)/24 * ( 1+cos(2*phi) ) + segLen*pow(distToCG.x(),2);
	}

	// Line segment center point
	cg     = (pnts[j-1] + pnts[0])/2;
	// Line segment vector
	segVec = pnts[j-1] - pnts[0];
	// Line segment length
	segLen = segVec.mag();
	segVec.normalize();
	// vector from line segment center point to xsec centroid
	distToCG = cg - CG;
	// angle between line segment and x-axis
	phi = angle(segVec, vec2d(1,0) );

	inertias[0] += segLen/24 * ( 1+cos(2*phi) );
	inertias[1] += pow(segLen,3)/24 * ( 1-cos(2*phi) ) + segLen*pow(distToCG.y(),2);

	inertias[2] += segLen/24 * ( 1-cos(2*phi) );
	inertias[3] += pow(segLen,3)/24 * ( 1+cos(2*phi) ) + segLen*pow(distToCG.x(),2);

	inertias[4] = inertias[0] + inertias[2];
	inertias[5] = inertias[1] + inertias[3];

	return inertias;
}

vector<double> Xsec_surf::calculate_refl_shell_inertias_in_plane(int ixs, int plane, float refl_mat[4][4])
{
	int j, platePnts = (num_pnts + 1) / 2;
	vector<double> inertias(6,0);
	vector<vec2d> pnts;
	vec2d tmpPnts;

	switch(plane)
	{
	case XY_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( refl_pnts_xsecs(ixs,i).transform(refl_mat).x(), refl_pnts_xsecs(ixs,i).transform(refl_mat).y() );
			pnts.push_back(tmpPnts);
		}
		break;
	case XZ_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( refl_pnts_xsecs(ixs,i).transform(refl_mat).x(), refl_pnts_xsecs(ixs,i).transform(refl_mat).z() );
			pnts.push_back(tmpPnts);
		}
		break;
	case YZ_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( refl_pnts_xsecs(ixs,i).transform(refl_mat).y(), refl_pnts_xsecs(ixs,i).transform(refl_mat).z() );
			pnts.push_back(tmpPnts);
		}
		break;
	}

	vec2d segVec, chordVec = pnts[0] - pnts[platePnts-1];
	vec2d distToCG, cg, CG = get_refl_xsec_centroid_in_plane(ixs, plane, refl_mat);
	chordVec.normalize();
	double segLen, phi;

	for ( j = 1; j < num_pnts; j++ )
	{
		// Line segment center point
		cg     = (pnts[j-1] + pnts[j])/2;
		// Line segment vector
		segVec = pnts[j-1] - pnts[j];
		// Line segment length
		segLen = segVec.mag();
		segVec.normalize();
		// vector from line segment center point to xsec centroid
		distToCG = cg - CG;
		// angle between line segment and x-axis
		phi = angle(segVec, vec2d(1,0) );

		inertias[0] += segLen/24 * ( 1+cos(2*phi) );
		inertias[1] += pow(segLen,3)/24 * ( 1-cos(2*phi) ) + segLen*pow(distToCG.y(),2);

		inertias[2] += segLen/24 * ( 1-cos(2*phi) );
		inertias[3] += pow(segLen,3)/24 * ( 1+cos(2*phi) ) + segLen*pow(distToCG.x(),2);
	}

	// Line segment center point
	cg     = (pnts[j-1] + pnts[0])/2;
	// Line segment vector
	segVec = pnts[j-1] - pnts[0];
	// Line segment length
	segLen = segVec.mag();
	segVec.normalize();
	// vector from line segment center point to xsec centroid
	distToCG = cg - CG;
	// angle between line segment and x-axis
	phi = angle(segVec, vec2d(1,0) );

	inertias[0] += segLen/24 * ( 1+cos(2*phi) );
	inertias[1] += pow(segLen,3)/24 * ( 1-cos(2*phi) ) + segLen*pow(distToCG.y(),2);

	inertias[2] += segLen/24 * ( 1-cos(2*phi) );
	inertias[3] += pow(segLen,3)/24 * ( 1+cos(2*phi) ) + segLen*pow(distToCG.x(),2);

	inertias[4] = inertias[0] + inertias[2];
	inertias[5] = inertias[1] + inertias[3];

	return inertias;
}

vector<double> Xsec_surf::calculate_solid_inertias( int ixs )
{
	vector<double> inertias(2,0), inertiasCG(3,0);
	int j;
	double area = get_xsec_area(ixs);
	vec3d xAxis(1,0,0), zAxis(0,0,1);
	vec3d areaNormal = get_area_normal(ixs);
	vec3d centroid   = get_xsec_centroid(ixs);

	vector<vec3d> pnts;
	// Load cross section points
	for( int i = 0; i < num_pnts; i++ )
	{
		pnts.push_back( pnts_xsecs(ixs,i) );
	}

	// Get rotation in xy plane to align areaNormal with yaxis
	vec2d an2(areaNormal.x(), areaNormal.y()), ax1(1,0), ax2(0,1);
	double theta = angle(an2, ax2);
	if (areaNormal.x() < 0) theta = -theta;
	// Rotate areaNormal around z to y axis
	areaNormal = RotateArbAxis( areaNormal, theta, zAxis );
	// Rotate centroid
	centroid   = RotateArbAxis( centroid, theta, zAxis );
	// Rotate points as well
	for ( int i = 0; i < num_pnts; i++ )
		pnts[i] = RotateArbAxis( pnts[i], theta, zAxis );

	// Get rotation in yz plane to align areaNormal with yaxis
	an2.set_xy(areaNormal.y(), areaNormal.z());
	theta = angle(an2, ax1);
	if(areaNormal.z() >= 0) theta = -theta;
	// Rotate points about x axis
	for ( int i = 0; i < num_pnts; i++ )
		pnts[i] = RotateArbAxis( pnts[i], theta, xAxis );
	// Rotate centroid
	centroid = RotateArbAxis( centroid, theta, xAxis );

	// Cross section should now be in x-z plane. Calculate inertias
	for ( j = 1; j < num_pnts; j++ )
	{
		inertias[0] += (pnts[j-1].x()*pnts[j].z() - pnts[j].x()*pnts[j-1].z()) \
				* (pow(pnts[j-1].z(),2) + pnts[j-1].z()*pnts[j].z() + pow(pnts[j].z(),2));
		inertias[1] += (pnts[j-1].x()*pnts[j].z() - pnts[j].x()*pnts[j-1].z()) \
				* (pow(pnts[j-1].x(),2) + pnts[j-1].x()*pnts[j].x() + pow(pnts[j].x(),2));
	}
	inertias[0] += (pnts[j-1].x()*pnts[0].z() - pnts[0].x()*pnts[j-1].z()) \
			* (pow(pnts[j-1].z(),2) + pnts[j-1].z()*pnts[0].z() + pow(pnts[0].z(),2));
	inertias[1] += (pnts[j-1].x()*pnts[0].z() - pnts[0].x()*pnts[j-1].z()) \
			* (pow(pnts[j-1].x(),2) + pnts[j-1].x()*pnts[0].x() + pow(pnts[0].x(),2));
	inertias[0] /= (12*area);
	inertias[1] /= (12*area);

	inertiasCG[0] = ( inertias[0] - pow(centroid.z(),2) ) * area;
	inertiasCG[1] = ( inertias[1] - pow(centroid.x(),2) ) * area;
	inertiasCG[2] = inertiasCG[0] + inertiasCG[1];

	return inertiasCG;
}

vector<double> Xsec_surf::calculate_solid_inertias_in_plane(int ixs, int plane, float mat[4][4])
{
	vector<double> inertias(2,0), inertiasCG(3,0);
	int j;
	double area = get_xsec_plane_area(ixs, plane, mat);
	vec2d tmpPnts, cg = get_xsec_centroid_in_plane(ixs, plane, mat);
	vector<vec2d> pnts;

	switch(plane)
	{
	case XY_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( pnts_xsecs(ixs,i).transform(mat).x(), pnts_xsecs(ixs,i).transform(mat).y() );
			pnts.push_back(tmpPnts);
		}
		break;
	case XZ_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( pnts_xsecs(ixs,i).transform(mat).x(), pnts_xsecs(ixs,i).transform(mat).z() );
			pnts.push_back(tmpPnts);
		}
		break;
	case YZ_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( pnts_xsecs(ixs,i).transform(mat).y(), pnts_xsecs(ixs,i).transform(mat).z() );
			pnts.push_back(tmpPnts);
		}
		break;
	}

	for ( j = 1; j < num_pnts; j++ )
	{
		inertias[0] += ( pnts[j-1].x()*pnts[j].y() - pnts[j].x()*pnts[j-1].y() ) \
				* ( pow( pnts[j-1].y(), 2) + pnts[j-1].y()*pnts[j].y() + pow( pnts[j].y(), 2) );
		inertias[1] += ( pnts[j-1].x()*pnts[j].y() - pnts[j].x()*pnts[j-1].y() ) \
				* ( pow( pnts[j-1].x(), 2) + pnts[j-1].x()*pnts[j].x() + pow( pnts[j].x(), 2) );
	}
	inertias[0] += (pnts[j-1].x()*pnts[0].y() - pnts[0].x()*pnts[j-1].y()) \
			* (pow(pnts[j-1].y(),2) + pnts[j-1].y()*pnts[0].y() + pow(pnts[0].y(),2));
	inertias[1] += (pnts[j-1].x()*pnts[0].y() - pnts[0].x()*pnts[j-1].y()) \
			* (pow(pnts[j-1].x(),2) + pnts[j-1].x()*pnts[0].x() + pow(pnts[0].x(),2));
	inertias[0] /= (12*area);
	inertias[1] /= (12*area);

	inertiasCG[0] = ( inertias[0] - pow(cg.y(),2) ) * abs(area);
	inertiasCG[1] = ( inertias[1] - pow(cg.x(),2) ) * abs(area);
	inertiasCG[2] = inertiasCG[0] + inertiasCG[1];

	return inertiasCG;
}

vector<double> Xsec_surf::calculate_refl_solid_inertias_in_plane(int ixs, int plane, float refl_mat[4][4])
{
	vector<double> inertias(2,0), inertiasCG(3,0);
	int j;
	double area = get_refl_xsec_plane_area(ixs, plane, refl_mat);
	vec2d tmpPnts, cg = get_refl_xsec_centroid_in_plane(ixs, plane, refl_mat);
	vector<vec2d> pnts;

	switch(plane)
	{
	case XY_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( refl_pnts_xsecs(ixs,i).transform(refl_mat).x(), refl_pnts_xsecs(ixs,i).transform(refl_mat).y() );
			pnts.push_back(tmpPnts);
		}
		break;
	case XZ_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( refl_pnts_xsecs(ixs,i).transform(refl_mat).x(), refl_pnts_xsecs(ixs,i).transform(refl_mat).z() );
			pnts.push_back(tmpPnts);
		}
		break;
	case YZ_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( refl_pnts_xsecs(ixs,i).transform(refl_mat).y(), refl_pnts_xsecs(ixs,i).transform(refl_mat).z() );
			pnts.push_back(tmpPnts);
		}
		break;
	}

	for ( j = 1; j < num_pnts; j++ )
	{
		inertias[0] += ( pnts[j-1].x()*pnts[j].y() - pnts[j].x()*pnts[j-1].y() ) \
				* ( pow( pnts[j-1].y(), 2) + pnts[j-1].y()*pnts[j].y() + pow( pnts[j].y(), 2) );
		inertias[1] += ( pnts[j-1].x()*pnts[j].y() - pnts[j].x()*pnts[j-1].y() ) \
				* ( pow( pnts[j-1].x(), 2) + pnts[j-1].x()*pnts[j].x() + pow( pnts[j].x(), 2) );
	}
	inertias[0] += (pnts[j-1].x()*pnts[0].y() - pnts[0].x()*pnts[j-1].y()) \
			* (pow(pnts[j-1].y(),2) + pnts[j-1].y()*pnts[0].y() + pow(pnts[0].y(),2));
	inertias[1] += (pnts[j-1].x()*pnts[0].y() - pnts[0].x()*pnts[j-1].y()) \
			* (pow(pnts[j-1].x(),2) + pnts[j-1].x()*pnts[0].x() + pow(pnts[0].x(),2));
	inertias[0] /= (12*area);
	inertias[1] /= (12*area);

	inertiasCG[0] = ( inertias[0] - pow(cg.y(),2) ) * area;
	inertiasCG[1] = ( inertias[1] - pow(cg.x(),2) ) * area;
	inertiasCG[2] = inertiasCG[0] + inertiasCG[1];

	return inertiasCG;
}
