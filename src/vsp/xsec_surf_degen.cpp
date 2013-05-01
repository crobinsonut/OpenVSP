
#include "xsec_surf.h"
#include "vector_util.h"
#include "geom.h"

#include "bezier_surf.h"

#include "tMesh.h"
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
