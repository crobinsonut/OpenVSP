
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
