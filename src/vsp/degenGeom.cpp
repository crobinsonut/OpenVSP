#include "degenGeom.h"
#include "geom.h"
#include <cmath>

//==== Get area normal vector for PLANAR cross section ====//
vec3d DegenGeom::get_area_normal( int ixs, const array_2d<vec3d> &pntsarr)
{
	vec3d areaNormal, refVec1, refVec2;
	int halfPnts = (num_pnts + 1) / 2;

	refVec1 = pntsarr(ixs, halfPnts - 1) - pntsarr(ixs, 0);
	refVec2 = pntsarr(ixs, 1) - pntsarr(ixs, num_pnts - 2);

	areaNormal = cross(refVec2, refVec1);
	areaNormal.normalize();

	return areaNormal;
}

//==== Get Cross Section Area ====//
// JBB: This only works for PLANAR cross sections.
double DegenGeom::get_xsec_area( int ixs, const array_2d<vec3d> &pntsarr )
{
	double area = 0;
	int j;
	vec3d xAxis(1,0,0), yAxis(0,1,0), zAxis(0,0,1), areaNormal = get_area_normal(ixs, pntsarr);
	vector<vec3d> pnts;

	// Load cross section points
	for( int i = 0; i < num_pnts; i++ )
	{
		pnts.push_back( pntsarr(ixs,i) );
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
double DegenGeom::get_xsec_plane_area( int ixs, int plane, float mat[4][4], const array_2d<vec3d> &pntsarr )
{
	double area = 0;
	int j;
	vector<vec3d> pnts;

	for ( int i = 0; i < num_pnts; i++ )
	{
		pnts.push_back( pntsarr(ixs,i).transform(mat) );
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
vec3d DegenGeom::get_xsec_centroid( int ixs, const array_2d<vec3d> &pntsarr )
{
	double xcg = 0, zcg = 0, area = get_xsec_area(ixs, pntsarr);
	int j;
	vec3d xAxis(1,0,0), zAxis(0,0,1), areaNormal = get_area_normal(ixs, pntsarr);

	vector<vec3d> pnts;
	// Load cross section points
	for( int i = 0; i < num_pnts; i++ )
	{
		pnts.push_back( pntsarr(ixs,i) );
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

vec2d DegenGeom::get_xsec_centroid_in_plane(int ixs, int plane, float mat[4][4], const array_2d<vec3d> &pntsarr)
{
	double xcg = 0, ycg = 0, area = get_xsec_plane_area(ixs, plane, mat, pntsarr);
	int j;
	vector<vec2d> pnts;
	vec2d tmpPnts;

	switch(plane)
	{
	case XY_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( pntsarr(ixs,i).transform(mat).x(), pntsarr(ixs,i).transform(mat).y() );
			pnts.push_back(tmpPnts);
		}
		break;
	case XZ_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( pntsarr(ixs,i).transform(mat).x(), pntsarr(ixs,i).transform(mat).z() );
			pnts.push_back(tmpPnts);
		}
		break;
	case YZ_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( pntsarr(ixs,i).transform(mat).y(), pntsarr(ixs,i).transform(mat).z() );
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


void DegenGeom::createDegenSurface(int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr, bool refl)
{
	int nLow = 0, nHigh = num_xsecs;

	vector< vector<vec3d> > xSurfMat = degenSurface.x;
	vector< vector<vec3d> > nSurfMat = degenSurface.nvec;

	vector<vec3d> xVec( num_pnts );
	vector<vec3d> nVec( num_pnts-1 );

	if ( getType() == DegenGeom::SURFACE_TYPE )
	{
		if ( getParentGeom()->getTypeStr() == "wing" || getParentGeom()->getTypeStr() == "prop" )
		{
			nLow  = 1;
			nHigh = num_xsecs - 1;
		}
	}

	for ( int i = nLow; i < nHigh-1; i++ )
	{
		for ( int j = 0; j < num_pnts-1; j++ )
		{
			vec3d sVec1 = pntsarr(i+1,j).transform(mat) - pntsarr(i,j).transform(mat);
			vec3d sVec2 = pntsarr(i,j+1).transform(mat) - pntsarr(i,j).transform(mat);

			if(!refl)
				nVec[j]     = cross(sVec1,sVec2);
			else
				nVec[j]     = cross(sVec2,sVec1);

			nVec[j].normalize();
		}

		nSurfMat.push_back(nVec);
	}

	for ( int i = nLow; i < nHigh; i++ )
	{
		for ( int j = 0; j < num_pnts; j++ )
		{
			xVec[j] = pntsarr(i,j).transform(mat);
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
}

void DegenGeom::createSurfDegenPlate(int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr)
{
	int nLow = 0, nHigh = num_xsecs;
	int startPnt = 0;

	if ( getParentGeom()->getTypeStr() == "wing" || getParentGeom()->getTypeStr() == "prop" )
	{
		// Keep only airfoil sections, discard endcap close-out lines
		nLow  = 1;
		nHigh = num_xsecs - 1;
	}

	degenPlates.resize(1);
	createDegenPlate(degenPlates[0], sym_code_in, mat, pntsarr, nLow, nHigh, startPnt);
}

void DegenGeom::createBodyDegenPlate(int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr)
{
	int nLow = 0, nHigh = num_xsecs;

	degenPlates.resize(2);

	int startPnt = 0;
	createDegenPlate(degenPlates[0], sym_code_in, mat, pntsarr, nLow, nHigh, startPnt);

	startPnt = (num_pnts - 1) / 4;
	createDegenPlate(degenPlates[1], sym_code_in, mat, pntsarr, nLow, nHigh, startPnt);
}

void DegenGeom::createDegenPlate(DegenPlate &degenPlate, int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr, int nLow, int nHigh, int startPnt)
{
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

	vec3d  lePnt, tePnt, topPnt, botPnt, chordVec, camberPnt, chordPnt, nPlate;

	for ( int i = nLow; i < nHigh; i++ )
	{
		lePnt = pntsarr(i,startPnt);
		tePnt = pntsarr(i, startPnt+platePnts-1);

		// Set first point (trailing edge)
		xVec[0]  = lePnt.transform(mat);
		nCVec[0] = vec3d(0,0,0); // on camber line
		tVec[0]  = 0;
		zVec[0]  = 0;

		// Set last point (leading edge)
		xVec[platePnts-1]  = tePnt.transform(mat);
		nCVec[platePnts-1] = vec3d(0,0,0); // on camber line
		tVec[platePnts-1]  = 0;
		zVec[platePnts-1]  = 0;

		//== Compute plate normal ==//
		// rotated chord vector
		vec3d rcv = xVec[platePnts-1] - xVec[0];
		// rotated area normal vector
		vec3d anv = get_area_normal(i, pntsarr).transform(mat) - vec3d(0,0,0).transform(mat);
		// plate normal vector
		nPlate = cross(rcv,anv);
		nPlate.normalize();
		degenPlate.nPlate.push_back( nPlate );

		// normalized, unrotated chord vector (te->le)
		chordVec = tePnt - lePnt;
		chordVec.normalize();

		for ( int j = 1, k = startPnt+num_pnts-2; j < platePnts-1; j++, k-- )
		{
			topPnt = pntsarr(i,startPnt+j);
			botPnt = pntsarr(i, k % (num_pnts-1) );

			camberPnt = (topPnt + botPnt) / 2;

			nCVec[j] = topPnt.transform(mat)-botPnt.transform(mat);  // vector from bottom point to top point.
			nCVec[j].normalize();

			tVec[j]  = dist(topPnt,botPnt);

			// Try finding least-squares minimum distance point.
			double s, t;
			bool intflag = line_line_intersect( botPnt, topPnt, lePnt, tePnt, &s, &t );
			if( intflag )
			{
				chordPnt = lePnt + (tePnt - lePnt) * t;
			}
			else  // If it doesn't work, project.  Projection isn't quite right, but should be close.  This should never happen anyway.
			{
				chordPnt  = lePnt + chordVec * dot(camberPnt-lePnt,chordVec);
			}

			xVec[j]  = chordPnt.transform(mat);
			zVec[j]  = dist(camberPnt,chordPnt);
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
}

void DegenGeom::createSurfDegenStick(int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr)
{
	int nLow = 0, nHigh = num_xsecs;
	int platePnts = (num_pnts + 1) / 2;
	vec3d chordVec, camberPnt, prevCamberPnt;

	if ( getParentGeom()->getTypeStr() == "wing" || getParentGeom()->getTypeStr() == "prop" )
	{
		// Keep only airfoil sections, discard endcap close-out lines
		nLow  = 1;
		nHigh = num_xsecs - 1;
	}

	degenSticks.resize(1);

	int startPnt = 0;
	createDegenStick(degenSticks[0], sym_code_in, mat, pntsarr, nLow, nHigh, startPnt);
	createDegenStickSweep(degenSticks[0], sym_code_in, mat, pntsarr, nLow, nHigh, startPnt);

}

void DegenGeom::createDegenStickSweep(DegenStick &degenStick, int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr, int nLow, int nHigh, int startPnt)
{
	int platePnts = (num_pnts + 1) / 2;

	// Calculate sweep angle
	for( int i = nLow; i < nHigh - 1; i++ )
	{
		vec3d  cvCurrent, cvNext, qcCurrent, qcNext;
		double chordCurrent, chordNext;

		// Get current section chord vector
		cvCurrent = pntsarr(i, platePnts-1) - pntsarr(i,0);
		chordCurrent = cvCurrent.mag();
		cvCurrent.normalize();
		// Get current section quarter chord point
		qcCurrent = (pntsarr(i,0) + cvCurrent*0.75*chordCurrent).transform(mat);

		// Get next section chord vector
		cvNext = pntsarr(i+1, platePnts-1) - pntsarr(i+1,0);
		chordNext = cvNext.mag();
		cvNext.normalize();
		// Get next section quarter chord point
		qcNext = (pntsarr(i+1,0) + cvNext*0.75*chordNext).transform(mat);

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
}

void DegenGeom::createBodyDegenStick(int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr)
{
	int nLow = 0, nHigh = num_xsecs;

	degenSticks.resize(2);

	int startPnt = 0;
	createDegenStick(degenSticks[0], sym_code_in, mat, pntsarr, nLow, nHigh, startPnt);

	startPnt = (num_pnts - 1) / 4;
	createDegenStick(degenSticks[1], sym_code_in, mat, pntsarr, nLow, nHigh, startPnt);
}

void DegenGeom::createDegenStick(DegenStick &degenStick, int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr, int nLow, int nHigh, int startPnt)
{
	int platePnts = (num_pnts + 1) / 2;
	vec3d chordVec, camberPnt, prevCamberPnt;

	for ( int i = nLow; i < nHigh; i++ )
	{
		vec3d topPnt, botPnt;
		double tempThickness = 0, perimTop = 0, perimBot = 0;
		int    maxThickIdx[2] = {0,0};

		// normalized, unrotated chord vector (te->le)
		chordVec = pntsarr(i, startPnt+platePnts-1) - pntsarr(i, startPnt);
		chordVec.normalize();

		degenStick.xle.push_back( pntsarr( i, startPnt+platePnts-1 ).transform(mat) );
		degenStick.xte.push_back( pntsarr( i, startPnt ).transform(mat) );
		degenStick.chord.push_back( dist(pntsarr(i, startPnt+platePnts-1), pntsarr(i, startPnt)) );
		degenStick.u.push_back( uArray[i] );
		degenStick.Ishell.push_back(calculate_shell_inertias_in_plane(i,YZ_PLANE, mat, pntsarr));
		degenStick.Isolid.push_back(calculate_solid_inertias_in_plane(i,YZ_PLANE, mat, pntsarr));
		degenStick.xcgSolid.push_back( get_xsec_centroid(i, pntsarr).transform(mat) );
		degenStick.xcgShell.push_back( get_xsec_shellCG(i, pntsarr).transform(mat)  );
		degenStick.area.push_back( get_xsec_area(i, pntsarr) );

		vec3d areaNormal = get_area_normal(i, pntsarr).transform(mat) - vec3d(0,0,0).transform(mat);
		areaNormal.normalize();
		degenStick.areaNormal.push_back( areaNormal );

		for ( int j = 1, k = startPnt+num_pnts-2; j < platePnts-1; j++, k-- )
		{
			topPnt = pntsarr(i,startPnt+j);
			botPnt = pntsarr(i, k % (num_pnts-1) );

			camberPnt = ( topPnt + botPnt ) / 2;

			if( dist(topPnt, botPnt) > tempThickness)
			{
				tempThickness  = dist(topPnt, botPnt);
				maxThickIdx[0] = j+startPnt;
				maxThickIdx[1] = k % (num_pnts-1);
			}
			perimTop += dist( pntsarr(i,startPnt+j), pntsarr(i,startPnt+j-1) );
			perimBot += dist( pntsarr(i, k % (num_pnts-1) ), pntsarr(i, k % (num_pnts-1)+1) );

			prevCamberPnt = camberPnt;
		}

		camberPnt = ( pntsarr(i,maxThickIdx[0]) + pntsarr(i,maxThickIdx[1]) ) / 2;
		degenStick.tLoc.push_back( 1 - (dot(camberPnt-pntsarr(i,startPnt),chordVec) / degenStick.chord.back()) );
		degenStick.toc.push_back( tempThickness / degenStick.chord.back() );

		perimTop += dist( pntsarr(i, startPnt+platePnts-1), pntsarr(i, startPnt+platePnts-2) );
		perimBot += dist( pntsarr(i, startPnt+platePnts), pntsarr(i, startPnt+platePnts-1) );
		degenStick.perimTop.push_back( perimTop );
		degenStick.perimBot.push_back( perimBot );
	}

	// Calculate sweep angle
	for( int i = nLow; i < nHigh; i++ )
	{
		degenStick.sweep.push_back( NAN );
	}
}

vec3d DegenGeom::get_xsec_shellCG( int ixs, const array_2d<vec3d> &pntsarr )
{
	double xcg = 0, zcg = 0, area = get_xsec_area(ixs, pntsarr);
	int j;
	vec3d xAxis(1,0,0), zAxis(0,0,1), areaNormal = get_area_normal(ixs, pntsarr);

	vector<vec3d> pnts;
	// Load cross section points
	for( int i = 0; i < num_pnts; i++ )
	{
		pnts.push_back( pntsarr(ixs,i) );
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
vector<double> DegenGeom::calculate_shell_inertias(int ixs, const array_2d<vec3d> &pntsarr )
{
	int j, platePnts = (num_pnts + 1) / 2;
	vector<double> inertias(6,0);
	vec3d xAxis(1,0,0), zAxis(0,0,1);
	vec3d areaNormal = get_area_normal(ixs, pntsarr), CG = get_xsec_centroid(ixs, pntsarr);

	vector<vec3d> pnts;
	// Load cross section points
	for( int i = 0; i < num_pnts; i++ )
	{
		pnts.push_back( pntsarr(ixs,i) );
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

vector<double> DegenGeom::calculate_shell_inertias_in_plane(int ixs, int plane, float mat[4][4], const array_2d<vec3d> &pntsarr)
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
			tmpPnts.set_xy( pntsarr(ixs,i).transform(mat).x(), pntsarr(ixs,i).transform(mat).y() );
			pnts.push_back(tmpPnts);
		}
		break;
	case XZ_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( pntsarr(ixs,i).transform(mat).x(), pntsarr(ixs,i).transform(mat).z() );
			pnts.push_back(tmpPnts);
		}
		break;
	case YZ_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( pntsarr(ixs,i).transform(mat).y(), pntsarr(ixs,i).transform(mat).z() );
			pnts.push_back(tmpPnts);
		}
		break;
	}

	vec2d segVec, chordVec = pnts[0] - pnts[platePnts-1];
	vec2d distToCG, cg, CG = get_xsec_centroid_in_plane(ixs, plane, mat, pntsarr);
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

vector<double> DegenGeom::calculate_solid_inertias( int ixs, const array_2d<vec3d> &pntsarr )
{
	vector<double> inertias(2,0), inertiasCG(3,0);
	int j;
	double area = get_xsec_area(ixs, pntsarr);
	vec3d xAxis(1,0,0), zAxis(0,0,1);
	vec3d areaNormal = get_area_normal(ixs, pntsarr);
	vec3d centroid   = get_xsec_centroid(ixs, pntsarr);

	vector<vec3d> pnts;
	// Load cross section points
	for( int i = 0; i < num_pnts; i++ )
	{
		pnts.push_back( pntsarr(ixs,i) );
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

vector<double> DegenGeom::calculate_solid_inertias_in_plane(int ixs, int plane, float mat[4][4], const array_2d<vec3d> &pntsarr)
{
	vector<double> inertias(2,0), inertiasCG(3,0);
	int j;
	double area = get_xsec_plane_area(ixs, plane, mat, pntsarr);
	vec2d tmpPnts, cg = get_xsec_centroid_in_plane(ixs, plane, mat, pntsarr);
	vector<vec2d> pnts;

	switch(plane)
	{
	case XY_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( pntsarr(ixs,i).transform(mat).x(), pntsarr(ixs,i).transform(mat).y() );
			pnts.push_back(tmpPnts);
		}
		break;
	case XZ_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( pntsarr(ixs,i).transform(mat).x(), pntsarr(ixs,i).transform(mat).z() );
			pnts.push_back(tmpPnts);
		}
		break;
	case YZ_PLANE:
		for ( int i = 0; i < num_pnts; i++ )
		{
			tmpPnts.set_xy( pntsarr(ixs,i).transform(mat).y(), pntsarr(ixs,i).transform(mat).z() );
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

const char* DegenGeom::makeCsvFmt( int n )
{
	char fmt[10];
	sprintf( fmt, "%%.%de", DBL_DIG + 3 );

	string fmtstring = "";
	for( int i=0; i < n; ++i)
	{
		fmtstring.append( fmt );
		if(i < n-1 )
			fmtstring.append(", ");
		else
			fmtstring.append("\n");
	}
	return fmtstring.c_str();
}

void DegenGeom::write_degenGeomSurfCsv_file(FILE* file_id, int nxsecs)
{
	fprintf(file_id, "# DegenGeom Type,nXsecs, nPnts/Xsec\n");
	fprintf(file_id, "SURFACE_NODE,%d,%d\n", nxsecs, num_pnts);
	fprintf(file_id, "# x,y,z,u,w\n");

	for ( int i = 0; i < nxsecs; i++ )
	{
		for ( int j = 0; j < num_pnts; j++ )
		{
			fprintf(file_id, makeCsvFmt(5),			\
					degenSurface.x[i][j].x(),		\
					degenSurface.x[i][j].y(),		\
					degenSurface.x[i][j].z(),		\
					degenSurface.u[i],				\
					degenSurface.w[j]				);
		}
	}

	fprintf(file_id, "SURFACE_FACE,%d,%d\n", nxsecs-1, num_pnts-1);
	fprintf(file_id, "# xn,yn,zn\n");

	for ( int i = 0; i < nxsecs-1; i++ )
	{
		for ( int j = 0; j < num_pnts-1; j++ )
		{
			fprintf(file_id, makeCsvFmt(3),			\
					degenSurface.nvec[i][j].x(),	\
					degenSurface.nvec[i][j].y(),	\
					degenSurface.nvec[i][j].z()		);
		}
	}


}

void DegenGeom::write_degenGeomPlateCsv_file(FILE* file_id, int nxsecs, DegenPlate &degenPlate)
{
	fprintf(file_id, "# DegenGeom Type,nXsecs,nPnts/Xsec\n");
	fprintf(file_id, "PLATE,%d,%d\n", nxsecs, (num_pnts+1)/2);
	fprintf(file_id,"# xn,yn,zn\n");
	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf(file_id, makeCsvFmt(3), degenPlate.nPlate[i].x(), \
				degenPlate.nPlate[i].y(), \
				degenPlate.nPlate[i].z()  );
	}

	fprintf(file_id, "# x,y,z,zCamb,t,nCambX,nCambY,nCambZ,u,wTop,wBot\n");
	for ( int i = 0; i < nxsecs; i++ )
	{
		for ( int j = 0; j < (num_pnts+1)/2; j++ )
		{
			fprintf(file_id, makeCsvFmt(11),	\
					degenPlate.x[i][j].x(),				\
					degenPlate.x[i][j].y(),				\
					degenPlate.x[i][j].z(),				\
					degenPlate.zcamber[i][j],			\
					degenPlate.t[i][j],					\
					degenPlate.nCamber[i][j].x(),		\
					degenPlate.nCamber[i][j].y(),		\
					degenPlate.nCamber[i][j].z(),		\
					degenPlate.u[i],					\
					degenPlate.wTop[j],					\
					degenPlate.wBot[j]					);
		}
	}
}

void DegenGeom::write_degenGeomStickCsv_file(FILE* file_id, int nxsecs, DegenStick &degenStick)
{

	fprintf(file_id, "# DegenGeom Type, nXsecs\nSTICK, %d\n# xle,yle,zle,xte,yte,zte,xcg_solid,ycg_solid,zcg_solid,"
			"xcg_shell,ycg_shell,zcg_shell,toc,tLoc,chord,sweep,Ixx_shell_A,Ixx_shell_B,Izz_shell_A,"
			"Izz_shell_B,J_shell_A,J_shell_B,Ixx_solid,Izz_solid,J_solid,area,areaNormalX,"
			"areaNormalY,areaNormalZ,perimTop,perimBot,u\n", nxsecs);

	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf(file_id, makeCsvFmt(32),	\
				degenStick.xle[i].x(),					\
				degenStick.xle[i].y(),					\
				degenStick.xle[i].z(),					\
				degenStick.xte[i].x(),					\
				degenStick.xte[i].y(),					\
				degenStick.xte[i].z(),					\
				degenStick.xcgSolid[i].x(),				\
				degenStick.xcgSolid[i].y(),				\
				degenStick.xcgSolid[i].z(),				\
				degenStick.xcgShell[i].x(),				\
				degenStick.xcgShell[i].y(),				\
				degenStick.xcgShell[i].z(),				\
				degenStick.toc[i],						\
				degenStick.tLoc[i],						\
				degenStick.chord[i],					\
				degenStick.sweep[i],					\
				degenStick.Ishell[i][0],				\
				degenStick.Ishell[i][1],				\
				degenStick.Ishell[i][2],				\
				degenStick.Ishell[i][3],				\
				degenStick.Ishell[i][4],				\
				degenStick.Ishell[i][5],				\
				degenStick.Isolid[i][0],				\
				degenStick.Isolid[i][1],				\
				degenStick.Isolid[i][2],				\
				degenStick.area[i],						\
				degenStick.areaNormal[i].x(),			\
				degenStick.areaNormal[i].y(),			\
				degenStick.areaNormal[i].z(),			\
				degenStick.perimTop[i],					\
				degenStick.perimBot[i],					\
				degenStick.u[i]							);
	}
}

void DegenGeom::write_degenGeomPointCsv_file(FILE* file_id, int nxsecs)
{
	fprintf(file_id, "# DegenGeom Type\n");
	fprintf(file_id, "POINT\n");
	fprintf(file_id, "# vol,volWet,area,areaWet,IxxShell,IyyShell,IzzShell,IxyShell,");
	fprintf(file_id, "IxzShell,IyzShell,IxxSolid,IyySolid,IzzSolid,IxySolid,IxzSolid,");
	fprintf(file_id, "IyzSolid,xcgShell,ycgShell,zcgShell,xcgSolid,ycgSolid,zcgSolid\n"); // Need newline below.  Omitting for diff-consistency.
	fprintf(file_id, makeCsvFmt(22),\
			degenPoint.vol[0],			\
			degenPoint.volWet[0],		\
			degenPoint.area[0],			\
			degenPoint.areaWet[0],		\
			degenPoint.Ishell[0][0],	\
			degenPoint.Ishell[0][1],	\
			degenPoint.Ishell[0][2],	\
			degenPoint.Ishell[0][3],	\
			degenPoint.Ishell[0][4],	\
			degenPoint.Ishell[0][5],	\
			degenPoint.Isolid[0][0],	\
			degenPoint.Isolid[0][1],	\
			degenPoint.Isolid[0][2],	\
			degenPoint.Isolid[0][3],	\
			degenPoint.Isolid[0][4],	\
			degenPoint.Isolid[0][5],	\
			degenPoint.xcgShell[0].x(),	\
			degenPoint.xcgShell[0].y(),	\
			degenPoint.xcgShell[0].z(),	\
			degenPoint.xcgSolid[0].x(),	\
			degenPoint.xcgSolid[0].y(),	\
			degenPoint.xcgSolid[0].z()	);
}

void DegenGeom::write_degenGeomPropCsv_file(FILE* file_id)
{
	char fmtstr[255];
	strcat(fmtstr, "%d, ");
	strcat(fmtstr, makeCsvFmt(7) );
	if ( parentGeom->getTypeStr() == "prop" )
	{
		fprintf(file_id, "# DegenGeom Type\n");
		fprintf(file_id, "PROP\n");
		fprintf(file_id, "# Num Blades, Diameter, xLoc, yLoc, zLoc, nRotX, nRotY, nRotZ\n");
		fprintf(file_id, fmtstr,\
				degenProp.nblade,	\
				degenProp.d,		\
				degenProp.x.x(),	\
				degenProp.x.y(),	\
				degenProp.x.z(),	\
				degenProp.nvec.x(),	\
				degenProp.nvec.y(),	\
				degenProp.nvec.z()	);
	}
}

void DegenGeom::write_degenGeomCsv_file(FILE* file_id)
{
	int nxsecs = num_xsecs;

	if( type == SURFACE_TYPE )
		fprintf(file_id, "\nLIFTING_SURFACE,%s\n", name.c_str() );
	else
		fprintf(file_id, "\nBODY,%s\n", name.c_str() );

	write_degenGeomPropCsv_file(file_id);

	if ( parentGeom->getTypeStr() == "wing" || parentGeom->getTypeStr() == "prop" )
		nxsecs -= 2;

	write_degenGeomSurfCsv_file( file_id, nxsecs);

	write_degenGeomPlateCsv_file( file_id, nxsecs, degenPlates[0]);

	if ( type == DegenGeom::BODY_TYPE )
		write_degenGeomPlateCsv_file( file_id, nxsecs, degenPlates[1]);

	write_degenGeomStickCsv_file(file_id, nxsecs, degenSticks[0]);

	if ( type == DegenGeom::BODY_TYPE )
		write_degenGeomStickCsv_file(file_id, nxsecs, degenSticks[1]);

	write_degenGeomPointCsv_file(file_id, nxsecs);
}

void DegenGeom::write_degenGeomSurfM_file(FILE* file_id, int nxsecs)
{
	//============================= DegenSurf =============================//
	fprintf(file_id, "\ndegenGeom(end).surf.x = [");
	for ( int i = 0; i < nxsecs; i++ )
	{
		for ( int j = 0; j < num_pnts-1; j++ )
			fprintf( file_id, "%f, ", degenSurface.x[i][j].x() );
		if ( i < nxsecs-1 )
			fprintf( file_id, "%f;\n", degenSurface.x[i][num_pnts-1].x() );
		else
			fprintf( file_id, "%f];\n", degenSurface.x[i][num_pnts-1].x() );
	}

	fprintf(file_id, "\ndegenGeom(end).surf.y = [");
	for ( int i = 0; i < nxsecs; i++ )
	{
		for ( int j = 0; j < num_pnts-1; j++ )
			fprintf( file_id, "%f, ", degenSurface.x[i][j].y() );
		if ( i < nxsecs-1 )
			fprintf( file_id, "%f;\n", degenSurface.x[i][num_pnts-1].y() );
		else
			fprintf( file_id, "%f];\n", degenSurface.x[i][num_pnts-1].y() );
	}

	fprintf(file_id, "\ndegenGeom(end).surf.z = [");
	for ( int i = 0; i < nxsecs; i++ )
	{
		for ( int j = 0; j < num_pnts-1; j++ )
			fprintf( file_id, "%f, ", degenSurface.x[i][j].z() );
		if ( i < nxsecs-1 )
			fprintf( file_id, "%f;\n", degenSurface.x[i][num_pnts-1].z() );
		else
			fprintf( file_id, "%f];\n", degenSurface.x[i][num_pnts-1].z() );
	}

	// Normal vectors
	fprintf(file_id, "\ndegenGeom(end).surf.nx = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		for ( int j = 0; j < num_pnts-2; j++ )
			fprintf( file_id, "%f, ", degenSurface.nvec[i][j].x() );
		if ( i < nxsecs-2 )
			fprintf( file_id, "%f;\n", degenSurface.nvec[i][num_pnts-2].x() );
		else
			fprintf( file_id, "%f];\n", degenSurface.nvec[i][num_pnts-2].x() );
	}

	fprintf(file_id, "\ndegenGeom(end).surf.ny = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		for ( int j = 0; j < num_pnts-2; j++ )
			fprintf( file_id, "%f, ", degenSurface.nvec[i][j].y() );
		if ( i < nxsecs-2 )
			fprintf( file_id, "%f;\n", degenSurface.nvec[i][num_pnts-2].y() );
		else
			fprintf( file_id, "%f];\n", degenSurface.nvec[i][num_pnts-2].y() );
	}

	fprintf(file_id, "\ndegenGeom(end).surf.nz = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		for ( int j = 0; j < num_pnts-2; j++ )
			fprintf( file_id, "%f, ", degenSurface.nvec[i][j].z() );
		if ( i < nxsecs-2 )
			fprintf( file_id, "%f;\n", degenSurface.nvec[i][num_pnts-2].z() );
		else
			fprintf( file_id, "%f];\n", degenSurface.nvec[i][num_pnts-2].z() );
	}

	fprintf(file_id, "\ndegenGeom(end).surf.u = [");
	for ( int i = 0; i < nxsecs - 1; i++ )
	{
		fprintf(file_id, "%f, ", degenSurface.u[i]);
	}
	fprintf(file_id, "%f];", degenSurface.u[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).surf.w = [");
	for ( int j = 0; j < num_pnts-1; j++ )
	{
		fprintf(file_id, "%f, ", degenSurface.w[j]);
	}
	fprintf(file_id, "%f];", degenSurface.w[num_pnts-1]);
}

void DegenGeom::write_degenGeomPlateM_file(FILE* file_id, int nxsecs, DegenPlate &degenPlate, int iplate)
{
	//============================= DegenPlate =============================//
	fprintf(file_id, "\ndegenGeom(end).plate(%d).nx = [", iplate);
	for ( int i = 0; i < nxsecs - 1; i++ )
	{
		fprintf(file_id, "%f, ", degenPlate.nPlate[i].x());
	}
	fprintf(file_id, "%f];", degenPlate.nPlate[nxsecs-1].x());

	fprintf(file_id, "\ndegenGeom(end).plate(%d).ny = [", iplate);
	for ( int i = 0; i < nxsecs - 1; i++ )
	{
		fprintf(file_id, "%f, ", degenPlate.nPlate[i].y());
	}
	fprintf(file_id, "%f];", degenPlate.nPlate[nxsecs-1].y());

	fprintf(file_id, "\ndegenGeom(end).plate(%d).nz = [", iplate);
	for ( int i = 0; i < nxsecs - 1; i++ )
	{
		fprintf(file_id, "%f, ", degenPlate.nPlate[i].z());
	}
	fprintf(file_id, "%f];", degenPlate.nPlate[nxsecs-1].z());

	fprintf(file_id, "\ndegenGeom(end).plate(%d).x = [", iplate);
	for ( int i = 0; i < nxsecs; i++ )
	{
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
			fprintf( file_id, "%f, ", degenPlate.x[i][j].x() );
		if ( i < nxsecs-1 )
			fprintf( file_id, "%f;\n", degenPlate.x[i][(num_pnts+1)/2-1].x() );
		else
			fprintf( file_id, "%f];\n", degenPlate.x[i][(num_pnts+1)/2-1].x() );
	}

	fprintf(file_id, "\ndegenGeom(end).plate(%d).y = [", iplate);
	for ( int i = 0; i < nxsecs; i++ )
	{
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
			fprintf( file_id, "%f, ", degenPlate.x[i][j].y() );
		if ( i < nxsecs-1 )
			fprintf( file_id, "%f;\n", degenPlate.x[i][(num_pnts+1)/2-1].y() );
		else
			fprintf( file_id, "%f];\n", degenPlate.x[i][(num_pnts+1)/2-1].y() );
	}

	fprintf(file_id, "\ndegenGeom(end).plate(%d).z = [", iplate);
	for ( int i = 0; i < nxsecs; i++ )
	{
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
			fprintf( file_id, "%f, ", degenPlate.x[i][j].z() );
		if ( i < nxsecs-1 )
			fprintf( file_id, "%f;\n", degenPlate.x[i][(num_pnts+1)/2-1].z() );
		else
			fprintf( file_id, "%f];\n", degenPlate.x[i][(num_pnts+1)/2-1].z() );
	}

	fprintf(file_id, "\ndegenGeom(end).plate(%d).zCamber = [", iplate);
	for ( int i = 0; i < nxsecs; i++ )
	{
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
			fprintf( file_id, "%f, ", degenPlate.zcamber[i][j] );
		if ( i < nxsecs-1 )
			fprintf( file_id, "%f;\n", degenPlate.zcamber[i][(num_pnts+1)/2-1] );
		else
			fprintf( file_id, "%f];\n", degenPlate.zcamber[i][(num_pnts+1)/2-1] );
	}

	fprintf(file_id, "\ndegenGeom(end).plate(%d).t = [", iplate);
	for ( int i = 0; i < nxsecs; i++ )
	{
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
			fprintf( file_id, "%f, ", degenPlate.t[i][j] );
		if ( i < nxsecs-1 )
			fprintf( file_id, "%f;\n", degenPlate.t[i][(num_pnts+1)/2-1] );
		else
			fprintf( file_id, "%f];\n", degenPlate.t[i][(num_pnts+1)/2-1] );
	}

	fprintf(file_id, "\ndegenGeom(end).plate(%d).nxCamber = [", iplate);
	for ( int i = 0; i < nxsecs; i++ )
	{
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
			fprintf( file_id, "%f, ", degenPlate.nCamber[i][j].x() );
		if ( i < nxsecs-1 )
			fprintf( file_id, "%f;\n", degenPlate.nCamber[i][(num_pnts+1)/2-1].x() );
		else
			fprintf( file_id, "%f];\n", degenPlate.nCamber[i][(num_pnts+1)/2-1].x() );
	}

	fprintf(file_id, "\ndegenGeom(end).plate(%d).nyCamber = [", iplate);
	for ( int i = 0; i < nxsecs; i++ )
	{
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
			fprintf( file_id, "%f, ", degenPlate.nCamber[i][j].y() );
		if ( i < nxsecs-1 )
			fprintf( file_id, "%f;\n", degenPlate.nCamber[i][(num_pnts+1)/2-1].y() );
		else
			fprintf( file_id, "%f];\n", degenPlate.nCamber[i][(num_pnts+1)/2-1].y() );
	}

	fprintf(file_id, "\ndegenGeom(end).plate(%d).nzCamber = [", iplate);
	for ( int i = 0; i < nxsecs; i++ )
	{
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
			fprintf( file_id, "%f, ", degenPlate.nCamber[i][j].z() );
		if ( i < nxsecs-1 )
			fprintf( file_id, "%f;\n", degenPlate.nCamber[i][(num_pnts+1)/2-1].z() );
		else
			fprintf( file_id, "%f];\n", degenPlate.nCamber[i][(num_pnts+1)/2-1].z() );
	}

	fprintf(file_id, "\ndegenGeom(end).plate(%d).u = [", iplate);
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenPlate.u[i]);
	}
	fprintf(file_id, "%f];", degenPlate.u[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).plate(%d).wTop = [", iplate);
	for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
	{
		fprintf(file_id, "%f, ", degenPlate.wTop[j]);
	}
	fprintf(file_id, "%f];", degenPlate.wTop[(num_pnts+1)/2-1]);

	fprintf(file_id, "\ndegenGeom(end).plate(%d).wBot = [", iplate);
	for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
	{
		fprintf(file_id, "%f, ", degenPlate.wBot[j]);
	}
	fprintf(file_id, "%f];", degenPlate.wBot[(num_pnts+1)/2-1]);

}

void DegenGeom::write_degenGeomStickM_file(FILE* file_id, int nxsecs, DegenStick &degenStick, int istick)
{
	//============================= DegenStick =============================//
	fprintf(file_id, "\ndegenGeom(end).stick(%d).Xle = [", istick);
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenStick.xle[i].x(),	\
				degenStick.xle[i].y(),	\
				degenStick.xle[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenStick.xle[nxsecs-1].x(),	\
			degenStick.xle[nxsecs-1].y(),	\
			degenStick.xle[nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick(%d).Xte = [", istick);
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenStick.xte[i].x(),	\
				degenStick.xte[i].y(),	\
				degenStick.xte[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenStick.xte[nxsecs-1].x(),	\
			degenStick.xte[nxsecs-1].y(),	\
			degenStick.xte[nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick(%d).XcgSolid = [", istick);
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenStick.xcgSolid[i].x(),	\
				degenStick.xcgSolid[i].y(),	\
				degenStick.xcgSolid[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenStick.xcgSolid[nxsecs-1].x(),	\
			degenStick.xcgSolid[nxsecs-1].y(),	\
			degenStick.xcgSolid[nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick(%d).XcgShell = [", istick);
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenStick.xcgShell[i].x(),	\
				degenStick.xcgShell[i].y(),	\
				degenStick.xcgShell[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenStick.xcgShell[nxsecs-1].x(),	\
			degenStick.xcgShell[nxsecs-1].y(),	\
			degenStick.xcgShell[nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick(%d).toc = [", istick);
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.toc[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.toc[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick(%d).tLoc = [", istick);
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.tLoc[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.tLoc[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick(%d).chord = [", istick);
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.chord[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.chord[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick(%d).sweep = [", istick);
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.sweep[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.sweep[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick(%d).Ishell = [", istick);
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(	file_id, "%f, %f, %f, %f, %f, %f;\n",		\
				degenStick.Ishell[i][0],	\
				degenStick.Ishell[i][1],	\
				degenStick.Ishell[i][2],	\
				degenStick.Ishell[i][3],	\
				degenStick.Ishell[i][4],	\
				degenStick.Ishell[i][5]		);
	}
	fprintf(	file_id, "%f, %f, %f, %f, %f, %f];\n",			\
			degenStick.Ishell[nxsecs-1][0],	\
			degenStick.Ishell[nxsecs-1][1],	\
			degenStick.Ishell[nxsecs-1][2],	\
			degenStick.Ishell[nxsecs-1][3],	\
			degenStick.Ishell[nxsecs-1][4],	\
			degenStick.Ishell[nxsecs-1][5]	);

	fprintf(file_id, "\ndegenGeom(end).stick(%d).Isolid = [", istick);
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(	file_id, "%f, %f, %f;\n",					\
				degenStick.Isolid[i][0],	\
				degenStick.Isolid[i][1],	\
				degenStick.Isolid[i][2]		);
	}
	fprintf(	file_id, "%f, %f, %f];\n",						\
			degenStick.Isolid[nxsecs-1][0],	\
			degenStick.Isolid[nxsecs-1][1],	\
			degenStick.Isolid[nxsecs-1][2]	);

	fprintf(file_id, "\ndegenGeom(end).stick(%d).area = [", istick);
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.area[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.area[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick(%d).areaNormal = [", istick);
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(	file_id, "%f, %f, %f;\n",						\
				degenStick.areaNormal[i].x(),	\
				degenStick.areaNormal[i].y(),	\
				degenStick.areaNormal[i].z()	);
	}
	fprintf(	file_id, "%f, %f, %f];\n",							\
			degenStick.areaNormal[nxsecs-1].x(),\
			degenStick.areaNormal[nxsecs-1].y(),\
			degenStick.areaNormal[nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick(%d).perimTop = [", istick);
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.perimTop[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.perimTop[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick(%d).perimBot = [", istick);
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.perimBot[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.perimBot[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick(%d).u = [", istick);
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.u[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.u[nxsecs-1]);
}

void DegenGeom::write_degenGeomPointM_file(FILE* file_id, int nxsecs)
{
	//============================= DegenPoint =============================//
	fprintf(file_id, "\ndegenGeom(end).point.vol = %f;\n", degenPoint.vol[0]);
	fprintf(file_id, "\ndegenGeom(end).point.volWet = %f;\n", degenPoint.volWet[0]);
	fprintf(file_id, "\ndegenGeom(end).point.area = %f;\n", degenPoint.area[0]);
	fprintf(file_id, "\ndegenGeom(end).point.areaWet = %f;\n", degenPoint.areaWet[0]);
	fprintf(file_id, "\ndegenGeom(end).point.Ishell = [%f, %f, %f, %f, %f, %f];\n",	\
			degenPoint.Ishell[0][0],	\
			degenPoint.Ishell[0][1],	\
			degenPoint.Ishell[0][2],	\
			degenPoint.Ishell[0][3],	\
			degenPoint.Ishell[0][4],	\
			degenPoint.Ishell[0][5]		);
	fprintf(file_id, "\ndegenGeom(end).point.Isolid = [%f, %f, %f, %f, %f, %f];\n",	\
			degenPoint.Isolid[0][0],	\
			degenPoint.Isolid[0][1],	\
			degenPoint.Isolid[0][2],	\
			degenPoint.Isolid[0][3],	\
			degenPoint.Isolid[0][4],	\
			degenPoint.Isolid[0][5]		);
	fprintf(file_id, "\ndegenGeom(end).point.xcgShell = [%f, %f, %f];\n",	\
			degenPoint.xcgShell[0].x(),	\
			degenPoint.xcgShell[0].y(),	\
			degenPoint.xcgShell[0].z()	);
	fprintf(file_id, "\ndegenGeom(end).point.xcgSolid = [%f, %f, %f];\n",	\
			degenPoint.xcgSolid[0].x(),	\
			degenPoint.xcgSolid[0].y(),	\
			degenPoint.xcgSolid[0].z()	);
}

void DegenGeom::write_degenGeomPropM_file(FILE* file_id)
{
	if ( parentGeom->getTypeStr() == "prop" )
	{
		fprintf(file_id, "\npropGeom(end).numBlades = %d;", degenProp.nblade);
		fprintf(file_id, "\npropGeom(end).diameter = %f;", degenProp.d);
		fprintf(file_id, "\npropGeom(end).rotCenter = [%f, %f, %f];",	degenProp.x.x(),	\
																		degenProp.x.y(),	\
																		degenProp.x.z()	);
		fprintf(file_id, "\npropGeom(end).rotVec = [%f, %f, %f];",		degenProp.nvec.x(),	\
																		degenProp.nvec.y(), \
																		degenProp.nvec.z()	);
	}
}

void DegenGeom::write_degenGeomM_file(FILE* file_id)
{
	int nxsecs = num_xsecs;

	if( type == SURFACE_TYPE )
	{
		fprintf(file_id, "\ndegenGeom(end+1).type = 'LIFTING_SURFACE';");
		fprintf(file_id, "\ndegenGeom(end).name = '%s';", name.c_str() );
	}
	else
	{
		fprintf(file_id, "\ndegenGeom(end+1).type = 'BODY';");
		fprintf(file_id, "\ndegenGeom(end).name = '%s';", name.c_str() );
	}

	write_degenGeomPropM_file(file_id);

	if ( parentGeom->getTypeStr() == "wing" || parentGeom->getTypeStr() == "prop" )
		nxsecs -= 2;

	write_degenGeomSurfM_file(file_id, nxsecs);

	write_degenGeomPlateM_file(file_id, nxsecs, degenPlates[0], 1);

	if ( type == DegenGeom::BODY_TYPE )
		write_degenGeomPlateM_file(file_id, nxsecs, degenPlates[1], 2);

	write_degenGeomStickM_file(file_id, nxsecs, degenSticks[0], 1);

	if ( type == DegenGeom::BODY_TYPE )
		write_degenGeomStickM_file(file_id, nxsecs, degenSticks[1], 2);

	write_degenGeomPointM_file(file_id, nxsecs);

}
