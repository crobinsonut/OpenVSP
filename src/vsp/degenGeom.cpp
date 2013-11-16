#include "degenGeom.h"
#include "geom.h"
#include <cmath>
#include "writeMatlab.h"

void DegenGeom::build_trans_mat( vec3d x, vec3d y, vec3d z, const vec3d &p, double mat[4][4], double invmat[4][4] )
{
	// Initialize transformation matrix as identity.
	for( int i = 0; i < 4; i++ )
	{
		for( int j = 0; j < 4; j++ )
			mat[i][j] = 0;
		mat[i][i] = 1.0;
	}

	// Input vectors are assumed to be unit vectors.
	// Their magnitude will equal one unless it is the zero vector.
	if ( x.mag() < 0.5 || y.mag() < 0.5 || z.mag() < 0.5 )  // At least one input vector is degenerate.
	{
		if( x.mag() > 0.5 && y.mag() < 0.5 ) // X is OK, Y is degenerate
		{
			int imin = x.minor_comp();
			y.v[imin] = 1.0;
			z = cross( x, y );
			z.normalize();
			y = cross( z, x );
			y.normalize();
		}
		else if( x.mag() < 0.5 && y.mag() > 0.5 ) // X is degenerate, Y is OK
		{
			int imin = y.minor_comp();
			x.v[imin] = 1.0;
			z = cross( x, y );
			z.normalize();
			x = cross( y, z );
			x.normalize();
		}
		else
		{
			x = vec3d( 1, 0, 0 );
			y = vec3d( 0, 1, 0 );
			z = vec3d( 0, 0, 1 );
		}
	}

	// Build rotation part of matrix.
	for( int i = 0; i < 3; i++ )
	{
		mat[i][0] = x.v[i];
		mat[i][1] = y.v[i];
		mat[i][2] = z.v[i];

		invmat[0][i] = x.v[i];
		invmat[1][i] = y.v[i];
		invmat[2][i] = z.v[i];
	}

	// Build translation part of matrix.
	vec3d offset = p.transform(mat);

	for( int i = 0; i < 3; i++ )
	{
		mat[3][i] = -offset.v[i];
		invmat[3][i] = p.v[i];
	}
}

void DegenGeom::build_basis( const int &startPnt, const vector < vec3d > &sect, vec3d &v1, vec3d &v2, vec3d &v3 )
{
	int n_halfpts = (num_pnts + 1) / 2;
	int i_half = n_halfpts - 1;
	int i_quarter = ( (num_pnts + 1) / 4 ) - 1;

	// Vector along chord from te to le.
	v1 = sect[ startPnt % ( num_pnts - 1 ) ] - sect[ ( startPnt + i_half ) % ( num_pnts - 1 ) ];
	v1.normalize();

	// Vector from 1/4 to 3/4 points around section.
	vec3d v = ( sect[ ( startPnt + i_quarter ) % ( num_pnts - 1 ) ] - sect[ ( startPnt + i_half + i_quarter ) % ( num_pnts - 1) ]  );
	v.normalize();

	v3 = cross( v1, v );
	v3.normalize();

	v2 = cross( v3, v1 );
	v2.normalize();
}


void DegenGeom::transform_section( const int &startPnt, vector < vec3d > &sect, double trans[4][4], double invtrans[4][4] )
{
	vec3d chordVec;
	vec3d areaNormal;
	vec3d up;

	int n_halfpts = (num_pnts + 1) / 2;
	int i_half = n_halfpts - 1;

	build_basis( startPnt, sect, chordVec, up, areaNormal );

	vec3d p = sect[ ( startPnt + i_half ) % ( num_pnts - 1 ) ];

	build_trans_mat( chordVec, up, areaNormal, p, trans, invtrans );

	for( int i = 0; i < num_pnts; i++ )
		sect[i] = sect[i].transform(trans);
}

void DegenGeom::calculate_section_prop( const vector < vec3d > &sect, double &len, double &area, vec3d &xcgshell, vec3d &xcgsolid, vector < double > &Ishell, vector < double > &Isolid )
{
	double dl, da;
	double c11, c22, c12;
	double I22shell = 0, I11shell = 0, I12shell = 0;
	double I22solid = 0, I11solid = 0, I12solid = 0;

	vec3d c1;

	len = 0;
	area = 0;
	xcgshell = vec3d( 0, 0, 0 );
	xcgsolid = vec3d( 0, 0, 0 );
	Ishell.clear();
	Isolid.clear();

	int n = sect.size();

	// Cross section should be planar.  Assume Z=0.

	for( int j = 0; j < n; j++ )
	{
		int jnext = j + 1;
		if( jnext >= n ) // Final segment wrap to first segment.
			jnext = 0;

		dl = dist( sect[j], sect[jnext] );
		da = sect[j].x()*sect[jnext].y() - sect[jnext].x()*sect[j].y();

		len += dl;
		area += da;

		c1 = sect[j] + sect[jnext];
		xcgshell = xcgshell + c1 * dl;
		xcgsolid = xcgsolid + c1 * da;

		c11 = ( (sect[j].y() * sect[j].y()) + (sect[j].y() * sect[jnext].y()) + (sect[jnext].y() * sect[jnext].y()) );
		c22 = ( (sect[j].x() * sect[j].x()) + (sect[j].x() * sect[jnext].x()) + (sect[jnext].x() * sect[jnext].x()) );
		c12 = ( (2.0 * sect[j].x() * sect[j].y()) + (sect[j].x() * sect[jnext].y()) + (sect[jnext].x() * sect[j].y()) + (2.0 * sect[jnext].x() * sect[jnext].y()) );

		I11shell += c11 * dl;
		I22shell += c22 * dl;
		I12shell += c12 * dl;

		I11solid += c11 * da;
		I22solid += c22 * da;
		I12solid += c12 * da;
	}

	if( abs(len) < 1e-6 )
	{
		xcgshell = vec3d( 0, 0, 0 );

		I11shell = 0.0;
		I22shell = 0.0;
		I12shell = 0.0;
	}
	else
	{
		xcgshell = xcgshell / (2.0*len);

		I11shell /= (3.0*len);
		I22shell /= (3.0*len);
		I12shell /= (6.0*len);
	}

	if( abs(area) < 1e-6 )
	{
		area = 0.0;

		xcgsolid = vec3d( 0, 0, 0 );

		I11solid = 0.0;
		I22solid = 0.0;
		I12solid = 0.0;
	}
	else
	{
		area /= 2.0;

		xcgsolid = xcgsolid / (6.0*area);

		I11solid /= (12.0*area);
		I22solid /= (12.0*area);
		I12solid /= (24.0*area);
	}

	I11shell = ( I11shell - (xcgshell.y() * xcgshell.y() ) ) * len;
	I22shell = ( I22shell - (xcgshell.x() * xcgshell.x() ) ) * len;
	I12shell = ( I12shell - (xcgshell.x() * xcgshell.y() ) ) * len;

	I11solid = ( I11solid - (xcgsolid.y() * xcgsolid.y() ) ) * area;
	I22solid = ( I22solid - (xcgsolid.x() * xcgsolid.x() ) ) * area;
	I12solid = ( I12solid - (xcgsolid.x() * xcgsolid.y() ) ) * area;

	Ishell.push_back( I11shell );
	Ishell.push_back( I22shell );
	Ishell.push_back( I12shell );

	Isolid.push_back( I11solid );
	Isolid.push_back( I22solid );
	Isolid.push_back( I12solid );
}

void DegenGeom::createDegenSurface(int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr, const vector< double > &uArray, const vector< double > &wArray, bool refl)
{
	int nLow = 0, nHigh = num_xsecs;

	vec3d nVec;

	if ( getType() == DegenGeom::SURFACE_TYPE )
	{
		if ( getParentGeom()->getTypeStr() == "wing" || getParentGeom()->getTypeStr() == "prop" )
		{
			nLow  = 1;
			nHigh = num_xsecs - 1;
		}
	}

	int nxs = nHigh - nLow;
	degenSurface.x.resize( nxs );
	degenSurface.u.resize( nxs );

	for ( int i = nLow; i < nHigh; i++ )
	{
		degenSurface.x[i-nLow].resize( num_pnts );
		for ( int j = 0; j < num_pnts; j++ )
		{
			degenSurface.x[i-nLow][j] = pntsarr(i,j).transform(mat);
		}
		degenSurface.u[i-nLow] = uArray[i];
	}

	degenSurface.nvec.resize( nxs - 1 );
	degenSurface.area.resize( nxs - 1 );
	for ( int i = 0; i < nxs-1; i++ )
	{
		degenSurface.nvec[i].resize( num_pnts - 1 );
		degenSurface.area[i].resize( num_pnts - 1 );
		for ( int j = 0; j < num_pnts-1; j++ )
		{
			vec3d sVec1 = degenSurface.x[i+1][j] - degenSurface.x[i][j];
			vec3d sVec2 = degenSurface.x[i][j+1] - degenSurface.x[i][j];

			vec3d sVec3 = degenSurface.x[i+1][j+1] - degenSurface.x[i+1][j];
			vec3d sVec4 = degenSurface.x[i+1][j+1] - degenSurface.x[i][j+1];

			// Take areas of both triangles for quad.
			degenSurface.area[i][j] = 0.5*( cross(sVec1, sVec2).mag() + cross(sVec3, sVec4).mag() );

			vec3d v1, v2;

			if( sVec1.mag() > sVec4.mag() )
				v1 = sVec1;
			else
				v1 = sVec4;

			if( sVec2.mag() > sVec3.mag() )
				v2 = sVec2;
			else
				v2 = sVec3;

			// Approximate normal as normal for two longest sides.
			nVec = cross( v1, v2 );
			if(refl)
				nVec = nVec*-1.0;

			nVec.normalize();

			degenSurface.nvec[i][j] = nVec;
		}
	}

	degenSurface.w.resize( num_pnts );
	for ( int j = 0; j < num_pnts; j++ )
	{
		degenSurface.w[j] = wArray[j];
	}
}

void DegenGeom::createSurfDegenPlate(int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr, const vector< double > &uArray, const vector< double > &wArray)
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
	createDegenPlate(degenPlates[0], sym_code_in, mat, pntsarr, uArray, wArray, nLow, nHigh, startPnt);
}

void DegenGeom::createBodyDegenPlate(int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr, const vector< double > &uArray, const vector< double > &wArray)
{
	int nLow = 0, nHigh = num_xsecs;

	degenPlates.resize(2);

	int startPnt = 0;
	createDegenPlate(degenPlates[0], sym_code_in, mat, pntsarr, uArray, wArray, nLow, nHigh, startPnt);

	startPnt = (num_pnts - 1) / 4;
	createDegenPlate(degenPlates[1], sym_code_in, mat, pntsarr, uArray, wArray, nLow, nHigh, startPnt);
}

void DegenGeom::createDegenPlate(DegenPlate &degenPlate, int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr, const vector< double > &uArray, const vector< double > &wArray, int nLow, int nHigh, int startPnt)
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


		vector < vec3d > section( num_pnts );
		for ( int j = 0; j < num_pnts; j++ )
			section[j] = pntsarr( i, j ).transform(mat);

		vec3d chordVec, anv, up;
		build_basis( startPnt, section, chordVec, up, anv );

		//== Compute plate normal ==//
		// rotated chord vector
		vec3d rcv = xVec[platePnts-1] - xVec[0];
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

void DegenGeom::createSurfDegenStick(int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr, const vector< double > &uArray, const vector< double > &wArray)
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
	createDegenStick(degenSticks[0], sym_code_in, mat, pntsarr, uArray, wArray, nLow, nHigh, startPnt);

}

void DegenGeom::createBodyDegenStick(int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr, const vector< double > &uArray, const vector< double > &wArray)
{
	int nLow = 0, nHigh = num_xsecs;

	degenSticks.resize(2);

	int startPnt = 0;
	createDegenStick(degenSticks[0], sym_code_in, mat, pntsarr, uArray, wArray, nLow, nHigh, startPnt);

	startPnt = (num_pnts - 1) / 4;
	createDegenStick(degenSticks[1], sym_code_in, mat, pntsarr, uArray, wArray, nLow, nHigh, startPnt);
}


void DegenGeom::createDegenStick(DegenStick &degenStick, int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr, const vector< double > &uArray, const vector< double > &wArray, int nLow, int nHigh, int startPnt)
{
	int platePnts = (num_pnts + 1) / 2;
	vec3d chordVec, camberPnt;

	for ( int i = nLow; i < nHigh; i++ )
	{
		vec3d topPnt, botPnt;
		double tempThickness = -1, perimTop = 0, perimBot = 0;
		int    maxThickIdx[2] = {0,0};

		vector < vec3d > section( num_pnts );
		for ( int j = 0; j < num_pnts; j++ )
			section[j] = pntsarr( i, j ).transform(mat);

		vec3d chordVec, areaNormal, up;
		build_basis( startPnt, section, chordVec, up, areaNormal );
		degenStick.sectnvec.push_back( areaNormal );

		// Transform section to XY plane, with LE at origin and TE on +X-axis.
		// Also return transformation matrix from 3D space to XY plane as sect_trans.
		double sect_trans[4][4], sect_inv_trans[4][4];
		transform_section( startPnt, section, sect_trans, sect_inv_trans );

		// Reshape matrix to vector and store in degenStick
		vector < double > tmatvec;
		for ( int j = 0; j < 4; j++ )
			for ( int k = 0; k < 4; k++ )
				tmatvec.push_back( sect_trans[j][k] );

		degenStick.transmat.push_back( tmatvec );

		// normalized, unrotated chord vector (te->le)
		chordVec = pntsarr(i, startPnt+platePnts-1) - pntsarr(i, startPnt);
		chordVec.normalize();

		degenStick.xle.push_back( pntsarr( i, startPnt+platePnts-1 ).transform(mat) );
		degenStick.xte.push_back( pntsarr( i, startPnt ).transform(mat) );

		double chord = dist(pntsarr(i, startPnt+platePnts-1), pntsarr(i, startPnt));
		degenStick.chord.push_back( chord );

		degenStick.u.push_back( uArray[i] );

		double len, area;
		vec3d xcgshell, xcgsolid;
		vector < double > Ishell, Isolid;
		calculate_section_prop( section, len, area, xcgshell, xcgsolid, Ishell, Isolid );

		xcgshell = xcgshell.transform(sect_inv_trans);
		xcgsolid = xcgsolid.transform(sect_inv_trans);

		degenStick.Ishell.push_back( Ishell );
		degenStick.Isolid.push_back( Isolid );
		degenStick.xcgShell.push_back( xcgshell );
		degenStick.xcgSolid.push_back( xcgsolid );
		degenStick.sectarea.push_back( area );

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
		}

		camberPnt = ( pntsarr(i,maxThickIdx[0]) + pntsarr(i,maxThickIdx[1]) ) / 2;

		if( chord > 0 )
		{
			degenStick.tLoc.push_back( 1 - (dot(camberPnt-pntsarr(i,startPnt),chordVec) / chord) );
			degenStick.toc.push_back( tempThickness / chord );
		}
		else
		{
			degenStick.tLoc.push_back( 0.0 );
			degenStick.toc.push_back( 0.0 );
		}

		perimTop += dist( pntsarr(i, startPnt+platePnts-1), pntsarr(i, startPnt+platePnts-2) );
		perimBot += dist( pntsarr(i, startPnt+platePnts), pntsarr(i, startPnt+platePnts-1) );
		degenStick.perimTop.push_back( perimTop );
		degenStick.perimBot.push_back( perimBot );
	}

	// Calculate sweep angle
	for( int i = nLow; i < nHigh-1; i++ )
	{
		vec3d xle0 = pntsarr( i, startPnt+platePnts-1 ).transform( mat );
		vec3d xte0 = pntsarr( i, startPnt ).transform( mat );

		vec3d xle1 = pntsarr( i+1, startPnt+platePnts-1 ).transform( mat );
		vec3d xte1 = pntsarr( i+1, startPnt ).transform( mat );

		vec3d vle = xle1 - xle0;
		vle.normalize();

		vec3d vte = xte1 - xte0;
		vte.normalize();

		vec3d vchd = xte0 - xle0;
		vchd.normalize();

		vec3d vdownstream( 1, 0, 0 );

		if( vchd.mag() > 0.5 )
		{
			vec3d n = cross( vchd, vle );
			n.normalize();

			vec3d downNormal = cross( n, vdownstream  );
			downNormal.normalize();

			if( downNormal.mag() > 0.5 && vle.mag() > 0.5 )
				degenStick.sweeple.push_back( RAD_2_DEG*signed_angle( downNormal, vle, n * -1.0 ) );
			else
				degenStick.sweeple.push_back( 0.0 );

			if( downNormal.mag() > 0.5 && vte.mag() > 0.5 )
				degenStick.sweepte.push_back( RAD_2_DEG*signed_angle( downNormal, vte, n * -1.0 ) );
			else
				degenStick.sweepte.push_back( 0.0 );
		}
		else
		{
			vec3d vref1 = cross( vdownstream, vle );
			vec3d vref2 = cross( vref1, vdownstream );
			degenStick.sweeple.push_back( RAD_2_DEG*signed_angle( vref2, vle, vref1 * -1.0 ) );

			vref1 = cross( vdownstream, vte );
			vref2 = cross( vref1, vdownstream );
			degenStick.sweepte.push_back( RAD_2_DEG*signed_angle( vref2, vte, vref1 * -1.0 ) );
		}



		double areaTop = 0, areaBot = 0;
		for ( int j = 0, k = startPnt+num_pnts-2; j < platePnts-1; j++, k-- )
		{
			areaTop += degenSurface.area[i-nLow][(startPnt+j) % (num_pnts-1)];
			areaBot += degenSurface.area[i-nLow][k % (num_pnts-1)];
		}
		degenStick.areaTop.push_back( areaTop );
		degenStick.areaBot.push_back( areaBot );
	}
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
	fprintf(file_id, "# nx,ny,nz,area\n");

	for ( int i = 0; i < nxsecs-1; i++ )
	{
		for ( int j = 0; j < num_pnts-1; j++ )
		{
			fprintf(file_id, makeCsvFmt(4),			\
					degenSurface.nvec[i][j].x(),	\
					degenSurface.nvec[i][j].y(),	\
					degenSurface.nvec[i][j].z(),	\
					degenSurface.area[i][j]			);
		}
	}


}

void DegenGeom::write_degenGeomPlateCsv_file(FILE* file_id, int nxsecs, DegenPlate &degenPlate)
{
	fprintf(file_id, "# DegenGeom Type,nXsecs,nPnts/Xsec\n");
	fprintf(file_id, "PLATE,%d,%d\n", nxsecs, (num_pnts+1)/2);
	fprintf(file_id,"# nx,ny,nz\n");
	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf(file_id, makeCsvFmt(3), degenPlate.nPlate[i].x(), \
				degenPlate.nPlate[i].y(), \
				degenPlate.nPlate[i].z()  );
	}

	fprintf(file_id, "# x,y,z,zCamber,t,nCamberx,nCambery,nCamberz,u,wTop,wBot\n");
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

	fprintf(file_id, "# DegenGeom Type, nXsecs\nSTICK_NODE, %d\n# lex,ley,lez,tex,tey,tez,cgShellx,cgShelly,cgShellz,"
			"cgSolidx,cgSolidy,cgSolidz,toc,tLoc,chord,Ishell11,Ishell22,"
			"Ishell12,Isolid11,Isolid22,Isolid12,sectArea,sectNormalx,"
			"sectNormaly,sectNormalz,perimTop,perimBot,u\n", nxsecs);

	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf(file_id, makeCsvFmt(28),	\
				degenStick.xle[i].x(),					\
				degenStick.xle[i].y(),					\
				degenStick.xle[i].z(),					\
				degenStick.xte[i].x(),					\
				degenStick.xte[i].y(),					\
				degenStick.xte[i].z(),					\
				degenStick.xcgShell[i].x(),				\
				degenStick.xcgShell[i].y(),				\
				degenStick.xcgShell[i].z(),				\
				degenStick.xcgSolid[i].x(),				\
				degenStick.xcgSolid[i].y(),				\
				degenStick.xcgSolid[i].z(),				\
				degenStick.toc[i],					\
				degenStick.tLoc[i],					\
				degenStick.chord[i],					\
				degenStick.Ishell[i][0],				\
				degenStick.Ishell[i][1],				\
				degenStick.Ishell[i][2],				\
				degenStick.Isolid[i][0],				\
				degenStick.Isolid[i][1],				\
				degenStick.Isolid[i][2],				\
				degenStick.sectarea[i],					\
				degenStick.sectnvec[i].x(),				\
				degenStick.sectnvec[i].y(),				\
				degenStick.sectnvec[i].z(),				\
				degenStick.perimTop[i],					\
				degenStick.perimBot[i],					\
				degenStick.u[i]						);
	}


	fprintf(file_id, "# DegenGeom Type, nXsecs\nSTICK_FACE, %d\n# sweeple,sweepte,areaTop,areaBot\n", nxsecs-1);

	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, makeCsvFmt(4),	\
				degenStick.sweeple[i],					\
				degenStick.sweepte[i],					\
				degenStick.areaTop[i],					\
				degenStick.areaBot[i]					);
	}
}

void DegenGeom::write_degenGeomPointCsv_file(FILE* file_id, int nxsecs)
{
	fprintf(file_id, "# DegenGeom Type\n");
	fprintf(file_id, "POINT\n");
	fprintf(file_id, "# vol,volWet,area,areaWet,Ishellxx,Ishellyy,Ishellzz,Ishellxy,");
	fprintf(file_id, "Ishellxz,Ishellyz,Isolidxx,Isolidyy,Isolidzz,Isolidxy,Isolidxz,");
	fprintf(file_id, "Isolidyz,cgShellx,cgShelly,cgShellz,cgSolidx,cgSolidy,cgSolidz\n");
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
		fprintf(file_id, "# numBlades,diameter,x,y,z,nx,ny,nz\n");
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
	string basename = string("degenGeom(end).surf.");

	WriteVecDoubleM writeVecDouble;
	WriteMatVec3dM writeMatVec3d;
	WriteMatDoubleM writeMatDouble;

	writeMatVec3d.write(  file_id, degenSurface.x,    basename, nxsecs, num_pnts );
	writeVecDouble.write( file_id, degenSurface.u,    basename + "u",   nxsecs );
	writeVecDouble.write( file_id, degenSurface.w,    basename + "w",   num_pnts );
	writeMatVec3d.write(  file_id, degenSurface.nvec, basename + "n",   nxsecs-1,    num_pnts-1 );
	writeMatDouble.write( file_id, degenSurface.area, basename + "area",nxsecs-1,    num_pnts-1 );
}

void DegenGeom::write_degenGeomPlateM_file(FILE* file_id, int nxsecs, DegenPlate &degenPlate, int iplate)
{
	char num[80];
	sprintf( num, "degenGeom(end).plate(%d).", iplate );
	string basename = string(num);

	WriteVecDoubleM writeVecDouble;
	WriteVecVec3dM writeVecVec3d;
	WriteMatDoubleM writeMatDouble;
	WriteMatVec3dM writeMatVec3d;

	writeVecVec3d.write(  file_id, degenPlate.nPlate,  basename + "n",       nxsecs );
	writeMatVec3d.write(  file_id, degenPlate.x,       basename,             nxsecs,           (num_pnts+1)/2 );
	writeMatDouble.write( file_id, degenPlate.zcamber, basename + "zCamber", nxsecs,           (num_pnts+1)/2 );
	writeMatDouble.write( file_id, degenPlate.t,       basename + "t",       nxsecs,           (num_pnts+1)/2 );
	writeMatVec3d.write(  file_id, degenPlate.nCamber, basename + "nCamber", nxsecs,           (num_pnts+1)/2 );
	writeVecDouble.write( file_id, degenPlate.u,       basename + "u",       nxsecs );
	writeVecDouble.write( file_id, degenPlate.wTop,    basename + "wTop",    (num_pnts+1)/2 );
	writeVecDouble.write( file_id, degenPlate.wBot,    basename + "wBot",    (num_pnts+1)/2 );
}

void DegenGeom::write_degenGeomStickM_file(FILE* file_id, int nxsecs, DegenStick &degenStick, int istick)
{
	char num[80];
	sprintf( num, "degenGeom(end).stick(%d).", istick );
	string basename = string(num);

	WriteVecDoubleM writeVecDouble;
	WriteVecVec3dM writeVecVec3d;
	WriteMatDoubleM writeMatDouble;

	writeVecVec3d.write(  file_id, degenStick.xle,        basename + "le",         nxsecs );
	writeVecVec3d.write(  file_id, degenStick.xte,        basename + "te",         nxsecs );
	writeVecVec3d.write(  file_id, degenStick.xcgShell,   basename + "cgShell",    nxsecs );
	writeVecVec3d.write(  file_id, degenStick.xcgSolid,   basename + "cgSolid",    nxsecs );
	writeVecDouble.write( file_id, degenStick.toc,        basename + "toc",        nxsecs );
	writeVecDouble.write( file_id, degenStick.tLoc,       basename + "tLoc",       nxsecs );
	writeVecDouble.write( file_id, degenStick.chord,      basename + "chord",      nxsecs );
	writeMatDouble.write( file_id, degenStick.Ishell,     basename + "Ishell",     nxsecs,        3 );
	writeMatDouble.write( file_id, degenStick.Isolid,     basename + "Isolid",     nxsecs,        3 );
	writeVecDouble.write( file_id, degenStick.sectarea,   basename + "sectArea",   nxsecs );
	writeVecVec3d.write(  file_id, degenStick.sectnvec,   basename + "sectNormal", nxsecs );
	writeVecDouble.write( file_id, degenStick.perimTop,   basename + "perimTop",   nxsecs );
	writeVecDouble.write( file_id, degenStick.perimBot,   basename + "perimBot",   nxsecs );
	writeVecDouble.write( file_id, degenStick.u,          basename + "u",          nxsecs );

	writeVecDouble.write( file_id, degenStick.sweeple,    basename + "sweeple",    nxsecs - 1 );
	writeVecDouble.write( file_id, degenStick.sweepte,    basename + "sweepte",    nxsecs - 1 );
	writeVecDouble.write( file_id, degenStick.areaTop,    basename + "areaTop",    nxsecs - 1 );
	writeVecDouble.write( file_id, degenStick.areaBot,    basename + "areaBot",    nxsecs - 1 );

}

void DegenGeom::write_degenGeomPointM_file(FILE* file_id, int nxsecs)
{
	string basename = string("degenGeom(end).point.");

	WriteDoubleM writeDouble;
	WriteVec3dM writeVec3d;
	WriteVecDoubleM writeVecDouble;

	writeDouble.write(    file_id, degenPoint.vol[0],      basename + "vol" );
	writeDouble.write(    file_id, degenPoint.volWet[0],   basename + "volWet" );
	writeDouble.write(    file_id, degenPoint.area[0],     basename + "area" );
	writeDouble.write(    file_id, degenPoint.areaWet[0],  basename + "areaWet" );
	writeVecDouble.write( file_id, degenPoint.Ishell[0],   basename + "Ishell",     6 );
	writeVecDouble.write( file_id, degenPoint.Isolid[0],   basename + "Isolid",     6 );
	writeVec3d.write(     file_id, degenPoint.xcgShell[0], basename + "cgShell");
	writeVec3d.write(     file_id, degenPoint.xcgSolid[0], basename + "cgSolid");
}

void DegenGeom::write_degenGeomPropM_file(FILE* file_id)
{
	if ( parentGeom->getTypeStr() == "prop" )
	{
		string basename = string("propGeom(end).");

		fprintf(file_id, "propGeom(end).numBlades = %d;\n", degenProp.nblade);

		WriteDoubleM writeDouble;
		WriteVec3dM writeVec3d;

		writeDouble.write( file_id, degenProp.d,    basename + "diameter" );
		writeVec3d.write(  file_id, degenProp.x,    basename);
		writeVec3d.write(  file_id, degenProp.nvec, basename + "n");
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
