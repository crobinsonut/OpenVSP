//
//  geom_degen.h
//  VSP
//
//  Created by Joel Belben on 8/5/12.
//

#ifndef VSP_geom_degen_h
#define VSP_geom_degen_h

#include <vector>
#include "vec3d.h"
#include "vec2d.h"
#include "array_2d.h"


class Geom;

typedef struct {
    vector< vector< vec3d > >	x;			//!
    vector< vector< vec3d > >	nvec;		//!
    vector< double >        	u;			//!
    vector< double >        	w;			//!
} DegenSurface;

typedef struct {
    vector< vector< vec3d > >	x;			//!
    vector< vector< double > >	zcamber;	//!
    vector< vector< vec3d > >	nCamber;	//!
    vector< vector< double > >	t;			//!
    vector< vec3d >			 	nPlate;		//!
    vector< double >         	u;			//!
    vector< double >         	wTop;		//!
    vector< double >         	wBot;		//!
} DegenPlate;

typedef struct {
    vector< vec3d >			 	xle;		//!
    vector< vec3d >			 	xte;		//!
    vector< double >		 	toc;		//!
    vector< double >		 	tLoc;		//!
    vector< double >		 	chord;		//!
    vector< double >		 	sweep;		//! c/4. Make function to give sweep anywhere
    vector< vector< double > >	Ishell;		//!
    vector< vector< double > >	Isolid;		//!
    vector< vec3d >			 	xcgSolid;	//!
    vector< vec3d >			 	xcgShell;	//!
    vector< double >         	area;		//!
    vector< vec3d >			 	areaNormal;	//!
    vector< double >         	perimTop;	//!
    vector< double >         	perimBot;	//!
    vector< double >         	u;			//!
} DegenStick;

typedef struct {
    vector< double >		 	vol;		//!
    vector< double >		 	volWet;		//!
    vector< double >		 	area;		//!
    vector< double >		 	areaWet;	//!
    vector< vector< double > >	Ishell;		//! Multiply by rho*t to get inertias
    vector< vector< double > >	Isolid;		//! Multiply by rho to get inertias
    vector< vec3d >			 	xcgShell;	//!
    vector< vec3d >			 	xcgSolid;	//!
} DegenPoint;


class DegenGeom
{
public:
	enum{ XY_PLANE, XZ_PLANE, YZ_PLANE };

	enum { SURFACE_TYPE, BODY_TYPE };

	DegenGeom(){};
	virtual ~DegenGeom(){};

	DegenSurface getDegenSurface()  { return degenSurface; }
	DegenPlate   getDegenPlate()    { return degenPlate;   }
	DegenStick   getDegenStick()    { return degenStick;   }
	DegenPoint   getDegenPoint()    { return degenPoint;   }

	void setDegenSurface( DegenSurface degenSurface )	{ this->degenSurface = degenSurface; }
	void setDegenPlate(   DegenPlate   degenPlate )	{ this->degenPlate   = degenPlate;   }
	void setDegenStick(   DegenStick   degenStick )	{ this->degenStick   = degenStick;   }
	void setDegenPoint(   DegenPoint   degenPoint )	{ this->degenPoint   = degenPoint;   }

	int getNumXSecs()	{ return num_xsecs; };
	int getNumPnts()	{ return num_pnts; };

	void setNumXSecs( const int &nxss )		{ num_xsecs = nxss; }
	void setNumPnts( const int &npts )		{ num_pnts = npts; }

	void setUarray( const vector< double > &uarr )	{ uArray = uarr; }
	void setWarray( const vector< double > &warr )	{ wArray = warr; }

	Geom* getParentGeom()				{ return parentGeom; }
	void  setParentGeom( Geom* geom )	{ parentGeom = geom; }

	int getType()				{ return type; }
	void setType( int geomType)	{ type = geomType; }

	void write_degenGeomCsv_file(FILE* file_id);
	void write_refl_degenGeomCsv_file(FILE* file_id);

	void write_degenGeomM_file(FILE* file_id);
	void write_refl_degenGeomM_file(FILE* file_id);


	vec3d  get_area_normal( int ixs, const array_2d<vec3d> &pntsarr );
	double get_xsec_area( int ixs, const array_2d<vec3d> &pntsarr );
	double get_xsec_plane_area( int ixs, int plane, float mat[4][4], const array_2d<vec3d> &pntsarr );
	vec3d  get_xsec_centroid( int ixs, const array_2d<vec3d> &pntsarr );
	vec2d  get_xsec_centroid_in_plane(int ixs, int plane, float mat[4][4], const array_2d<vec3d> &pntsarr);

	void createDegenSurface(int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr, bool refl);
	void createSurfDegenPlate(int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr);
	void createBodyDegenPlate(int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr);
	void createSurfDegenStick(int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr);
	void createBodyDegenStick(int sym_code_in, float mat[4][4], const array_2d<vec3d> &pntsarr);

	vec3d get_xsec_shellCG( int ixs, const array_2d<vec3d> &pntsarr );

	vector<double> calculate_shell_inertias(int ixs, const array_2d<vec3d> &pntsarr);
	vector<double> calculate_shell_inertias_in_plane(int ixs, int plane, float mat[4][4], const array_2d<vec3d> &pntsarr);
	vector<double> calculate_solid_inertias(int ixs, const array_2d<vec3d> &pntsarr);
	vector<double> calculate_solid_inertias_in_plane(int ixs, int plane, float mat[4][4], const array_2d<vec3d> &pntsarr);

protected:

	DegenSurface degenSurface;
	DegenPlate   degenPlate;
	DegenStick   degenStick;
	DegenPoint   degenPoint;

	int num_xsecs;
	int num_pnts;

	vector< double > uArray;
	vector< double > wArray;

	Geom* parentGeom;
	int   type;

};



#endif
