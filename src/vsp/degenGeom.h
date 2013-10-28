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

	void setNumXSecs( int nxss )		{ num_xsecs = nxss; }
	void setNumPnts( int npts )		{ num_pnts = npts; }

	Geom* getParentGeom()				{ return parentGeom; }
	void  setParentGeom( Geom* geom )	{ parentGeom = geom; }

	int getType()				{ return type; }
	void setType( int geomType)	{ type = geomType; }

	void write_degenGeomCsv_file(DegenGeom* degenGeom, FILE* file_id);
	void write_refl_degenGeomCsv_file(DegenGeom* degenGeom, FILE* file_id);

	void write_degenGeomM_file(DegenGeom* degenGeom, FILE* file_id);
	void write_refl_degenGeomM_file(DegenGeom* degenGeom, FILE* file_id);

protected:

	DegenSurface degenSurface;
	DegenPlate   degenPlate;
	DegenStick   degenStick;
	DegenPoint   degenPoint;

	int num_xsecs;
	int num_pnts;

	Geom* parentGeom;
	int   type;

};

#endif
