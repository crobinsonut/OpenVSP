
#include "xsec_surf.h"
#include "geom.h"
#include <cmath>


DegenGeom* Xsec_surf::createSurfDegenGeom(Geom* parentGeom, int sym_code_in, float mat[4][4], float refl_mat[4][4])
{
	DegenGeom*	degenGeom = new DegenGeom();
	degenGeom->setParentGeom(parentGeom);
	degenGeom->setType(DegenGeom::SURFACE_TYPE);

	degenGeom->setNumXSecs( num_xsecs );
	degenGeom->setNumPnts( num_pnts );
	degenGeom->setUarray( uArray );
	degenGeom->setWarray( wArray );
	degenGeom->setName( parentGeom->getName() );

	if ( sym_code_in != refl_pnts_xsecs_code )
	{
		refl_pnts_xsecs_code = sym_code_in;
		load_refl_pnts_xsecs();
	}

	degenGeom->createDegenSurface(sym_code_in, mat, pnts_xsecs, false);
	if(sym_code_in != NO_SYM)
		degenGeom->createDegenSurface(sym_code_in, refl_mat, refl_pnts_xsecs, true);

	degenGeom->createSurfDegenPlate(sym_code_in, mat, pnts_xsecs);
	if(sym_code_in != NO_SYM)
		degenGeom->createSurfDegenPlate(sym_code_in, refl_mat, refl_pnts_xsecs);

	degenGeom->createSurfDegenStick(sym_code_in, mat, pnts_xsecs);
	if ( sym_code_in != NO_SYM )
		degenGeom->createSurfDegenStick(sym_code_in, refl_mat, refl_pnts_xsecs);

	return degenGeom;
}

DegenGeom* Xsec_surf::createBodyDegenGeom(Geom* parentGeom, int sym_code_in, float mat[4][4], float refl_mat[4][4])
{
	DegenGeom*	degenGeom = new DegenGeom();
	degenGeom->setParentGeom(parentGeom);
	degenGeom->setType(DegenGeom::BODY_TYPE);

	degenGeom->setNumXSecs( num_xsecs );
	degenGeom->setNumPnts( num_pnts );
	degenGeom->setUarray( uArray );
	degenGeom->setWarray( wArray );
	degenGeom->setName( parentGeom->getName() );

	if ( sym_code_in != refl_pnts_xsecs_code )
	{
		refl_pnts_xsecs_code = sym_code_in;
		load_refl_pnts_xsecs();
	}

	degenGeom->createDegenSurface(sym_code_in, mat, pnts_xsecs, false);
	if(sym_code_in != NO_SYM)
		degenGeom->createDegenSurface(sym_code_in, refl_mat, refl_pnts_xsecs, true);

	degenGeom->createBodyDegenPlate(sym_code_in, mat, pnts_xsecs);
	if ( sym_code_in != NO_SYM )
		degenGeom->createBodyDegenPlate(sym_code_in, refl_mat, refl_pnts_xsecs);

	degenGeom->createBodyDegenStick(sym_code_in, mat, pnts_xsecs);
	if ( sym_code_in != NO_SYM )
		degenGeom->createBodyDegenStick(sym_code_in, refl_mat, refl_pnts_xsecs);

	return degenGeom;
}
