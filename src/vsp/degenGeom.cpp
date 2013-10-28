#include "degenGeom.h"
#include "geom.h"
#include <cmath>

void DegenGeom::write_degenGeomCsv_file(DegenGeom* degenGeom, FILE* file_id)
{
	int nxsecs = num_xsecs;

	if ( degenGeom->getParentGeom()->getTypeStr() == "wing" || degenGeom->getParentGeom()->getTypeStr() == "prop" )
		nxsecs -= 2;
	int nxsecsOrig = nxsecs;

	fprintf(file_id, "# DegenGeom Type,nXsecs, nPnts/Xsec\n");
	fprintf(file_id, "FULL_SURFACE,%d,%d\n", nxsecs, num_pnts);
	fprintf(file_id, "# x,y,z,xn,yn,zn,u,w\n");

	for ( int i = 0; i < nxsecs; i++ )
	{
		for ( int j = 0; j < num_pnts; j++ )
		{
			fprintf(	file_id, "%f,%f,%f,%f,%f,%f,%f,%f\n",			\
					degenGeom->getDegenSurface().x[i][j].x(),		\
					degenGeom->getDegenSurface().x[i][j].y(),		\
					degenGeom->getDegenSurface().x[i][j].z(),		\
					degenGeom->getDegenSurface().nvec[i][j].x(),	\
					degenGeom->getDegenSurface().nvec[i][j].y(),	\
					degenGeom->getDegenSurface().nvec[i][j].z(),	\
					degenGeom->getDegenSurface().u[i],				\
					degenGeom->getDegenSurface().w[j]				);
		}
	}

	// JBB: Twice as many cross sections for BODY type
	if ( degenGeom->getType() == DegenGeom::BODY_TYPE ) nxsecs *= 2;

	fprintf(file_id, "# DegenGeom Type,nXsecs,nPnts/Xsec\n");
	fprintf(file_id, "PLATE,%d,%d\n", nxsecs, (num_pnts+1)/2);
	fprintf(file_id,"# xn,yn,zn\n");
	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf( file_id, "%f,%f,%f\n", degenGeom->getDegenPlate().nPlate[i].x(), \
				degenGeom->getDegenPlate().nPlate[i].y(), \
				degenGeom->getDegenPlate().nPlate[i].z()  );
	}

	fprintf(file_id, "# x,y,z,zCamb,t,nCamrX,nCambY,nCambZ,u,wTop,wBot\n");
	for ( int i = 0; i < nxsecsOrig; i++ )
	{
		for ( int j = 0; j < (num_pnts+1)/2; j++ )
		{
			fprintf(	file_id, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",	\
					degenGeom->getDegenPlate().x[i][j].x(),				\
					degenGeom->getDegenPlate().x[i][j].y(),				\
					degenGeom->getDegenPlate().x[i][j].z(),				\
					degenGeom->getDegenPlate().zcamber[i][j],			\
					degenGeom->getDegenPlate().t[i][j],					\
					degenGeom->getDegenPlate().nCamber[i][j].x(),		\
					degenGeom->getDegenPlate().nCamber[i][j].y(),		\
					degenGeom->getDegenPlate().nCamber[i][j].z(),		\
					degenGeom->getDegenPlate().u[i],					\
					degenGeom->getDegenPlate().wTop[j],					\
					degenGeom->getDegenPlate().wBot[j]					);
		}
	}
	for ( int i = nxsecsOrig; i < nxsecs; i++ )
	{
		for ( int j = 0; j < (num_pnts+1)/2; j++ )
		{
			fprintf(	file_id, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",	\
					degenGeom->getDegenPlate().x[i][j].x(),				\
					degenGeom->getDegenPlate().x[i][j].y(),				\
					degenGeom->getDegenPlate().x[i][j].z(),				\
					degenGeom->getDegenPlate().zcamber[i][j],			\
					degenGeom->getDegenPlate().t[i][j],					\
					degenGeom->getDegenPlate().nCamber[i][j].x(),		\
					degenGeom->getDegenPlate().nCamber[i][j].y(),		\
					degenGeom->getDegenPlate().nCamber[i][j].z(),		\
					degenGeom->getDegenPlate().u[i],					\
					degenGeom->getDegenPlate().wTop[(num_pnts+1)/2+j],	\
					degenGeom->getDegenPlate().wBot[(num_pnts+1)/2+j]	);
		}
	}

	fprintf(file_id, "# DegenGeom Type, nXsecs\nSTICK, %d\n# xle,yle,zle,xte,yte,zte,xcg_solid,ycg_solid,zcg_solid,"
			"xcg_shell,ycg_shell,zcg_shell,toc,tLoc,chord,sweep,Ixx_shell_A,Ixx_shell_B,Izz_shell_A,"
			"Izz_shell_B,J_shell_A,J_shell_B,Ixx_solid,Izz_solid,J_solid,area,areaNormalX,"
			"areaNormalY,areaNormalZ,perimTop,perimBot,u\n", nxsecs);

	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf(	file_id, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",	\
				degenGeom->getDegenStick().xle[i].x(),					\
				degenGeom->getDegenStick().xle[i].y(),					\
				degenGeom->getDegenStick().xle[i].z(),					\
				degenGeom->getDegenStick().xte[i].x(),					\
				degenGeom->getDegenStick().xte[i].y(),					\
				degenGeom->getDegenStick().xte[i].z(),					\
				degenGeom->getDegenStick().xcgSolid[i].x(),				\
				degenGeom->getDegenStick().xcgSolid[i].y(),				\
				degenGeom->getDegenStick().xcgSolid[i].z(),				\
				degenGeom->getDegenStick().xcgShell[i].x(),				\
				degenGeom->getDegenStick().xcgShell[i].y(),				\
				degenGeom->getDegenStick().xcgShell[i].z(),				\
				degenGeom->getDegenStick().toc[i],						\
				degenGeom->getDegenStick().tLoc[i],						\
				degenGeom->getDegenStick().chord[i],					\
				degenGeom->getDegenStick().sweep[i],					\
				degenGeom->getDegenStick().Ishell[i][0],				\
				degenGeom->getDegenStick().Ishell[i][1],				\
				degenGeom->getDegenStick().Ishell[i][2],				\
				degenGeom->getDegenStick().Ishell[i][3],				\
				degenGeom->getDegenStick().Ishell[i][4],				\
				degenGeom->getDegenStick().Ishell[i][5],				\
				degenGeom->getDegenStick().Isolid[i][0],				\
				degenGeom->getDegenStick().Isolid[i][1],				\
				degenGeom->getDegenStick().Isolid[i][2],				\
				degenGeom->getDegenStick().area[i],						\
				degenGeom->getDegenStick().areaNormal[i].x(),			\
				degenGeom->getDegenStick().areaNormal[i].y(),			\
				degenGeom->getDegenStick().areaNormal[i].z(),			\
				degenGeom->getDegenStick().perimTop[i],					\
				degenGeom->getDegenStick().perimBot[i],					\
				degenGeom->getDegenStick().u[i]							);
	}

	nxsecs = nxsecsOrig;
	fprintf(file_id, "# DegenGeom Type\n");
	fprintf(file_id, "POINT\n");
	fprintf(file_id, "# vol,volWet,area,areaWet,IxxShell,IyyShell,IzzShell,IxyShell,");
	fprintf(file_id, "IxzShell,IyzShell,IxxSolid,IyySolid,IzzSolid,IxySolid,IxzSolid,");
	fprintf(file_id, "IyzSolid,xcgShell,ycgShell,zcgShell,xcgSolid,ycgSolid,zcgSolid\n");
	fprintf(	file_id, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",\
			degenGeom->getDegenPoint().vol[0],			\
			degenGeom->getDegenPoint().volWet[0],		\
			degenGeom->getDegenPoint().area[0],			\
			degenGeom->getDegenPoint().areaWet[0],		\
			degenGeom->getDegenPoint().Ishell[0][0],	\
			degenGeom->getDegenPoint().Ishell[0][1],	\
			degenGeom->getDegenPoint().Ishell[0][2],	\
			degenGeom->getDegenPoint().Ishell[0][3],	\
			degenGeom->getDegenPoint().Ishell[0][4],	\
			degenGeom->getDegenPoint().Ishell[0][5],	\
			degenGeom->getDegenPoint().Isolid[0][0],	\
			degenGeom->getDegenPoint().Isolid[0][1],	\
			degenGeom->getDegenPoint().Isolid[0][2],	\
			degenGeom->getDegenPoint().Isolid[0][3],	\
			degenGeom->getDegenPoint().Isolid[0][4],	\
			degenGeom->getDegenPoint().Isolid[0][5],	\
			degenGeom->getDegenPoint().xcgShell[0].x(),	\
			degenGeom->getDegenPoint().xcgShell[0].y(),	\
			degenGeom->getDegenPoint().xcgShell[0].z(),	\
			degenGeom->getDegenPoint().xcgSolid[0].x(),	\
			degenGeom->getDegenPoint().xcgSolid[0].y(),	\
			degenGeom->getDegenPoint().xcgSolid[0].z()	);
}

void DegenGeom::write_refl_degenGeomCsv_file(DegenGeom* degenGeom, FILE* file_id)
{
	int nxsecs = num_xsecs;

	if ( degenGeom->getParentGeom()->getTypeStr() == "wing" || degenGeom->getParentGeom()->getTypeStr() == "prop" )
		nxsecs -= 2;
	int nxsecsOrig = nxsecs;

	fprintf(file_id, "# DegenGeom Type,nXsecs,nPnts/Xsec\n");
	fprintf(file_id, "FULL_SURFACE,%d,%d\n", nxsecs, num_pnts);
	fprintf(file_id, "# x,y,z,xn,yn,zn,u,w\n");

	for ( int i = nxsecs; i < 2 * nxsecs; i++ )
	{
		for ( int j = 0; j < num_pnts; j++ )
		{
			fprintf(	file_id, "%f,%f,%f,%f,%f,%f,%f,%f\n",		\
					degenGeom->getDegenSurface().x[i][j].x(),		\
					degenGeom->getDegenSurface().x[i][j].y(),		\
					degenGeom->getDegenSurface().x[i][j].z(),		\
					degenGeom->getDegenSurface().nvec[i][j].x(),	\
					degenGeom->getDegenSurface().nvec[i][j].y(),	\
					degenGeom->getDegenSurface().nvec[i][j].z(),	\
					degenGeom->getDegenSurface().u[i],				\
					degenGeom->getDegenSurface().w[j]				);
		}
	}

	// Twice as many cross sections for BODY type
	if ( degenGeom->getType() == DegenGeom::BODY_TYPE ) nxsecs *= 2;
	fprintf(file_id, "# DegenGeom Type,nXsecs,nPnts/Xsec\n");
	fprintf(file_id, "PLATE,%d,%d\n", nxsecs, (num_pnts+1)/2);
	fprintf(file_id,"# xn,yn,zn\n");
	for ( int i = nxsecs; i < 2 * nxsecs; i++ )
	{
		fprintf( file_id, "%f,%f,%f\n", degenGeom->getDegenPlate().nPlate[i].x(), \
				degenGeom->getDegenPlate().nPlate[i].y(), \
				degenGeom->getDegenPlate().nPlate[i].z()  );
	}

	fprintf(file_id, "# x,y,z,zCamb,t,nCambX,nCambY,nCambZ,u,wTop,wBot\n");
	for ( int i = nxsecs; i < nxsecsOrig + nxsecs; i++ )
	{
		for ( int j = 0; j < (num_pnts+1)/2; j++ )
		{
			fprintf(	file_id, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",	\
					degenGeom->getDegenPlate().x[i][j].x(),				\
					degenGeom->getDegenPlate().x[i][j].y(),				\
					degenGeom->getDegenPlate().x[i][j].z(),				\
					degenGeom->getDegenPlate().zcamber[i][j],			\
					degenGeom->getDegenPlate().t[i][j],					\
					degenGeom->getDegenPlate().nCamber[i][j].x(),		\
					degenGeom->getDegenPlate().nCamber[i][j].y(),		\
					degenGeom->getDegenPlate().nCamber[i][j].z(),		\
					degenGeom->getDegenPlate().u[i],					\
					degenGeom->getDegenPlate().wTop[j],					\
					degenGeom->getDegenPlate().wBot[j]					);
		}
	}
	for ( int i = nxsecsOrig + nxsecs; i < 2*nxsecs; i++ )
	{
		for ( int j = 0; j < (num_pnts+1)/2; j++ )
		{
			fprintf(	file_id, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",	\
					degenGeom->getDegenPlate().x[i][j].x(),				\
					degenGeom->getDegenPlate().x[i][j].y(),				\
					degenGeom->getDegenPlate().x[i][j].z(),				\
					degenGeom->getDegenPlate().zcamber[i][j],			\
					degenGeom->getDegenPlate().t[i][j],					\
					degenGeom->getDegenPlate().nCamber[i][j].x(),		\
					degenGeom->getDegenPlate().nCamber[i][j].y(),		\
					degenGeom->getDegenPlate().nCamber[i][j].z(),		\
					degenGeom->getDegenPlate().u[i],					\
					degenGeom->getDegenPlate().wTop[(num_pnts+1)/2+j],	\
					degenGeom->getDegenPlate().wBot[(num_pnts+1)/2+j]	);
		}
	}

	fprintf(file_id, "# DegenGeom Type, nXsecs\nSTICK, %d\n# xle,yle,zle,xte,yte,zte,xcg_solid,ycg_solid,zcg_solid,"
			"xcg_shell,ycg_shell,zcg_shell,toc,tLoc,chord,sweep,Ixx_shell_A,Ixx_shell_B,Izz_shell_A,"
			"Izz_shell_B,J_shell_A,J_shell_B,Ixx_solid,Izz_solid,J_solid,area,areaNormalX,"
			"areaNormalY,areaNormalZ,perimTop,perimBot,u\n", nxsecs);

	for ( int i = nxsecs; i < 2 * nxsecs; i++ )
	{
		fprintf(	file_id, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",	\
				degenGeom->getDegenStick().xle[i].x(),					\
				degenGeom->getDegenStick().xle[i].y(),					\
				degenGeom->getDegenStick().xle[i].z(),					\
				degenGeom->getDegenStick().xte[i].x(),					\
				degenGeom->getDegenStick().xte[i].y(),					\
				degenGeom->getDegenStick().xte[i].z(),					\
				degenGeom->getDegenStick().xcgSolid[i].x(),				\
				degenGeom->getDegenStick().xcgSolid[i].y(),				\
				degenGeom->getDegenStick().xcgSolid[i].z(),				\
				degenGeom->getDegenStick().xcgShell[i].x(),				\
				degenGeom->getDegenStick().xcgShell[i].y(),				\
				degenGeom->getDegenStick().xcgShell[i].z(),				\
				degenGeom->getDegenStick().toc[i],						\
				degenGeom->getDegenStick().tLoc[i],						\
				degenGeom->getDegenStick().chord[i],					\
				degenGeom->getDegenStick().sweep[i],					\
				degenGeom->getDegenStick().Ishell[i][0],				\
				degenGeom->getDegenStick().Ishell[i][1],				\
				degenGeom->getDegenStick().Ishell[i][2],				\
				degenGeom->getDegenStick().Ishell[i][3],				\
				degenGeom->getDegenStick().Ishell[i][4],				\
				degenGeom->getDegenStick().Ishell[i][5],				\
				degenGeom->getDegenStick().Isolid[i][0],				\
				degenGeom->getDegenStick().Isolid[i][1],				\
				degenGeom->getDegenStick().Isolid[i][2],				\
				degenGeom->getDegenStick().area[i],						\
				degenGeom->getDegenStick().areaNormal[i].x(),			\
				degenGeom->getDegenStick().areaNormal[i].y(),			\
				degenGeom->getDegenStick().areaNormal[i].z(),			\
				degenGeom->getDegenStick().perimTop[i],					\
				degenGeom->getDegenStick().perimBot[i],					\
				degenGeom->getDegenStick().u[i]							);
	}

	nxsecs = nxsecsOrig;
	fprintf(file_id, "# DegenGeom Type\n");
	fprintf(file_id, "POINT\n");
	fprintf(file_id, "# vol,volWet,area,areaWet,IxxShell,IyyShell,IzzShell,IxyShell,");
	fprintf(file_id, "IxzShell,IyzShell,IxxSolid,IyySolid,IzzSolid,IxySolid,IxzSolid,");
	fprintf(file_id, "IyzSolid,xcgShell,ycgShell,zcgShell,xcgSolid,ycgSolid,zcgSolid\n");
	fprintf(	file_id, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",\
			degenGeom->getDegenPoint().vol[1],			\
			degenGeom->getDegenPoint().volWet[1],		\
			degenGeom->getDegenPoint().area[1],			\
			degenGeom->getDegenPoint().areaWet[1],		\
			degenGeom->getDegenPoint().Ishell[1][0],	\
			degenGeom->getDegenPoint().Ishell[1][1],	\
			degenGeom->getDegenPoint().Ishell[1][2],	\
			degenGeom->getDegenPoint().Ishell[1][3],	\
			degenGeom->getDegenPoint().Ishell[1][4],	\
			degenGeom->getDegenPoint().Ishell[1][5],	\
			degenGeom->getDegenPoint().Isolid[1][0],	\
			degenGeom->getDegenPoint().Isolid[1][1],	\
			degenGeom->getDegenPoint().Isolid[1][2],	\
			degenGeom->getDegenPoint().Isolid[1][3],	\
			degenGeom->getDegenPoint().Isolid[1][4],	\
			degenGeom->getDegenPoint().Isolid[1][5],	\
			degenGeom->getDegenPoint().xcgShell[1].x(),	\
			degenGeom->getDegenPoint().xcgShell[1].y(),	\
			degenGeom->getDegenPoint().xcgShell[1].z(),	\
			degenGeom->getDegenPoint().xcgSolid[1].x(),	\
			degenGeom->getDegenPoint().xcgSolid[1].y(),	\
			degenGeom->getDegenPoint().xcgSolid[1].z()	);
}

void DegenGeom::write_degenGeomM_file(DegenGeom* degenGeom, FILE* file_id)
{
	int nxsecs = num_xsecs;

	if ( degenGeom->getParentGeom()->getTypeStr() == "wing" || degenGeom->getParentGeom()->getTypeStr() == "prop" )
		nxsecs -= 2;
	int nxsecsOrig = nxsecs;

	//============================= DegenSurf =============================//
	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).surf.sect(%d).X = [", (i+1));
		for ( int j = 0; j < num_pnts-1; j++ )
		{
			fprintf(	file_id, "%f, %f, %f;\n",						\
					degenGeom->getDegenSurface().x[i][j].x(),		\
					degenGeom->getDegenSurface().x[i][j].y(),		\
					degenGeom->getDegenSurface().x[i][j].z()		);
		}
		fprintf(	file_id, "%f, %f, %f];",							\
				degenGeom->getDegenSurface().x[i][num_pnts-1].x(),	\
				degenGeom->getDegenSurface().x[i][num_pnts-1].y(),	\
				degenGeom->getDegenSurface().x[i][num_pnts-1].z()	);
	}

	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).surf.sect(%d).Xn = [", (i+1));
		for ( int j = 0; j < num_pnts-1; j++ )
		{
			fprintf(	file_id, "%f, %f, %f;\n",						\
					degenGeom->getDegenSurface().nvec[i][j].x(),	\
					degenGeom->getDegenSurface().nvec[i][j].y(),	\
					degenGeom->getDegenSurface().nvec[i][j].z()		);
		}
		fprintf(	file_id, "%f, %f, %f];",								\
				degenGeom->getDegenSurface().nvec[i][num_pnts-1].x(),	\
				degenGeom->getDegenSurface().nvec[i][num_pnts-1].y(),	\
				degenGeom->getDegenSurface().nvec[i][num_pnts-1].z()	);
	}

	fprintf(file_id, "\ndegenGeom(end).surf.u = [");
	for ( int i = 0; i < nxsecs - 1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenSurface().u[i]);
	}
	fprintf(file_id, "%f];", degenGeom->getDegenSurface().u[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).surf.w = [");
	for ( int j = 0; j < num_pnts-1; j++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenSurface().w[j]);
	}
	fprintf(file_id, "%f];", degenGeom->getDegenSurface().w[num_pnts-1]);

	// Twice as many cross sections for BODY type
	if ( degenGeom->getType() == DegenGeom::BODY_TYPE ) nxsecs *= 2;
	//============================= DegenPlate =============================//
	fprintf(file_id, "\ndegenGeom(end).plate.nPlate = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf( file_id, "%f, %f, %f;\n",  degenGeom->getDegenPlate().nPlate[i].x(), \
				degenGeom->getDegenPlate().nPlate[i].y(), \
				degenGeom->getDegenPlate().nPlate[i].z()  );
	}
	fprintf(file_id, "%f, %f, %f];", degenGeom->getDegenPlate().nPlate[nxsecs-1].x(), \
			degenGeom->getDegenPlate().nPlate[nxsecs-1].y(), \
			degenGeom->getDegenPlate().nPlate[nxsecs-1].z()  );

	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).plate.sect(%d).X = [", (i+1));
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
		{
			fprintf(	file_id, "%f, %f, %f;\n",				\
					degenGeom->getDegenPlate().x[i][j].x(),	\
					degenGeom->getDegenPlate().x[i][j].y(),	\
					degenGeom->getDegenPlate().x[i][j].z()	);
		}
		fprintf(	file_id, "%f, %f, %f];\n",				\
				degenGeom->getDegenPlate().x[i][(num_pnts+1)/2-1].x(),	\
				degenGeom->getDegenPlate().x[i][(num_pnts+1)/2-1].y(),	\
				degenGeom->getDegenPlate().x[i][(num_pnts+1)/2-1].z()	);
	}

	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).plate.sect(%d).zCamber = [", (i+1));
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
		{
			fprintf( file_id, "%f, ", degenGeom->getDegenPlate().zcamber[i][j] );
		}
		fprintf(	file_id, "%f];\n", degenGeom->getDegenPlate().zcamber[i][(num_pnts+1)/2-1]);
	}

	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).plate.sect(%d).t = [", (i+1));
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
		{
			fprintf( file_id, "%f, ", degenGeom->getDegenPlate().t[i][j] );
		}
		fprintf(	file_id, "%f];\n", degenGeom->getDegenPlate().t[i][(num_pnts+1)/2-1]);
	}

	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).plate.sect(%d).nCamber = [", (i+1));
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
		{
			fprintf(	file_id, "%f, %f, %f;\n",						\
					degenGeom->getDegenPlate().nCamber[i][j].x(),	\
					degenGeom->getDegenPlate().nCamber[i][j].y(),	\
					degenGeom->getDegenPlate().nCamber[i][j].z()	);
		}
		fprintf(	file_id, "%f, %f, %f];\n",										\
				degenGeom->getDegenPlate().nCamber[i][(num_pnts+1)/2-1].x(),	\
				degenGeom->getDegenPlate().nCamber[i][(num_pnts+1)/2-1].y(),	\
				degenGeom->getDegenPlate().nCamber[i][(num_pnts+1)/2-1].z()		);
	}

	fprintf(file_id, "\ndegenGeom(end).plate.u = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenPlate().u[i]);
	}
	fprintf(file_id, "%f];", degenGeom->getDegenPlate().u[nxsecs-1]);

	int wCnt = (num_pnts+1)/2;
	if ( nxsecs != nxsecsOrig ) wCnt *= 2;
	fprintf(file_id, "\ndegenGeom(end).plate.wTop = [");
	for ( int j = 0; j < wCnt-1; j++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenPlate().wTop[j]);
	}
	fprintf(file_id, "%f];", degenGeom->getDegenPlate().wTop[wCnt-1]);

	fprintf(file_id, "\ndegenGeom(end).plate.wBot = [");
	for ( int j = 0; j < wCnt-1; j++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenPlate().wBot[j]);
	}
	fprintf(file_id, "%f];", degenGeom->getDegenPlate().wBot[wCnt-1]);

	//============================= DegenStick =============================//
	fprintf(file_id, "\ndegenGeom(end).stick.Xle = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenGeom->getDegenStick().xle[i].x(),	\
				degenGeom->getDegenStick().xle[i].y(),	\
				degenGeom->getDegenStick().xle[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenGeom->getDegenStick().xle[nxsecs-1].x(),	\
			degenGeom->getDegenStick().xle[nxsecs-1].y(),	\
			degenGeom->getDegenStick().xle[nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick.Xte = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenGeom->getDegenStick().xte[i].x(),	\
				degenGeom->getDegenStick().xte[i].y(),	\
				degenGeom->getDegenStick().xte[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenGeom->getDegenStick().xte[nxsecs-1].x(),	\
			degenGeom->getDegenStick().xte[nxsecs-1].y(),	\
			degenGeom->getDegenStick().xte[nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick.XcgSolid = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenGeom->getDegenStick().xcgSolid[i].x(),	\
				degenGeom->getDegenStick().xcgSolid[i].y(),	\
				degenGeom->getDegenStick().xcgSolid[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenGeom->getDegenStick().xcgSolid[nxsecs-1].x(),	\
			degenGeom->getDegenStick().xcgSolid[nxsecs-1].y(),	\
			degenGeom->getDegenStick().xcgSolid[nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick.XcgShell = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenGeom->getDegenStick().xcgShell[i].x(),	\
				degenGeom->getDegenStick().xcgShell[i].y(),	\
				degenGeom->getDegenStick().xcgShell[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenGeom->getDegenStick().xcgShell[nxsecs-1].x(),	\
			degenGeom->getDegenStick().xcgShell[nxsecs-1].y(),	\
			degenGeom->getDegenStick().xcgShell[nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick.toc = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenStick().toc[i]);
	}
	fprintf(file_id, "%f];\n", degenGeom->getDegenStick().toc[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.tLoc = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenStick().tLoc[i]);
	}
	fprintf(file_id, "%f];\n", degenGeom->getDegenStick().tLoc[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.chord = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenStick().chord[i]);
	}
	fprintf(file_id, "%f];\n", degenGeom->getDegenStick().chord[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.sweep = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenStick().sweep[i]);
	}
	fprintf(file_id, "%f];\n", degenGeom->getDegenStick().sweep[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.Ishell = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(	file_id, "%f, %f, %f, %f, %f, %f;\n",		\
				degenGeom->getDegenStick().Ishell[i][0],	\
				degenGeom->getDegenStick().Ishell[i][1],	\
				degenGeom->getDegenStick().Ishell[i][2],	\
				degenGeom->getDegenStick().Ishell[i][3],	\
				degenGeom->getDegenStick().Ishell[i][4],	\
				degenGeom->getDegenStick().Ishell[i][5]		);
	}
	fprintf(	file_id, "%f, %f, %f, %f, %f, %f];\n",			\
			degenGeom->getDegenStick().Ishell[nxsecs-1][0],	\
			degenGeom->getDegenStick().Ishell[nxsecs-1][1],	\
			degenGeom->getDegenStick().Ishell[nxsecs-1][2],	\
			degenGeom->getDegenStick().Ishell[nxsecs-1][3],	\
			degenGeom->getDegenStick().Ishell[nxsecs-1][4],	\
			degenGeom->getDegenStick().Ishell[nxsecs-1][5]	);

	fprintf(file_id, "\ndegenGeom(end).stick.Isolid = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(	file_id, "%f, %f, %f;\n",					\
				degenGeom->getDegenStick().Isolid[i][0],	\
				degenGeom->getDegenStick().Isolid[i][1],	\
				degenGeom->getDegenStick().Isolid[i][2]		);
	}
	fprintf(	file_id, "%f, %f, %f];\n",						\
			degenGeom->getDegenStick().Isolid[nxsecs-1][0],	\
			degenGeom->getDegenStick().Isolid[nxsecs-1][1],	\
			degenGeom->getDegenStick().Isolid[nxsecs-1][2]	);

	fprintf(file_id, "\ndegenGeom(end).stick.area = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenStick().area[i]);
	}
	fprintf(file_id, "%f];\n", degenGeom->getDegenStick().area[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.areaNormal = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(	file_id, "%f, %f, %f;\n",						\
				degenGeom->getDegenStick().areaNormal[i].x(),	\
				degenGeom->getDegenStick().areaNormal[i].y(),	\
				degenGeom->getDegenStick().areaNormal[i].z()	);
	}
	fprintf(	file_id, "%f, %f, %f];\n",							\
			degenGeom->getDegenStick().areaNormal[nxsecs-1].x(),\
			degenGeom->getDegenStick().areaNormal[nxsecs-1].y(),\
			degenGeom->getDegenStick().areaNormal[nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick.perimTop = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenStick().perimTop[i]);
	}
	fprintf(file_id, "%f];\n", degenGeom->getDegenStick().perimTop[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.perimBot = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenStick().perimBot[i]);
	}
	fprintf(file_id, "%f];\n", degenGeom->getDegenStick().perimBot[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.u = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenStick().u[i]);
	}
	fprintf(file_id, "%f];\n", degenGeom->getDegenStick().u[nxsecs-1]);

	nxsecs = nxsecsOrig;
	//============================= DegenPoint =============================//
	fprintf(file_id, "\ndegenGeom(end).point.vol = %f;\n", degenGeom->getDegenPoint().vol[0]);
	fprintf(file_id, "\ndegenGeom(end).point.volWet = %f;\n", degenGeom->getDegenPoint().volWet[0]);
	fprintf(file_id, "\ndegenGeom(end).point.area = %f;\n", degenGeom->getDegenPoint().area[0]);
	fprintf(file_id, "\ndegenGeom(end).point.areaWet = %f;\n", degenGeom->getDegenPoint().areaWet[0]);
	fprintf(file_id, "\ndegenGeom(end).point.Ishell = [%f, %f, %f, %f, %f, %f];\n",	\
			degenGeom->getDegenPoint().Ishell[0][0],	\
			degenGeom->getDegenPoint().Ishell[0][1],	\
			degenGeom->getDegenPoint().Ishell[0][2],	\
			degenGeom->getDegenPoint().Ishell[0][3],	\
			degenGeom->getDegenPoint().Ishell[0][4],	\
			degenGeom->getDegenPoint().Ishell[0][5]		);
	fprintf(file_id, "\ndegenGeom(end).point.Isolid = [%f, %f, %f, %f, %f, %f];\n",	\
			degenGeom->getDegenPoint().Isolid[0][0],	\
			degenGeom->getDegenPoint().Isolid[0][1],	\
			degenGeom->getDegenPoint().Isolid[0][2],	\
			degenGeom->getDegenPoint().Isolid[0][3],	\
			degenGeom->getDegenPoint().Isolid[0][4],	\
			degenGeom->getDegenPoint().Isolid[0][5]		);
	fprintf(file_id, "\ndegenGeom(end).point.xcgShell = [%f, %f, %f];\n",	\
			degenGeom->getDegenPoint().xcgShell[0].x(),	\
			degenGeom->getDegenPoint().xcgShell[0].y(),	\
			degenGeom->getDegenPoint().xcgShell[0].z()	);
	fprintf(file_id, "\ndegenGeom(end).point.xcgSolid = [%f, %f, %f];\n",	\
			degenGeom->getDegenPoint().xcgSolid[0].x(),	\
			degenGeom->getDegenPoint().xcgSolid[0].y(),	\
			degenGeom->getDegenPoint().xcgSolid[0].z()	);
}

void DegenGeom::write_refl_degenGeomM_file(DegenGeom* degenGeom, FILE* file_id)
{
	int nxsecs = num_xsecs;

	if ( degenGeom->getParentGeom()->getTypeStr() == "wing" || degenGeom->getParentGeom()->getTypeStr() == "prop" )
		nxsecs -= 2;
	int nxsecsOrig = nxsecs;

	//============================= DegenSurf =============================//
	for ( int i = nxsecs; i < 2*nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).surf.sect(%d).X = [", (i-nxsecs+1));
		for ( int j = 0; j < num_pnts-1; j++ )
		{
			fprintf(	file_id, "%f, %f, %f;\n",						\
					degenGeom->getDegenSurface().x[i][j].x(),		\
					degenGeom->getDegenSurface().x[i][j].y(),		\
					degenGeom->getDegenSurface().x[i][j].z()		);
		}
		fprintf(	file_id, "%f, %f, %f];",							\
				degenGeom->getDegenSurface().x[i][num_pnts-1].x(),	\
				degenGeom->getDegenSurface().x[i][num_pnts-1].y(),	\
				degenGeom->getDegenSurface().x[i][num_pnts-1].z()	);
	}

	for ( int i = nxsecs; i < 2*nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).surf.sect(%d).Xn = [", (i-nxsecs+1));
		for ( int j = 0; j < num_pnts-1; j++ )
		{
			fprintf(	file_id, "%f, %f, %f;\n",						\
					degenGeom->getDegenSurface().nvec[i][j].x(),	\
					degenGeom->getDegenSurface().nvec[i][j].y(),	\
					degenGeom->getDegenSurface().nvec[i][j].z()		);
		}
		fprintf(	file_id, "%f, %f, %f];",								\
				degenGeom->getDegenSurface().nvec[i][num_pnts-1].x(),	\
				degenGeom->getDegenSurface().nvec[i][num_pnts-1].y(),	\
				degenGeom->getDegenSurface().nvec[i][num_pnts-1].z()	);
	}

	fprintf(file_id, "\ndegenGeom(end).surf.u = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenSurface().u[i]);
	}
	fprintf(file_id, "%f];", degenGeom->getDegenSurface().u[2*nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).surf.w = [");
	for ( int j = 0; j < num_pnts-1; j++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenSurface().w[j]);
	}
	fprintf(file_id, "%f];", degenGeom->getDegenSurface().w[num_pnts-1]);

	// Twice as many cross sections for BODY type
	if ( degenGeom->getType() == DegenGeom::BODY_TYPE ) nxsecs *= 2;
	//============================= DegenPlate =============================//
	fprintf(file_id, "\ndegenGeom(end).plate.nPlate = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf( file_id, "%f, %f, %f;\n",  degenGeom->getDegenPlate().nPlate[i].x(), \
				degenGeom->getDegenPlate().nPlate[i].y(), \
				degenGeom->getDegenPlate().nPlate[i].z()  );
	}
	fprintf(file_id, "%f, %f, %f];", degenGeom->getDegenPlate().nPlate[2*nxsecs-1].x(), \
			degenGeom->getDegenPlate().nPlate[2*nxsecs-1].y(), \
			degenGeom->getDegenPlate().nPlate[2*nxsecs-1].z()  );

	for ( int i = nxsecs; i < 2*nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).plate.sect(%d).X = [", (i-nxsecs+1));
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
		{
			fprintf(	file_id, "%f, %f, %f;\n",				\
					degenGeom->getDegenPlate().x[i][j].x(),	\
					degenGeom->getDegenPlate().x[i][j].y(),	\
					degenGeom->getDegenPlate().x[i][j].z()	);
		}
		fprintf(	file_id, "%f, %f, %f];\n",				\
				degenGeom->getDegenPlate().x[i][(num_pnts+1)/2-1].x(),	\
				degenGeom->getDegenPlate().x[i][(num_pnts+1)/2-1].y(),	\
				degenGeom->getDegenPlate().x[i][(num_pnts+1)/2-1].z()	);
	}

	for ( int i = nxsecs; i < 2*nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).plate.sect(%d).zCamber = [", (i-nxsecs+1));
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
		{
			fprintf( file_id, "%f, ", degenGeom->getDegenPlate().zcamber[i][j] );
		}
		fprintf(	file_id, "%f];\n", degenGeom->getDegenPlate().zcamber[i][(num_pnts+1)/2-1]);
	}

	for ( int i = nxsecs; i < 2*nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).plate.sect(%d).t = [", (i-nxsecs+1));
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
		{
			fprintf( file_id, "%f, ", degenGeom->getDegenPlate().t[i][j] );
		}
		fprintf(	file_id, "%f];\n", degenGeom->getDegenPlate().t[i][(num_pnts+1)/2-1]);
	}

	for ( int i = nxsecs; i < 2*nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).plate.sect(%d).nCamber = [", (i-nxsecs+1));
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
		{
			fprintf(	file_id, "%f, %f, %f;\n",						\
					degenGeom->getDegenPlate().nCamber[i][j].x(),	\
					degenGeom->getDegenPlate().nCamber[i][j].y(),	\
					degenGeom->getDegenPlate().nCamber[i][j].z()	);
		}
		fprintf(	file_id, "%f, %f, %f];\n",										\
				degenGeom->getDegenPlate().nCamber[i][(num_pnts+1)/2-1].x(),	\
				degenGeom->getDegenPlate().nCamber[i][(num_pnts+1)/2-1].y(),	\
				degenGeom->getDegenPlate().nCamber[i][(num_pnts+1)/2-1].z()		);
	}

	fprintf(file_id, "\ndegenGeom(end).plate.u = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenPlate().u[i]);
	}
	fprintf(file_id, "%f];", degenGeom->getDegenPlate().u[2*nxsecs-1]);

	int wCnt = (num_pnts+1)/2;
	if (nxsecs != nxsecsOrig ) wCnt *= 2;
	fprintf(file_id, "\ndegenGeom(end).plate.wTop = [");
	for ( int j = wCnt; j < 2*wCnt-1; j++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenPlate().wTop[j]);
	}
	fprintf(file_id, "%f];", degenGeom->getDegenPlate().wTop[2*wCnt-1]);

	fprintf(file_id, "\ndegenGeom(end).plate.wBot = [");
	for ( int j = wCnt; j < 2*wCnt-1; j++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenPlate().wBot[j]);
	}
	fprintf(file_id, "%f];", degenGeom->getDegenPlate().wBot[2*wCnt-1]);

	//============================= DegenStick =============================//
	fprintf(file_id, "\ndegenGeom(end).stick.Xle = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenGeom->getDegenStick().xle[i].x(),	\
				degenGeom->getDegenStick().xle[i].y(),	\
				degenGeom->getDegenStick().xle[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenGeom->getDegenStick().xle[2*nxsecs-1].x(),	\
			degenGeom->getDegenStick().xle[2*nxsecs-1].y(),	\
			degenGeom->getDegenStick().xle[2*nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick.Xte = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenGeom->getDegenStick().xte[i].x(),	\
				degenGeom->getDegenStick().xte[i].y(),	\
				degenGeom->getDegenStick().xte[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenGeom->getDegenStick().xte[2*nxsecs-1].x(),	\
			degenGeom->getDegenStick().xte[2*nxsecs-1].y(),	\
			degenGeom->getDegenStick().xte[2*nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick.XcgSolid = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenGeom->getDegenStick().xcgSolid[i].x(),	\
				degenGeom->getDegenStick().xcgSolid[i].y(),	\
				degenGeom->getDegenStick().xcgSolid[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenGeom->getDegenStick().xcgSolid[2*nxsecs-1].x(),	\
			degenGeom->getDegenStick().xcgSolid[2*nxsecs-1].y(),	\
			degenGeom->getDegenStick().xcgSolid[2*nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick.XcgShell = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenGeom->getDegenStick().xcgShell[i].x(),	\
				degenGeom->getDegenStick().xcgShell[i].y(),	\
				degenGeom->getDegenStick().xcgShell[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenGeom->getDegenStick().xcgShell[2*nxsecs-1].x(),	\
			degenGeom->getDegenStick().xcgShell[2*nxsecs-1].y(),	\
			degenGeom->getDegenStick().xcgShell[2*nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick.toc = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenStick().toc[i]);
	}
	fprintf(file_id, "%f];\n", degenGeom->getDegenStick().toc[2*nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.tLoc = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenStick().tLoc[i]);
	}
	fprintf(file_id, "%f];\n", degenGeom->getDegenStick().tLoc[2*nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.chord = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenStick().chord[i]);
	}
	fprintf(file_id, "%f];\n", degenGeom->getDegenStick().chord[2*nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.sweep = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenStick().sweep[i]);
	}
	fprintf(file_id, "%f];\n", degenGeom->getDegenStick().sweep[2*nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.Ishell = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(	file_id, "%f, %f, %f, %f, %f, %f;\n",		\
				degenGeom->getDegenStick().Ishell[i][0],	\
				degenGeom->getDegenStick().Ishell[i][1],	\
				degenGeom->getDegenStick().Ishell[i][2],	\
				degenGeom->getDegenStick().Ishell[i][3],	\
				degenGeom->getDegenStick().Ishell[i][4],	\
				degenGeom->getDegenStick().Ishell[i][5]		);
	}
	fprintf(	file_id, "%f, %f, %f, %f, %f, %f];\n",			\
			degenGeom->getDegenStick().Ishell[2*nxsecs-1][0],	\
			degenGeom->getDegenStick().Ishell[2*nxsecs-1][1],	\
			degenGeom->getDegenStick().Ishell[2*nxsecs-1][2],	\
			degenGeom->getDegenStick().Ishell[2*nxsecs-1][3],	\
			degenGeom->getDegenStick().Ishell[2*nxsecs-1][4],	\
			degenGeom->getDegenStick().Ishell[2*nxsecs-1][5]	);

	fprintf(file_id, "\ndegenGeom(end).stick.Isolid = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(	file_id, "%f, %f, %f;\n",					\
				degenGeom->getDegenStick().Isolid[i][0],	\
				degenGeom->getDegenStick().Isolid[i][1],	\
				degenGeom->getDegenStick().Isolid[i][2]		);
	}
	fprintf(	file_id, "%f, %f, %f];\n",						\
			degenGeom->getDegenStick().Isolid[2*nxsecs-1][0],	\
			degenGeom->getDegenStick().Isolid[2*nxsecs-1][1],	\
			degenGeom->getDegenStick().Isolid[2*nxsecs-1][2]	);

	fprintf(file_id, "\ndegenGeom(end).stick.area = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenStick().area[i]);
	}
	fprintf(file_id, "%f];\n", degenGeom->getDegenStick().area[2*nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.areaNormal = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(	file_id, "%f, %f, %f;\n",						\
				degenGeom->getDegenStick().areaNormal[i].x(),	\
				degenGeom->getDegenStick().areaNormal[i].y(),	\
				degenGeom->getDegenStick().areaNormal[i].z()	);
	}
	fprintf(	file_id, "%f, %f, %f];\n",							\
			degenGeom->getDegenStick().areaNormal[2*nxsecs-1].x(),\
			degenGeom->getDegenStick().areaNormal[2*nxsecs-1].y(),\
			degenGeom->getDegenStick().areaNormal[2*nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick.perimTop = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenStick().perimTop[i]);
	}
	fprintf(file_id, "%f];\n", degenGeom->getDegenStick().perimTop[2*nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.perimBot = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenStick().perimBot[i]);
	}
	fprintf(file_id, "%f];\n", degenGeom->getDegenStick().perimBot[2*nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.u = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenGeom->getDegenStick().u[i]);
	}
	fprintf(file_id, "%f];\n", degenGeom->getDegenStick().u[2*nxsecs-1]);

	nxsecs = nxsecsOrig;
	//============================= DegenPoint =============================//
	fprintf(file_id, "\ndegenGeom(end).point.vol = %f;\n", degenGeom->getDegenPoint().vol[1]);
	fprintf(file_id, "\ndegenGeom(end).point.volWet = %f;\n", degenGeom->getDegenPoint().volWet[1]);
	fprintf(file_id, "\ndegenGeom(end).point.area = %f;\n", degenGeom->getDegenPoint().area[1]);
	fprintf(file_id, "\ndegenGeom(end).point.areaWet = %f;\n", degenGeom->getDegenPoint().areaWet[1]);
	fprintf(file_id, "\ndegenGeom(end).point.Ishell = [%f, %f, %f, %f, %f, %f];\n",	\
			degenGeom->getDegenPoint().Ishell[1][0],	\
			degenGeom->getDegenPoint().Ishell[1][1],	\
			degenGeom->getDegenPoint().Ishell[1][2],	\
			degenGeom->getDegenPoint().Ishell[1][3],	\
			degenGeom->getDegenPoint().Ishell[1][4],	\
			degenGeom->getDegenPoint().Ishell[1][5]		);
	fprintf(file_id, "\ndegenGeom(end).point.Isolid = [%f, %f, %f, %f, %f, %f];\n",	\
			degenGeom->getDegenPoint().Isolid[1][0],	\
			degenGeom->getDegenPoint().Isolid[1][1],	\
			degenGeom->getDegenPoint().Isolid[1][2],	\
			degenGeom->getDegenPoint().Isolid[1][3],	\
			degenGeom->getDegenPoint().Isolid[1][4],	\
			degenGeom->getDegenPoint().Isolid[1][5]		);
	fprintf(file_id, "\ndegenGeom(end).point.xcgShell = [%f, %f, %f];\n",	\
			degenGeom->getDegenPoint().xcgShell[1].x(),	\
			degenGeom->getDegenPoint().xcgShell[1].y(),	\
			degenGeom->getDegenPoint().xcgShell[1].z()	);
	fprintf(file_id, "\ndegenGeom(end).point.xcgSolid = [%f, %f, %f];\n",	\
			degenGeom->getDegenPoint().xcgSolid[1].x(),	\
			degenGeom->getDegenPoint().xcgSolid[1].y(),	\
			degenGeom->getDegenPoint().xcgSolid[1].z()	);
}
