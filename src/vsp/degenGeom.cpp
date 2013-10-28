#include "degenGeom.h"
#include "geom.h"
#include <cmath>

void DegenGeom::write_degenGeomCsv_file(FILE* file_id)
{
	int nxsecs = num_xsecs;

	if ( parentGeom->getTypeStr() == "wing" || parentGeom->getTypeStr() == "prop" )
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
					degenSurface.x[i][j].x(),		\
					degenSurface.x[i][j].y(),		\
					degenSurface.x[i][j].z(),		\
					degenSurface.nvec[i][j].x(),	\
					degenSurface.nvec[i][j].y(),	\
					degenSurface.nvec[i][j].z(),	\
					degenSurface.u[i],				\
					degenSurface.w[j]				);
		}
	}

	// JBB: Twice as many cross sections for BODY type
	if ( type == DegenGeom::BODY_TYPE ) nxsecs *= 2;

	fprintf(file_id, "# DegenGeom Type,nXsecs,nPnts/Xsec\n");
	fprintf(file_id, "PLATE,%d,%d\n", nxsecs, (num_pnts+1)/2);
	fprintf(file_id,"# xn,yn,zn\n");
	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf( file_id, "%f,%f,%f\n", degenPlate.nPlate[i].x(), \
				degenPlate.nPlate[i].y(), \
				degenPlate.nPlate[i].z()  );
	}

	fprintf(file_id, "# x,y,z,zCamb,t,nCamrX,nCambY,nCambZ,u,wTop,wBot\n");
	for ( int i = 0; i < nxsecsOrig; i++ )
	{
		for ( int j = 0; j < (num_pnts+1)/2; j++ )
		{
			fprintf(	file_id, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",	\
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
	for ( int i = nxsecsOrig; i < nxsecs; i++ )
	{
		for ( int j = 0; j < (num_pnts+1)/2; j++ )
		{
			fprintf(	file_id, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",	\
					degenPlate.x[i][j].x(),				\
					degenPlate.x[i][j].y(),				\
					degenPlate.x[i][j].z(),				\
					degenPlate.zcamber[i][j],			\
					degenPlate.t[i][j],					\
					degenPlate.nCamber[i][j].x(),		\
					degenPlate.nCamber[i][j].y(),		\
					degenPlate.nCamber[i][j].z(),		\
					degenPlate.u[i],					\
					degenPlate.wTop[(num_pnts+1)/2+j],	\
					degenPlate.wBot[(num_pnts+1)/2+j]	);
		}
	}

	fprintf(file_id, "# DegenGeom Type, nXsecs\nSTICK, %d\n# xle,yle,zle,xte,yte,zte,xcg_solid,ycg_solid,zcg_solid,"
			"xcg_shell,ycg_shell,zcg_shell,toc,tLoc,chord,sweep,Ixx_shell_A,Ixx_shell_B,Izz_shell_A,"
			"Izz_shell_B,J_shell_A,J_shell_B,Ixx_solid,Izz_solid,J_solid,area,areaNormalX,"
			"areaNormalY,areaNormalZ,perimTop,perimBot,u\n", nxsecs);

	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf(	file_id, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",	\
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

	nxsecs = nxsecsOrig;
	fprintf(file_id, "# DegenGeom Type\n");
	fprintf(file_id, "POINT\n");
	fprintf(file_id, "# vol,volWet,area,areaWet,IxxShell,IyyShell,IzzShell,IxyShell,");
	fprintf(file_id, "IxzShell,IyzShell,IxxSolid,IyySolid,IzzSolid,IxySolid,IxzSolid,");
	fprintf(file_id, "IyzSolid,xcgShell,ycgShell,zcgShell,xcgSolid,ycgSolid,zcgSolid\n");
	fprintf(	file_id, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",\
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

void DegenGeom::write_refl_degenGeomCsv_file(FILE* file_id)
{
	int nxsecs = num_xsecs;

	if ( parentGeom->getTypeStr() == "wing" || parentGeom->getTypeStr() == "prop" )
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
					degenSurface.x[i][j].x(),		\
					degenSurface.x[i][j].y(),		\
					degenSurface.x[i][j].z(),		\
					degenSurface.nvec[i][j].x(),	\
					degenSurface.nvec[i][j].y(),	\
					degenSurface.nvec[i][j].z(),	\
					degenSurface.u[i],				\
					degenSurface.w[j]				);
		}
	}

	// Twice as many cross sections for BODY type
	if ( type == DegenGeom::BODY_TYPE ) nxsecs *= 2;
	fprintf(file_id, "# DegenGeom Type,nXsecs,nPnts/Xsec\n");
	fprintf(file_id, "PLATE,%d,%d\n", nxsecs, (num_pnts+1)/2);
	fprintf(file_id,"# xn,yn,zn\n");
	for ( int i = nxsecs; i < 2 * nxsecs; i++ )
	{
		fprintf( file_id, "%f,%f,%f\n", degenPlate.nPlate[i].x(), \
				degenPlate.nPlate[i].y(), \
				degenPlate.nPlate[i].z()  );
	}

	fprintf(file_id, "# x,y,z,zCamb,t,nCambX,nCambY,nCambZ,u,wTop,wBot\n");
	for ( int i = nxsecs; i < nxsecsOrig + nxsecs; i++ )
	{
		for ( int j = 0; j < (num_pnts+1)/2; j++ )
		{
			fprintf(	file_id, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",	\
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
	for ( int i = nxsecsOrig + nxsecs; i < 2*nxsecs; i++ )
	{
		for ( int j = 0; j < (num_pnts+1)/2; j++ )
		{
			fprintf(	file_id, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",	\
					degenPlate.x[i][j].x(),				\
					degenPlate.x[i][j].y(),				\
					degenPlate.x[i][j].z(),				\
					degenPlate.zcamber[i][j],			\
					degenPlate.t[i][j],					\
					degenPlate.nCamber[i][j].x(),		\
					degenPlate.nCamber[i][j].y(),		\
					degenPlate.nCamber[i][j].z(),		\
					degenPlate.u[i],					\
					degenPlate.wTop[(num_pnts+1)/2+j],	\
					degenPlate.wBot[(num_pnts+1)/2+j]	);
		}
	}

	fprintf(file_id, "# DegenGeom Type, nXsecs\nSTICK, %d\n# xle,yle,zle,xte,yte,zte,xcg_solid,ycg_solid,zcg_solid,"
			"xcg_shell,ycg_shell,zcg_shell,toc,tLoc,chord,sweep,Ixx_shell_A,Ixx_shell_B,Izz_shell_A,"
			"Izz_shell_B,J_shell_A,J_shell_B,Ixx_solid,Izz_solid,J_solid,area,areaNormalX,"
			"areaNormalY,areaNormalZ,perimTop,perimBot,u\n", nxsecs);

	for ( int i = nxsecs; i < 2 * nxsecs; i++ )
	{
		fprintf(	file_id, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",	\
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

	nxsecs = nxsecsOrig;
	fprintf(file_id, "# DegenGeom Type\n");
	fprintf(file_id, "POINT\n");
	fprintf(file_id, "# vol,volWet,area,areaWet,IxxShell,IyyShell,IzzShell,IxyShell,");
	fprintf(file_id, "IxzShell,IyzShell,IxxSolid,IyySolid,IzzSolid,IxySolid,IxzSolid,");
	fprintf(file_id, "IyzSolid,xcgShell,ycgShell,zcgShell,xcgSolid,ycgSolid,zcgSolid\n");
	fprintf(	file_id, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",\
			degenPoint.vol[1],			\
			degenPoint.volWet[1],		\
			degenPoint.area[1],			\
			degenPoint.areaWet[1],		\
			degenPoint.Ishell[1][0],	\
			degenPoint.Ishell[1][1],	\
			degenPoint.Ishell[1][2],	\
			degenPoint.Ishell[1][3],	\
			degenPoint.Ishell[1][4],	\
			degenPoint.Ishell[1][5],	\
			degenPoint.Isolid[1][0],	\
			degenPoint.Isolid[1][1],	\
			degenPoint.Isolid[1][2],	\
			degenPoint.Isolid[1][3],	\
			degenPoint.Isolid[1][4],	\
			degenPoint.Isolid[1][5],	\
			degenPoint.xcgShell[1].x(),	\
			degenPoint.xcgShell[1].y(),	\
			degenPoint.xcgShell[1].z(),	\
			degenPoint.xcgSolid[1].x(),	\
			degenPoint.xcgSolid[1].y(),	\
			degenPoint.xcgSolid[1].z()	);
}

void DegenGeom::write_degenGeomM_file(FILE* file_id)
{
	int nxsecs = num_xsecs;

	if ( parentGeom->getTypeStr() == "wing" || parentGeom->getTypeStr() == "prop" )
		nxsecs -= 2;
	int nxsecsOrig = nxsecs;

	//============================= DegenSurf =============================//
	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).surf.sect(%d).X = [", (i+1));
		for ( int j = 0; j < num_pnts-1; j++ )
		{
			fprintf(	file_id, "%f, %f, %f;\n",						\
					degenSurface.x[i][j].x(),		\
					degenSurface.x[i][j].y(),		\
					degenSurface.x[i][j].z()		);
		}
		fprintf(	file_id, "%f, %f, %f];",							\
				degenSurface.x[i][num_pnts-1].x(),	\
				degenSurface.x[i][num_pnts-1].y(),	\
				degenSurface.x[i][num_pnts-1].z()	);
	}

	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).surf.sect(%d).Xn = [", (i+1));
		for ( int j = 0; j < num_pnts-1; j++ )
		{
			fprintf(	file_id, "%f, %f, %f;\n",						\
					degenSurface.nvec[i][j].x(),	\
					degenSurface.nvec[i][j].y(),	\
					degenSurface.nvec[i][j].z()		);
		}
		fprintf(	file_id, "%f, %f, %f];",								\
				degenSurface.nvec[i][num_pnts-1].x(),	\
				degenSurface.nvec[i][num_pnts-1].y(),	\
				degenSurface.nvec[i][num_pnts-1].z()	);
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

	// Twice as many cross sections for BODY type
	if ( type == DegenGeom::BODY_TYPE ) nxsecs *= 2;
	//============================= DegenPlate =============================//
	fprintf(file_id, "\ndegenGeom(end).plate.nPlate = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf( file_id, "%f, %f, %f;\n",  degenPlate.nPlate[i].x(), \
				degenPlate.nPlate[i].y(), \
				degenPlate.nPlate[i].z()  );
	}
	fprintf(file_id, "%f, %f, %f];", degenPlate.nPlate[nxsecs-1].x(), \
			degenPlate.nPlate[nxsecs-1].y(), \
			degenPlate.nPlate[nxsecs-1].z()  );

	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).plate.sect(%d).X = [", (i+1));
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
		{
			fprintf(	file_id, "%f, %f, %f;\n",				\
					degenPlate.x[i][j].x(),	\
					degenPlate.x[i][j].y(),	\
					degenPlate.x[i][j].z()	);
		}
		fprintf(	file_id, "%f, %f, %f];\n",				\
				degenPlate.x[i][(num_pnts+1)/2-1].x(),	\
				degenPlate.x[i][(num_pnts+1)/2-1].y(),	\
				degenPlate.x[i][(num_pnts+1)/2-1].z()	);
	}

	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).plate.sect(%d).zCamber = [", (i+1));
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
		{
			fprintf( file_id, "%f, ", degenPlate.zcamber[i][j] );
		}
		fprintf(	file_id, "%f];\n", degenPlate.zcamber[i][(num_pnts+1)/2-1]);
	}

	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).plate.sect(%d).t = [", (i+1));
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
		{
			fprintf( file_id, "%f, ", degenPlate.t[i][j] );
		}
		fprintf(	file_id, "%f];\n", degenPlate.t[i][(num_pnts+1)/2-1]);
	}

	for ( int i = 0; i < nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).plate.sect(%d).nCamber = [", (i+1));
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
		{
			fprintf(	file_id, "%f, %f, %f;\n",						\
					degenPlate.nCamber[i][j].x(),	\
					degenPlate.nCamber[i][j].y(),	\
					degenPlate.nCamber[i][j].z()	);
		}
		fprintf(	file_id, "%f, %f, %f];\n",										\
				degenPlate.nCamber[i][(num_pnts+1)/2-1].x(),	\
				degenPlate.nCamber[i][(num_pnts+1)/2-1].y(),	\
				degenPlate.nCamber[i][(num_pnts+1)/2-1].z()		);
	}

	fprintf(file_id, "\ndegenGeom(end).plate.u = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenPlate.u[i]);
	}
	fprintf(file_id, "%f];", degenPlate.u[nxsecs-1]);

	int wCnt = (num_pnts+1)/2;
	if ( nxsecs != nxsecsOrig ) wCnt *= 2;
	fprintf(file_id, "\ndegenGeom(end).plate.wTop = [");
	for ( int j = 0; j < wCnt-1; j++ )
	{
		fprintf(file_id, "%f, ", degenPlate.wTop[j]);
	}
	fprintf(file_id, "%f];", degenPlate.wTop[wCnt-1]);

	fprintf(file_id, "\ndegenGeom(end).plate.wBot = [");
	for ( int j = 0; j < wCnt-1; j++ )
	{
		fprintf(file_id, "%f, ", degenPlate.wBot[j]);
	}
	fprintf(file_id, "%f];", degenPlate.wBot[wCnt-1]);

	//============================= DegenStick =============================//
	fprintf(file_id, "\ndegenGeom(end).stick.Xle = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenStick.xle[i].x(),	\
				degenStick.xle[i].y(),	\
				degenStick.xle[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenStick.xle[nxsecs-1].x(),	\
			degenStick.xle[nxsecs-1].y(),	\
			degenStick.xle[nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick.Xte = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenStick.xte[i].x(),	\
				degenStick.xte[i].y(),	\
				degenStick.xte[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenStick.xte[nxsecs-1].x(),	\
			degenStick.xte[nxsecs-1].y(),	\
			degenStick.xte[nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick.XcgSolid = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenStick.xcgSolid[i].x(),	\
				degenStick.xcgSolid[i].y(),	\
				degenStick.xcgSolid[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenStick.xcgSolid[nxsecs-1].x(),	\
			degenStick.xcgSolid[nxsecs-1].y(),	\
			degenStick.xcgSolid[nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick.XcgShell = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenStick.xcgShell[i].x(),	\
				degenStick.xcgShell[i].y(),	\
				degenStick.xcgShell[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenStick.xcgShell[nxsecs-1].x(),	\
			degenStick.xcgShell[nxsecs-1].y(),	\
			degenStick.xcgShell[nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick.toc = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.toc[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.toc[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.tLoc = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.tLoc[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.tLoc[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.chord = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.chord[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.chord[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.sweep = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.sweep[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.sweep[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.Ishell = [");
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

	fprintf(file_id, "\ndegenGeom(end).stick.Isolid = [");
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

	fprintf(file_id, "\ndegenGeom(end).stick.area = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.area[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.area[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.areaNormal = [");
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

	fprintf(file_id, "\ndegenGeom(end).stick.perimTop = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.perimTop[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.perimTop[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.perimBot = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.perimBot[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.perimBot[nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.u = [");
	for ( int i = 0; i < nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.u[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.u[nxsecs-1]);

	nxsecs = nxsecsOrig;
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

void DegenGeom::write_refl_degenGeomM_file(FILE* file_id)
{
	int nxsecs = num_xsecs;

	if ( parentGeom->getTypeStr() == "wing" || parentGeom->getTypeStr() == "prop" )
		nxsecs -= 2;
	int nxsecsOrig = nxsecs;

	//============================= DegenSurf =============================//
	for ( int i = nxsecs; i < 2*nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).surf.sect(%d).X = [", (i-nxsecs+1));
		for ( int j = 0; j < num_pnts-1; j++ )
		{
			fprintf(	file_id, "%f, %f, %f;\n",						\
					degenSurface.x[i][j].x(),		\
					degenSurface.x[i][j].y(),		\
					degenSurface.x[i][j].z()		);
		}
		fprintf(	file_id, "%f, %f, %f];",							\
				degenSurface.x[i][num_pnts-1].x(),	\
				degenSurface.x[i][num_pnts-1].y(),	\
				degenSurface.x[i][num_pnts-1].z()	);
	}

	for ( int i = nxsecs; i < 2*nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).surf.sect(%d).Xn = [", (i-nxsecs+1));
		for ( int j = 0; j < num_pnts-1; j++ )
		{
			fprintf(	file_id, "%f, %f, %f;\n",						\
					degenSurface.nvec[i][j].x(),	\
					degenSurface.nvec[i][j].y(),	\
					degenSurface.nvec[i][j].z()		);
		}
		fprintf(	file_id, "%f, %f, %f];",								\
				degenSurface.nvec[i][num_pnts-1].x(),	\
				degenSurface.nvec[i][num_pnts-1].y(),	\
				degenSurface.nvec[i][num_pnts-1].z()	);
	}

	fprintf(file_id, "\ndegenGeom(end).surf.u = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenSurface.u[i]);
	}
	fprintf(file_id, "%f];", degenSurface.u[2*nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).surf.w = [");
	for ( int j = 0; j < num_pnts-1; j++ )
	{
		fprintf(file_id, "%f, ", degenSurface.w[j]);
	}
	fprintf(file_id, "%f];", degenSurface.w[num_pnts-1]);

	// Twice as many cross sections for BODY type
	if ( type == DegenGeom::BODY_TYPE ) nxsecs *= 2;
	//============================= DegenPlate =============================//
	fprintf(file_id, "\ndegenGeom(end).plate.nPlate = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf( file_id, "%f, %f, %f;\n",  degenPlate.nPlate[i].x(), \
				degenPlate.nPlate[i].y(), \
				degenPlate.nPlate[i].z()  );
	}
	fprintf(file_id, "%f, %f, %f];", degenPlate.nPlate[2*nxsecs-1].x(), \
			degenPlate.nPlate[2*nxsecs-1].y(), \
			degenPlate.nPlate[2*nxsecs-1].z()  );

	for ( int i = nxsecs; i < 2*nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).plate.sect(%d).X = [", (i-nxsecs+1));
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
		{
			fprintf(	file_id, "%f, %f, %f;\n",				\
					degenPlate.x[i][j].x(),	\
					degenPlate.x[i][j].y(),	\
					degenPlate.x[i][j].z()	);
		}
		fprintf(	file_id, "%f, %f, %f];\n",				\
				degenPlate.x[i][(num_pnts+1)/2-1].x(),	\
				degenPlate.x[i][(num_pnts+1)/2-1].y(),	\
				degenPlate.x[i][(num_pnts+1)/2-1].z()	);
	}

	for ( int i = nxsecs; i < 2*nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).plate.sect(%d).zCamber = [", (i-nxsecs+1));
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
		{
			fprintf( file_id, "%f, ", degenPlate.zcamber[i][j] );
		}
		fprintf(	file_id, "%f];\n", degenPlate.zcamber[i][(num_pnts+1)/2-1]);
	}

	for ( int i = nxsecs; i < 2*nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).plate.sect(%d).t = [", (i-nxsecs+1));
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
		{
			fprintf( file_id, "%f, ", degenPlate.t[i][j] );
		}
		fprintf(	file_id, "%f];\n", degenPlate.t[i][(num_pnts+1)/2-1]);
	}

	for ( int i = nxsecs; i < 2*nxsecs; i++ )
	{
		fprintf(file_id, "\ndegenGeom(end).plate.sect(%d).nCamber = [", (i-nxsecs+1));
		for ( int j = 0; j < (num_pnts+1)/2-1; j++ )
		{
			fprintf(	file_id, "%f, %f, %f;\n",						\
					degenPlate.nCamber[i][j].x(),	\
					degenPlate.nCamber[i][j].y(),	\
					degenPlate.nCamber[i][j].z()	);
		}
		fprintf(	file_id, "%f, %f, %f];\n",										\
				degenPlate.nCamber[i][(num_pnts+1)/2-1].x(),	\
				degenPlate.nCamber[i][(num_pnts+1)/2-1].y(),	\
				degenPlate.nCamber[i][(num_pnts+1)/2-1].z()		);
	}

	fprintf(file_id, "\ndegenGeom(end).plate.u = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenPlate.u[i]);
	}
	fprintf(file_id, "%f];", degenPlate.u[2*nxsecs-1]);

	int wCnt = (num_pnts+1)/2;
	if (nxsecs != nxsecsOrig ) wCnt *= 2;
	fprintf(file_id, "\ndegenGeom(end).plate.wTop = [");
	for ( int j = wCnt; j < 2*wCnt-1; j++ )
	{
		fprintf(file_id, "%f, ", degenPlate.wTop[j]);
	}
	fprintf(file_id, "%f];", degenPlate.wTop[2*wCnt-1]);

	fprintf(file_id, "\ndegenGeom(end).plate.wBot = [");
	for ( int j = wCnt; j < 2*wCnt-1; j++ )
	{
		fprintf(file_id, "%f, ", degenPlate.wBot[j]);
	}
	fprintf(file_id, "%f];", degenPlate.wBot[2*wCnt-1]);

	//============================= DegenStick =============================//
	fprintf(file_id, "\ndegenGeom(end).stick.Xle = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenStick.xle[i].x(),	\
				degenStick.xle[i].y(),	\
				degenStick.xle[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenStick.xle[2*nxsecs-1].x(),	\
			degenStick.xle[2*nxsecs-1].y(),	\
			degenStick.xle[2*nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick.Xte = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenStick.xte[i].x(),	\
				degenStick.xte[i].y(),	\
				degenStick.xte[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenStick.xte[2*nxsecs-1].x(),	\
			degenStick.xte[2*nxsecs-1].y(),	\
			degenStick.xte[2*nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick.XcgSolid = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenStick.xcgSolid[i].x(),	\
				degenStick.xcgSolid[i].y(),	\
				degenStick.xcgSolid[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenStick.xcgSolid[2*nxsecs-1].x(),	\
			degenStick.xcgSolid[2*nxsecs-1].y(),	\
			degenStick.xcgSolid[2*nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick.XcgShell = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(	file_id,"%f, %f, %f;\n", degenStick.xcgShell[i].x(),	\
				degenStick.xcgShell[i].y(),	\
				degenStick.xcgShell[i].z()	);
	}
	fprintf(file_id,"%f, %f, %f];", degenStick.xcgShell[2*nxsecs-1].x(),	\
			degenStick.xcgShell[2*nxsecs-1].y(),	\
			degenStick.xcgShell[2*nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick.toc = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.toc[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.toc[2*nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.tLoc = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.tLoc[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.tLoc[2*nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.chord = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.chord[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.chord[2*nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.sweep = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.sweep[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.sweep[2*nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.Ishell = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
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
			degenStick.Ishell[2*nxsecs-1][0],	\
			degenStick.Ishell[2*nxsecs-1][1],	\
			degenStick.Ishell[2*nxsecs-1][2],	\
			degenStick.Ishell[2*nxsecs-1][3],	\
			degenStick.Ishell[2*nxsecs-1][4],	\
			degenStick.Ishell[2*nxsecs-1][5]	);

	fprintf(file_id, "\ndegenGeom(end).stick.Isolid = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(	file_id, "%f, %f, %f;\n",					\
				degenStick.Isolid[i][0],	\
				degenStick.Isolid[i][1],	\
				degenStick.Isolid[i][2]		);
	}
	fprintf(	file_id, "%f, %f, %f];\n",						\
			degenStick.Isolid[2*nxsecs-1][0],	\
			degenStick.Isolid[2*nxsecs-1][1],	\
			degenStick.Isolid[2*nxsecs-1][2]	);

	fprintf(file_id, "\ndegenGeom(end).stick.area = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.area[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.area[2*nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.areaNormal = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(	file_id, "%f, %f, %f;\n",						\
				degenStick.areaNormal[i].x(),	\
				degenStick.areaNormal[i].y(),	\
				degenStick.areaNormal[i].z()	);
	}
	fprintf(	file_id, "%f, %f, %f];\n",							\
			degenStick.areaNormal[2*nxsecs-1].x(),\
			degenStick.areaNormal[2*nxsecs-1].y(),\
			degenStick.areaNormal[2*nxsecs-1].z()	);

	fprintf(file_id, "\ndegenGeom(end).stick.perimTop = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.perimTop[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.perimTop[2*nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.perimBot = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.perimBot[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.perimBot[2*nxsecs-1]);

	fprintf(file_id, "\ndegenGeom(end).stick.u = [");
	for ( int i = nxsecs; i < 2*nxsecs-1; i++ )
	{
		fprintf(file_id, "%f, ", degenStick.u[i]);
	}
	fprintf(file_id, "%f];\n", degenStick.u[2*nxsecs-1]);

	nxsecs = nxsecsOrig;
	//============================= DegenPoint =============================//
	fprintf(file_id, "\ndegenGeom(end).point.vol = %f;\n", degenPoint.vol[1]);
	fprintf(file_id, "\ndegenGeom(end).point.volWet = %f;\n", degenPoint.volWet[1]);
	fprintf(file_id, "\ndegenGeom(end).point.area = %f;\n", degenPoint.area[1]);
	fprintf(file_id, "\ndegenGeom(end).point.areaWet = %f;\n", degenPoint.areaWet[1]);
	fprintf(file_id, "\ndegenGeom(end).point.Ishell = [%f, %f, %f, %f, %f, %f];\n",	\
			degenPoint.Ishell[1][0],	\
			degenPoint.Ishell[1][1],	\
			degenPoint.Ishell[1][2],	\
			degenPoint.Ishell[1][3],	\
			degenPoint.Ishell[1][4],	\
			degenPoint.Ishell[1][5]		);
	fprintf(file_id, "\ndegenGeom(end).point.Isolid = [%f, %f, %f, %f, %f, %f];\n",	\
			degenPoint.Isolid[1][0],	\
			degenPoint.Isolid[1][1],	\
			degenPoint.Isolid[1][2],	\
			degenPoint.Isolid[1][3],	\
			degenPoint.Isolid[1][4],	\
			degenPoint.Isolid[1][5]		);
	fprintf(file_id, "\ndegenGeom(end).point.xcgShell = [%f, %f, %f];\n",	\
			degenPoint.xcgShell[1].x(),	\
			degenPoint.xcgShell[1].y(),	\
			degenPoint.xcgShell[1].z()	);
	fprintf(file_id, "\ndegenGeom(end).point.xcgSolid = [%f, %f, %f];\n",	\
			degenPoint.xcgSolid[1].x(),	\
			degenPoint.xcgSolid[1].y(),	\
			degenPoint.xcgSolid[1].z()	);
}
