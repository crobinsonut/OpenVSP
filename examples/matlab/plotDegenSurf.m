function h = plotDegenSurf( dgfile )

% Load VSP output degenerate geometry.
run(dgfile);

% Initialize flags to prevent undefined access.
plotnormals = false; plotuparm = false; plotwparm = false;

plotnormals = true;
%plotuparm = true;
%plotwparm = true;

ngeom = length(degenGeom);

figure(1)
clf
hold on

for i=1:ngeom
  disp(['Component ' num2str(i) ' Name: ' degenGeom(i).name]);
  
  [nxsecs, npt]=size(degenGeom(i).surf.x);

  x = degenGeom(i).surf.x;
  y = degenGeom(i).surf.y;
  z = degenGeom(i).surf.z;
  
  if(plotuparm)
    surf(x,y,z,repmat(degenGeom(i).surf.u,1,npt));
  elseif(plotwparm)
    surf(x,y,z,repmat(degenGeom(i).surf.w,1,nxsecs)');
  else
    surf(x,y,z);
  end

  if(plotnormals)
    % Find approximate center point
    cx = (x(1:end-1,1:end-1)+x(2:end,2:end))/2;
    cy = (y(1:end-1,1:end-1)+y(2:end,2:end))/2;
    cz = (z(1:end-1,1:end-1)+z(2:end,2:end))/2;
    
    quiver3(cx,cy,cz,degenGeom(i).surf.nx, degenGeom(i).surf.ny, degenGeom(i).surf.nz)
  end
end
hold off

axis equal
axis off







