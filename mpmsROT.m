function [z0,Mz,ROTr,ROTp,ROTmr,ROTmp] = mpmsROT(fname,fdir,r,Mz)
% INFO
%  Fits a set of rotation data (longitudinal + transverse) to determine the actual rotation path the sample is making.
%
% INPUT
%  fname : row cell with the filenames of the longitudinal and transverse rotation scan
%  fdir  : directory containing the datafiles
%  r     : initial value for r
%  Mz    : initial/fix value for Mz // set Mz=1i to estimate initial value from raw data, else Mz is fixed at set value
%
% OUTPUT
%  z0     : obtained z0
%  Mz     : obtained Mz (from the longitudinal data with common Mz and Mr=0)
%  ROTr   : 2-column matrix with [rot r]
%  ROTp   : 2-column matrix with [rot phi0]
%  ROTmr  : 2-column matrix with [rot Mr]
%  ROTmp  : 2-column matrix with [rot Mp]
% ----------------------------------------------------------------------------------------------------------------------

% get the filenames
fileLO = fname{1};
fileTR = fname{2};

% create Pnfo and set standard fitparameter limits
Pnfo = zeros(4,8); % 4x8 matrix containing all required fitparameter infos for r,phi0,z0 Mr,Mphi,Mz,c0,c1
Pnfo(3:4,:) = mpmsfitlim();

% first make a fit to find z0 (in case centering was far off), to prevent r from jumping to limits
FitAll = 1;
Pnfo(1,:) = [1 0 1 0 0 1 1 1]; % fit z0,Mz, with z0,Mz common
Pnfo(2,:) = [0 0 1 0 0 1 0 0];
Pinit = [2E-3 0 1i 0 0 1i 0 0];
squidfit(fileLO,fdir,'mpmsfun',Pinit,Pnfo,FitAll,[],[],[],[]);
close(gcf); pause(0.001);
load('FitResults','aA');
z0 = aA(1,3); %#ok<*NODEF>

% fit the longitudinal data to find r(rot),z0,Mz // set Mr=0 (required Mz > Mr)
if isequal(Mz,1i)
  Pnfo(1,:) = [1 0 1 0 0 1 1 1]; % fit r,z0,Mz, with z0,Mz common
  Pnfo(2,:) = [0 0 1 0 0 1 0 0];
else
  Pnfo(1,:) = [1 0 1 0 0 0 1 1]; % fit r,z0 with z0 common (and Mz fixed at input value)
  Pnfo(2,:) = [0 0 1 0 0 0 0 0];
end
Pinit = [r 0 z0 0 0 Mz 0 0];
squidfit(fileLO,fdir,'mpmsfun',Pinit,Pnfo,FitAll,[],[],[],[]);
close(gcf); pause(0.001);
load('FitResults','Wqnty','aA');
z0 = aA(1,3);
Mz = aA(1,6);
ROTr = [Wqnty(:,4) aA(:,1)];

% fit the transverse data to find phi0(rot) // set Mr=Mphi=0 (requires Mz >> Mr,Mphi)
%  > If the dipole would rotate along a perfect circle, centered on the z-axis, then phi0 would be the same value for
%    all scans.
FitAll = 0;
Pnfo(1,:) = [0 1 0 0 0 0 1 1]; % fit phi0 per scan
Pnfo(2,:) = [0 0 0 0 0 0 0 0];
Pinit = [0 0 z0 0 0 Mz 0 0]; % fix the obtained z0,Mz
squidfit(fileTR,fdir,'mpmsfun',Pinit,Pnfo,FitAll,ROTr,[],[],[]); % apply r(rot)
close(gcf); pause(0.001);
load('FitResults','Wqnty','aA');

ROTp = [Wqnty(:,4) aA(:,2)]; % the obtained phi0 for each scan
phi = ROTp(:,1) - ROTp(:,2); % the obtained phi = rot - phi0, which will have values in the range [-360 to +360]
% take out any integers of 360 from phi = rot - phi0
phi = phi + 360; % now phi has values in the range [0 720] so we can subtract higher orders
phi = rem(phi,360); % subtract higher orders of 360

% now try to correct degenerate solutions by reflection in phi=0 axis and to keep phi(rot) smooth
%  > If Mr=Mphi=0, the voltage response is symmetric around phi=0 (i.e. phi=20 and phi=340 give same result). 
%  > Scan the obtained phi(rot) for large jumps (which is likely due to the degeneracy) and apply symmetry to correct.

% reflect phi for the first scan into the "positive" side (the [0-180] part)
if phi(1) > 180
  phi(1) = 360-phi(1);
end
% keep all subsequent values "smooth" (i.e. check which of the two solutions, phi(k) and 360-phi(k), is closest to
% phi(k-1) and keep that one as the correct solution).
for k=2:length(phi)
  if abs(phi(k)-phi(k-1)) > abs(360-phi(k)-phi(k-1))  
    phi(k) = 360-phi(k);
  end
end
% apply the changes to phi0(rot) = rot - phi(rot)
ROTp(:,2) = ROTp(:,1) - phi;

% ----------------------------------------------------------------------------------------------------------------------

% now fit data fixing the obtained z0,r(rot),phi0(rot) to obtain Mr(rot),Mp(rot) (and check that Mz~constant)
%  > Before they were kept at zero and as long as they are small compared to Mz it should not have made a change in the
%    results for z0,r(rot),phi0(rot) since Mz has opposite symmetry compared to Mr,Mphi. For long Mz is even, the others
%    are odd, and for trans Mz is odd and the others are even.

% fit longitudinal data to obtain Mr(rot)
FitAll = 0;
Pnfo(1,:) = [0 0 0 1 0 1 1 1]; % fit Mr,Mz per scan with Mz common
Pnfo(2,:) = [0 0 0 0 0 0 0 0];
Pinit = [0 0 z0 0 0 Mz 0 0];
squidfit(fileLO,fdir,'mpmsfun',Pinit,Pnfo,FitAll,ROTr,[],[],[]); % apply ROTr (ROTp no influence on long)
close(gcf); pause(0.001);
% load data to merge it with next
load('FitResults','fname','ZV','ZVi','Wqnty','Wreps','FNS','NSi','Vfit','Pchi','a','ZVA','ZVAi','VAfit','PAchi','aA');
ROTmr = [Wqnty(:,4) aA(:,4)];

% fit transverse data to obtain Mphi(rot) => does it makes sense? // does the total moment remains ~ constant?
FitAll = 0;
Pnfo(1,:) = [0 0 0 0 1 1 1 1]; % fit Mphi per scan
Pnfo(2,:) = [0 0 0 0 0 0 0 0];
Pinit = [0 0 z0 0 0 aA(1,6) 0 0]; % use for Mz the value obtained from previous long fits
squidfit(fileTR,fdir,'mpmsfun',Pinit,Pnfo,FitAll,ROTr,ROTp,ROTmr,[]); % apply ROTr,ROTp,ROTmr
close(gcf); pause(0.001);
mergedata(fname,ZV,ZVi,Wqnty,Wreps,FNS,NSi,Vfit,Pchi,a,ZVA,ZVAi,VAfit,PAchi,aA)

% ----------------------------------------------------------------------------------------------------------------------

% now plot the obtained results for the rotation path and the "transverse" magnetization vector
load('FitResults','fname','ZV','ZVi','Wqnty','Wreps','FNS','NSi','Vfit','Pchi','a','ZVA','ZVAi','VAfit','PAchi','aA');

% calculate transverse component (those normal to z axis) in cartesian // phi = ROTp(:,1) - ROTp(:,2)
dataf = 2; % need the fit results of the transverse data file to obtain ROTmp (ROTp and ROTmr are already known)
ri = sum(FNS(1:(dataf-1)))+(1:FNS(dataf));
ROTmp = [Wqnty(ri,4) aA(ri,5)];
Mx = ROTmr(:,2).*cos(phi*pi/180) - ROTmp(:,2).*sin(phi*pi/180);
My = ROTmr(:,2).*sin(phi*pi/180) + ROTmp(:,2).*cos(phi*pi/180);
[u,v] = pol2cart(phi*pi/180,ROTr(:,2)); % convert polar coordinates to cartesian

% create single figure
h = figure;
set(h,'NumberTitle','off','Name',['Rotation Data: ' fileLO],'Position',[100 100 1200 420]);
pan = uipanel(h,'position',[0 0 1/3 1]);
ha1 = polaraxes(pan); % axes to plot rotation path (polar)
pan = uipanel(h,'position',[1/3 0 1/3 1]);
ha2 = axes(pan); % axes to plot Mxy along path (cartesian)
pan = uipanel(h,'position',[2/3 0 1/3 1]);
ha3 = axes(pan); % axes to plot Mxy vs rotation

% plot the rotation path on polar axes
polarplot(ha1,phi*pi/180,ROTr(:,2),'-k');
ha1.Title.String = 'Rotation Path';
% plot the rotation path on cartesian axes to add transverse vectors
axes(ha2);
plot(ha2,1E3*u,1E3*v,'k'); % plot in mm
hold on; quiver(1E3*u,1E3*v,1E3*Mx,1E3*My,'b'); % add vectors (in emu)
ha2.Title.String = 'Mxy along rotation path';
ha2.XLabel.String = 'x (mm) [phi=0 line]';
ha2.YLabel.String = 'y (mm) [phi=90 line]';
% plot the length of the Mxy vector for each set rotation
quiver(ha3,Wqnty(1:FNS(1),4),zeros(length(Mx),1),zeros(length(Mx),1),1E3*sqrt(Mx.^2+My.^2),'b')
ha3.XLim = [min(Wqnty(1:FNS(1),4)) max(Wqnty(1:FNS(1),4))];
ha3.XLabel.String = 'set rotation (deg)';
ha3.YLabel.String = '|Mxy| (emu)';
end