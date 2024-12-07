function [yfit,J,par] = mpmsfun(z,rot,reps,step,type,par,nfo,jac)
% INFO
%  Evaluates the fitfunction for a series of mpms data.
%
% INPUT
%  z      : row vector with z data of all scans in series
%  rot    : column vector with rotation per "new" scan // rot values are always in range [0 to 360]
%  reps   : column vector with number of scans per "new" scan
%  step   : column vector with number of steps per "new" scan
%  type   : column vector with the coil-system used per "new" scan (0=longitudinal, 1=transverse)
%  par    : row vector with fitparameter values of all scans in series
%  nfo    : 4-row matrix with for each fitparameter: fit flag / common flag / lower limit / upper limit
%  jac    : flag to calculate the Jacobian
%
% OUTPUT
%  yfit   : value of the fitfunction for the given set of fitparameter values
%  J      : Jaccobian matrix
%  par    : same as (input) par, but with the common fitparameter values set correctly
%
% FITPARAMETERS
%  a1 = r (m)
%  a2 = phi0 (deg) // give phi0 freedom in range [0 to 360], then phi=rot-phi0 stays in range [-360 to +360]
%  a3 = z0 (m)
%  a4 = Mr (Am^2)
%  a5 = Mp (Am^2)
%  a6 = Mz (Am^2)
%  a7 = c0 (V) // the background offset
%  a8 = c1 (V/m) // the background slope
%
% SPECIAL DIPOLE POSITIONS
%  There are several particular values for r and phi where some of the component of M the net contribution to the
%  measurement signal should be zero (they may have non-zero contribution from each individual coil, but the sum should
%  be sero). However, due to numerical accuracy they are not and as function of z appear to generate a relative smooth
%  curve which the fitting routine will amplify by setting extremely large values for M such that the smooth curve
%  becomes comparable to the measurement signal. For these cases, the derivatives should be set to zero (one could also
%  set the net contribution of the particular M component to the fitfunction to zero). Note that for the longitudinal
%  coils, the fitfunction does not depend on phi and Mp (and there won't be an accuracy issue between yfit and yfith).
%
%  LONGITUDINAL COILS
%   r=0           : Mr contribution cancels out         => set d/dMr components to zero
%
%  TRANSVERSE COILS
%   r=0           : Mz contribition cancels out         => set d/dMz components to zero
%   phi=[0,180]   : Mp contribution cancels out         => set d/dMp components to zero
%   phi=[90,270]  : Mz and Mr contributions cancel out  => set d/dMz and d/dMr components to zero
% ----------------------------------------------------------------------------------------------------------------------

% constants
mu0 = 4*pi*1E-7; % Vs/(Am)
Nseg = 15; % number of angle-segments used for numerical integration
INCR = 1E-3; % the small step for determining numerical derivatives
Omin = 1E-10; % the smallest order of magnitude for any non-zero fitparameter (detection limit for M ~ 1E-10 Am^2)
VLx = 1.633E7; % conversion factor (Hz) from magnetic flux (Vs) into scaled voltage (V) for longitudinal pickup coils
VTx = 1.46*2.828E8; % conversion factor (Hz) from magnetic flux (Vs) into scaled voltage (V) for transverse pickup coils
 % note1: using VLx = sqrt(8/3)*1E7 and VTx = sqrt(8)*1E8 gives near identical results compared to the mpms
 % note2: the two pickup coil systems may not yield the same result when comparing a moment aligned along z (at r=0) and
 %        making a longitudinal scan, compared to the same moment aligned along r (phi=0) and making a transverse scan. 
 % note3: the prefactor 1.46 is for now a working factor to have mathcing results for long vs trans
 
% transverse coil system (the longitudinal values are known from Quantum Desisgn manuals/documentation)
CR = 1.01E-2; % coil radius (m) // the added 0.1mm is likely the thickness of the insulation of the wires
Ccir = 185; % coil "circumference" (deg)
Cgap = 90; % angle around which the gap between the paired coils is situated (second gap at +180)
Ch = 1E-2; % bot/top coil height (m) // the center coil pair has double the height
Cs = 0.2E-3; % vertical spacing between coils (m) // 2x the 0.1mm insulation is the minimum distance between the coils

% determine coil centers (bot to top)
TC_pos = [-(1.5*Ch+Cs) 0 (1.5*Ch+Cs)]; % transverse coil system
LC_pos = [-1.52 -0.05 0.05 1.52]*1E-2; % longitudinal coil system

% determine outer angles of the transverse coils (R coil "centered" around phi=0, L coil around phi=180)
phi1R = (pi/180)*(Cgap-90-0.5*Ccir); % lower angle limit for R coil
phi2R = (pi/180)*(Cgap-90+0.5*Ccir); % upper angle limit for R coil
phi1L = (pi/180)*(Cgap+90-0.5*Ccir); % lower angle limit for L coil
phi2L = (pi/180)*(Cgap+90+0.5*Ccir); % upper angle limit for L coil

% allocate memory
a = zeros(1,size(nfo,2)); % these will become a(1), a(2), ... a(8) for the scan being evaluated
yfit = zeros(1,length(z));
yfith = zeros(1,length(z));
J = zeros(length(z),length(par));

% correct par of each scan according to the common fitparameters
osa = 0; % offset counter for par to mark start of current scan
for k=1:length(reps)
  for j=1:reps(k)
    for ai=1:size(nfo,2)
      if isequal(nfo(2,ai),1)
        par(osa+ai) = par(ai); % common parameter: take value from first run
      end
    end
    osa = osa + size(nfo,2); % set offset to next scan
  end
end

% build the total signal, scan by scan
osa = 0; % offset counter for par to mark start of current scan
for k=1:length(reps)
  for j=1:reps(k)
    t = sum(reps(1:k-1).*step(1:k-1)) + (j-1)*step(k) + (1:step(k)); % index range for current scan
    a(1:size(nfo,2)) = par(osa+(1:size(nfo,2))); % fitparameters for current scan
    osa = osa + size(nfo,2); % set offset to next scan
    % calculate the fitfunction
    if isequal(type(k),0)
      F1 = -LC_RP(a(4),a(6),CR,z(t)-a(3)-LC_pos(1),a(1),Nseg);
      F2 = +LC_RP(a(4),a(6),CR,z(t)-a(3)-LC_pos(2),a(1),Nseg);
      F3 = +LC_RP(a(4),a(6),CR,z(t)-a(3)-LC_pos(3),a(1),Nseg);
      F4 = -LC_RP(a(4),a(6),CR,z(t)-a(3)-LC_pos(4),a(1),Nseg);
      yfit(1,t) = VLx*((mu0*CR)/(4*pi))*(F1+F2+F3+F4) + a(7) + a(8)*z(t);
    else
      F1R = +TC_RP(a(4),a(5),a(6),CR,1*Ch,phi1R,phi2R,z(t)-a(3)-TC_pos(1),a(1),(rot(k)-a(2))*pi/180,Nseg);
      F1L = -TC_RP(a(4),a(5),a(6),CR,1*Ch,phi1L,phi2L,z(t)-a(3)-TC_pos(1),a(1),(rot(k)-a(2))*pi/180,Nseg);
      F2R = -TC_RP(a(4),a(5),a(6),CR,2*Ch,phi1R,phi2R,z(t)-a(3)-TC_pos(2),a(1),(rot(k)-a(2))*pi/180,Nseg);
      F2L = +TC_RP(a(4),a(5),a(6),CR,2*Ch,phi1L,phi2L,z(t)-a(3)-TC_pos(2),a(1),(rot(k)-a(2))*pi/180,Nseg);
      F3R = +TC_RP(a(4),a(5),a(6),CR,1*Ch,phi1R,phi2R,z(t)-a(3)-TC_pos(3),a(1),(rot(k)-a(2))*pi/180,Nseg);
      F3L = -TC_RP(a(4),a(5),a(6),CR,1*Ch,phi1L,phi2L,z(t)-a(3)-TC_pos(3),a(1),(rot(k)-a(2))*pi/180,Nseg);
      yfit(1,t) = VTx*((mu0*CR)/(4*pi))*(F1R+F1L+F2R+F2L+F3R+F3L) + a(7) + a(8)*z(t);
    end
  end
end

% calculate Jaccobian = dyda ~ (y(a+h)-y(a))/(a+h - a) = (y(a+h)-y(a))/h
if isequal(jac,1)
  % the only non-zero elements are the "blocks" inside J that belong to a particular scan with its fitparameters, with
  % the exception that yfith for common fitparameters are stored in the block belonging to the first scan. It is thus
  % sufficient to only calculate dyda for the fitparameters belonging to the scan under evaluation where for common
  % fitparameters the value of the fitparameter is taken from the first scan and the result is stored under the first
  % scan as well (which creates the link between the first and current scan).
  osa = 0; % offset counter for par to mark start of current scan
  for k=1:length(reps)
    for j=1:reps(k)
      t = sum(reps(1:k-1).*step(1:k-1)) + (j-1)*step(k) + (1:step(k)); % index range for current scan
      % check for special conditions (first reset all flags)
      LMr = 0; % flag to set d/dMr to zero for longitudinal coil
      TMr = 0; % flag to set d/dMr to zero for transverse coil
      TMp = 0; % flag to set d/dMp to zero for transverse coil
      TMz = 0; % flag to set d/dMz to zero for transverse coil
      if isequal(par(osa+1),0) && isequal(type(k),0) % check r=0 (long)
        LMr = 1;
      end
      if isequal(par(osa+1),0) && isequal(type(k),1) % check r=0 (trans)
        TMz = 1;
      end
      if ~isequal(sum(eq(rot(k)-par(osa+2),[-360 -180 0 180 360])),0) && isequal(type(k),1) % check phi=[0,180]
        TMp = 1;
      end
      if ~isequal(sum(eq(rot(k)-par(osa+2),[-270 -90 90 270])),0) && isequal(type(k),1) % check phi=[90,270]
        TMz = 1;
        TMr = 1;
      end
      % loop over the fitparameters (for a single scan)
      for m=1:size(nfo,2)
        % only consider free fitparameters (for fixed fitparameters dyda=0 and nothing needs to be calculated)
        if isequal(nfo(1,m),1)
          % first make a copy of the original fitparameters (as used to calculate yfit)
          par_store = par;
          % now determine, and add, the small step to add to fitparameter a(m) of the current scan
          if isequal(par(osa+m),0)
            % add a fixed number rather then taking a fraction of the parameter value where care needs to be taken that
            % this value is sensible (i.e. that it is indeed a small value for the fitparameter)
            hstep = Omin * INCR;
          else
            hstep = par(osa+m)*INCR;
          end
          par(osa+m) = par(osa+m) + hstep;
          % get the fitparameters for current scan (which now include the small change to parameter m)
          a(1:size(nfo,2)) = par(osa+(1:size(nfo,2)));
                    
          % calculate the fitfunction
          if isequal(type(k),0)
            F1 = -LC_RP(a(4),a(6),CR,z(t)-a(3)-LC_pos(1),a(1),Nseg);
            F2 = +LC_RP(a(4),a(6),CR,z(t)-a(3)-LC_pos(2),a(1),Nseg);
            F3 = +LC_RP(a(4),a(6),CR,z(t)-a(3)-LC_pos(3),a(1),Nseg);
            F4 = -LC_RP(a(4),a(6),CR,z(t)-a(3)-LC_pos(4),a(1),Nseg);
            yfith(1,t) = VLx*((mu0*CR)/(4*pi))*(F1+F2+F3+F4) + a(7) + a(8)*z(t);
          else
            F1R = +TC_RP(a(4),a(5),a(6),CR,1*Ch,phi1R,phi2R,z(t)-a(3)-TC_pos(1),a(1),(rot(k)-a(2))*pi/180,Nseg);
            F1L = -TC_RP(a(4),a(5),a(6),CR,1*Ch,phi1L,phi2L,z(t)-a(3)-TC_pos(1),a(1),(rot(k)-a(2))*pi/180,Nseg);
            F2R = -TC_RP(a(4),a(5),a(6),CR,2*Ch,phi1R,phi2R,z(t)-a(3)-TC_pos(2),a(1),(rot(k)-a(2))*pi/180,Nseg);
            F2L = +TC_RP(a(4),a(5),a(6),CR,2*Ch,phi1L,phi2L,z(t)-a(3)-TC_pos(2),a(1),(rot(k)-a(2))*pi/180,Nseg);
            F3R = +TC_RP(a(4),a(5),a(6),CR,1*Ch,phi1R,phi2R,z(t)-a(3)-TC_pos(3),a(1),(rot(k)-a(2))*pi/180,Nseg);
            F3L = -TC_RP(a(4),a(5),a(6),CR,1*Ch,phi1L,phi2L,z(t)-a(3)-TC_pos(3),a(1),(rot(k)-a(2))*pi/180,Nseg);
            yfith(1,t) = VTx*((mu0*CR)/(4*pi))*(F1R+F1L+F2R+F2L+F3R+F3L) + a(7) + a(8)*z(t);
          end
          
          % check for special cases (first reset flag)
          Jzero = 0; % flag to set derivatives dy/da(m) to zero
          if isequal(type(k),0) && isequal(LMr,1) && isequal(m,4) % long & LMr flagged & m=4(Mr)
            Jzero = 1;
          end
          if isequal(type(k),1) && isequal(TMr,1) && isequal(m,4) % trans & TMr flagged & m=4(Mr)
            Jzero = 1;
          end
          if isequal(type(k),1) && isequal(TMp,1) && isequal(m,5) % trans & TMp flagged & m=5(Mp)
            Jzero = 1;
          end
          if isequal(type(k),1) && isequal(TMz,1) && isequal(m,6) % trans & TMz flagged & m=6(Mz)
            Jzero = 1;
          end
          
          % determine dyda
          if isequal(Jzero,0)
            dyda = (yfith(1,t) - yfit(1,t)) / hstep;
          else
            dyda = 0; % don't need to bother setting the correct length of vector here, can leave it single valued
          end
          
          % store in the Jaccobian under the correct scan
          if isequal(nfo(2,m),1)
            J(t,m) = dyda; % common: store under fitparameters of first scan
          else
            J(t,osa+m) = dyda; % unique: store under fitparameters of current scan
          end
          
          % restore the original parameter values
          par = par_store;
        end
      end % this closes the loop over m (the fitparameters)
      osa = osa + size(nfo,2); % set offset to next scan
    end
  end
end
end

% ======================================================================================================================
function F = TC_RP(Mr,Mp,Mz,CR,L,phi1,phi2,z,r,phi,Nseg)
% INFO
%  F*(mu0*CR)/(4*pi) is the enclosed flux (in Tm^2 = Vs) of a single transverse/saddle coil, due to a dipole with moment
%  M=(Mr,Mp,Mz) and z-location z. So when z<0 the dipole is located below the coil central xy-plane and for z>0 above.
%
%  The path around a transverse coil (a cylinder cut in half) consist of two (horizontal) ring segments connected by two
%  (vertical) straight lines. The horizontal paths are p1 and p3, the vertical paths are p2 and p4. For p1 and p3 an
%  integration over the coil circumference needs to be made at each value of z.
%
%  The direction of the path is taken counter-clockwise when looking at the coil from inside the cylinder (i.e. when
%  viewing from the origin). The sign of the measured voltage depends on the current direction. It is taken positive for
%  a counter-clockwise current and negative for a clockwise current. At this stage it doesn't matter if it is defined in
%  this way or with opposite sign since there is a global sign for the total signal (combined voltages of all coils)
%  still to decide on.
%
% INPUT
%  Mr, Mp, Mz : components of the moment M along r, phi and z
%  CR         : coil radius
%  L          : coil height
%  phi1, phi2 : lower angle and upper angle of the coil
%  z          : vector with z-data (should have been corrected by z0 and coil position)
%  r          : radial distance from coil center axis (m)
%  phi        : angular coordinate (azimuth) (rad)
%  Nseg       : number of angle-segments used for numerical integration
%
% OUTPUT
%  F          : F*(mu0*CR)/(4*pi) is the enclosed flux (Tm^2 = Vs)
% ----------------------------------------------------------------------------------------------------------------------
dQ = abs((phi2-phi1)/Nseg); % small ring segment (Nseg of these segments precisely fill the integration range)
Q = phi1 + 0.5*dQ + (0:(Nseg-1))*dQ; % midpoints of all the segments
 % important to set it at the midpoint, not at the endpoints which leads to larger errors (i.e. requires smaller
 % segments to be as accurate)
p1 = zeros(1,length(z));
p3 = zeros(1,length(z));
for k=1:length(z)
  % path1
  denom = (r^2 + (-z(k)+L/2)^2 + CR^2 - 2*r*CR*cos(Q-phi)).^(3/2);
  p1(k) = sum(dQ*(-Mr*(-z(k)+L/2)*cos(Q-phi) - Mp*(-z(k)+L/2)*sin(Q-phi) + Mz*(CR-r*cos(Q-phi)))./denom);
  % path3
  denom = (r^2 + (-z(k)-L/2)^2 + CR^2 - 2*r*CR*cos(Q-phi)).^(3/2);
  p3(k) = -sum(dQ*(-Mr*(-z(k)-L/2)*cos(Q-phi) - Mp*(-z(k)-L/2)*sin(Q-phi) + Mz*(CR-r*cos(Q-phi)))./denom);
end
% path2
temp = r^2 + CR^2 - 2*r*CR*cos(phi2-phi);
INT = (-z+L/2)./(temp*((-z+L/2).^2+temp).^(1/2)) - (-z-L/2)./(temp*((-z-L/2).^2+temp).^(1/2));
p2 = -(Mr*sin(phi2-phi)-Mp*cos(phi2-phi)+Mp*r/CR) * INT;
% path4
temp = r^2 + CR^2 - 2*r*CR*cos(phi1-phi);
INT = (-z+L/2)./(temp*((-z+L/2).^2+temp).^(1/2)) - (-z-L/2)./(temp*((-z-L/2).^2+temp).^(1/2));
p4 = (Mr*sin(phi1-phi)-Mp*cos(phi1-phi)+Mp*r/CR) * INT;
% total
F = p1+p2+p3+p4;
end

% ======================================================================================================================
function F = LC_RP(Mr,Mz,CR,z,r,Nseg)
% INFO
%  F*(mu0*CR)/(4*pi) is the enclosed flux (in Tm^2 = Vs) of a single longitudinal coil, due to a dipole with moment
%  M=(Mx,My,Mz) and z-location z. So when z<0 the dipole is located below the coil central xy-plane and for z>0 above.
%
%  The path around a longitudinal coil is a circle. An integration over Q needs to be made at each value of z (similar
%  to the ring segments of the transverse coils).
%
% INPUT
%  Mr, Mz : components of the moment M along r, phi and z
%  CR     : coil radius
%  z      : vector with z-data (should have been corrected by z0 and coil position)
%  r      : radial distance from coil center axis (m)
%  Nseg   : number of angle-segments used for numerical integration
%
% OUTPUT
%  F      : F*(mu0*CR)/(4*pi) is the enclosed flux (Tm^2 = Vs)
% ----------------------------------------------------------------------------------------------------------------------
dQ = 2*pi/Nseg;
Q = 0.5*dQ + (0:(Nseg-1))*dQ;
F = zeros(1,length(z));
for k=1:length(z)
  denom = (r^2 + z(k)^2 + CR^2 - 2*r*CR*cos(Q)).^(3/2);
  F(k) = sum(dQ*(Mr*z(k)*cos(Q) + Mz*(CR-r*cos(Q)))./denom); % the Mp part gives zero contribution, and phi no effect
  % can redefine Q as Q - phi, which leaves dQ as dQ and integration is effectively still over the 0 to 2pi range
end
end