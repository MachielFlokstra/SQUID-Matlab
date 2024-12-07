function [yfit,J,par] = mpmsfun2(z,rot,reps,step,~,par,nfo,jac)
% INFO
%  Evaluates the fitfunction for a series of linked longitudinal and transverse mpms data (i.e. when data was collected
%  using both the pickup coils for each scan). The data needs to be ordered starting with all the longitudinal scans and
%  followed by their accompanying transverse scans (i.e. when using mpmsload(fname) make sure that the first N files in
%  fname are the longitudinal scans and the next N files their corresponding transverse scans). 
%
% INPUT
%  z      : row vector with z data of all scans in series
%  rot    : column vector with rotation per "new" scan // rot values are always in range [0 to 360]
%  reps   : column vector with number of scans per "new" scan
%  step   : column vector with number of steps per "new" scan
%  ~      : input not used (its where the type info goes when making the call from LMfit_mpms.m)
%  par    : row vector with fitparameter values of all scans in series
%  nfo    : 4-row matrix with for each fitparameter: fit flag / common flag / lower limit / upper limit
%  jac    : flag to calculate the Jacobian
%
% OUTPUT
%  yfit   : value of the fitfunction for the given set of fitparameter values
%  J      : Jaccobian matrix
%  par    : same as par, but with the common fitparameter values set correctly
%
% FITPARAMETERS
%  a1 = r (m)
%  a2 = phi0 (deg) // give phi0 freedom in range [0 to 360], then phi=rot-phi0 stays in range [-360 to +360]
%  a3 = z0 (m)
%  a4 = Mr (Am^2)
%  a5 = Mp (Am^2)
%  a6 = Mz (Am^2)
%  a7 = long c0 (V)     // background offset for longitudinal scan
%  a8 = long c1 (V/m)   // background slope for longitudinal scan
%  a9 = trans c0 (V)    // background offset for transverse scan
%  a10 = trans c1 (V/m) // background slope for transverse scan
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

nfo(1:4,9:10) = nfo(1:4,7:8); % add nfo for a9 and a10 (copies of the settings for a7 and a8)
scans = length(reps)/2; % the number of single logitudinal scans (= same for transverse)
zoff = length(z)/2; % offset for z to mark the start of the trans part

% allocate memory
%a = zeros(1,size(nfo,2)); % these will become a(1), a(2), ... a(10) for the scan being evaluated
yfit = zeros(1,length(z));
yfith = zeros(1,length(z));
J = zeros(length(z),length(par));

% correct par for common and linked fitparameters
osa = 0; % offset counter for par to mark start of current scan
for k=1:scans
  for j=1:reps(k)
    % check for common fitparameters
    for ai=1:size(nfo,2)
      if isequal(nfo(2,ai),1)
        % common parameter: take value from first scan of long block for a1-a8, and of trans block for a9-a10
        if ai < 9
          par(osa+ai) = par(ai); % take from first scan in long block
        else
          par(8*scans+osa+ai-2) = par(8*scans+ai-2); % take from first scan in trans block
        end
      end
    end
    % set linked parameters (trans block a1-a6 should be set to long block a1-a6)
    par(8*scans+osa+(1:6)) = par(osa+(1:6));
    osa = osa + 8; % set offset to next scan - only add 8 since par is structured using 8 fitparameters per scan
  end
end

% build the total signal, scan by scan
osa = 0; % offset counter for par to mark start of current scan
for k=1:scans
  for j=1:reps(k)
    tL = sum(reps(1:k-1).*step(1:k-1)) + (j-1)*step(k) + (1:step(k)); % index range for long part of current scan
    tT = zoff + tL; % index range for trans part of current scan
    a = [par(osa+(1:8)) par(8*scans+osa+(7:8))]; % fitparameters for current scan, take a9,a10 from trans part
    osa = osa + 8; % set offset to next scan
    % calculate the fitfunction for the long part
    F1 = -LC_RP(a(4),a(6),CR,z(tL)-a(3)-LC_pos(1),a(1),Nseg);
    F2 = +LC_RP(a(4),a(6),CR,z(tL)-a(3)-LC_pos(2),a(1),Nseg);
    F3 = +LC_RP(a(4),a(6),CR,z(tL)-a(3)-LC_pos(3),a(1),Nseg);
    F4 = -LC_RP(a(4),a(6),CR,z(tL)-a(3)-LC_pos(4),a(1),Nseg);
    yfit(1,tL) = VLx*((mu0*CR)/(4*pi))*(F1+F2+F3+F4) + a(7) + a(8)*z(tL);
    % calculate the fitfunction for the trans part
    F1R = +TC_RP(a(4),a(5),a(6),CR,1*Ch,phi1R,phi2R,z(tT)-a(3)-TC_pos(1),a(1),(rot(k)-a(2))*pi/180,Nseg);
    F1L = -TC_RP(a(4),a(5),a(6),CR,1*Ch,phi1L,phi2L,z(tT)-a(3)-TC_pos(1),a(1),(rot(k)-a(2))*pi/180,Nseg);
    F2R = -TC_RP(a(4),a(5),a(6),CR,2*Ch,phi1R,phi2R,z(tT)-a(3)-TC_pos(2),a(1),(rot(k)-a(2))*pi/180,Nseg);
    F2L = +TC_RP(a(4),a(5),a(6),CR,2*Ch,phi1L,phi2L,z(tT)-a(3)-TC_pos(2),a(1),(rot(k)-a(2))*pi/180,Nseg);
    F3R = +TC_RP(a(4),a(5),a(6),CR,1*Ch,phi1R,phi2R,z(tT)-a(3)-TC_pos(3),a(1),(rot(k)-a(2))*pi/180,Nseg);
    F3L = -TC_RP(a(4),a(5),a(6),CR,1*Ch,phi1L,phi2L,z(tT)-a(3)-TC_pos(3),a(1),(rot(k)-a(2))*pi/180,Nseg);
    yfit(1,tT) = VTx*((mu0*CR)/(4*pi))*(F1R+F1L+F2R+F2L+F3R+F3L) + a(9) + a(10)*z(tT);
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
  for k=1:scans
    for j=1:reps(k)
      tL = sum(reps(1:k-1).*step(1:k-1)) + (j-1)*step(k) + (1:step(k)); % index range for long part of current scan
      tT = zoff + tL; % index range for trans part of current scan
      % check for special conditions (first reset all flags)
      LMr = 0; % flag to set d/dMr to zero for longitudinal coil
      TMr = 0; % flag to set d/dMr to zero for transverse coil
      TMp = 0; % flag to set d/dMp to zero for transverse coil
      TMz = 0; % flag to set d/dMz to zero for transverse coil
      if isequal(par(osa+1),0) % check r=0
        LMr = 1;
        TMz = 1;
      end
      if ~isequal(sum(eq(rot(k)-par(osa+2),[-360 -180 0 180 360])),0) % check phi=[0,180]
        TMp = 1;
      end
      if ~isequal(sum(eq(rot(k)-par(osa+2),[-270 -90 90 270 450])),0) % check phi=[90,270]
        TMz = 1;
        TMr = 1;
      end
      % loop over the fitparameters
      for m=1:size(nfo,2)
        % only consider free fitparameters (for fixed fitparameters dyda=0 and nothing needs to be calculated)
        if isequal(nfo(1,m),1)
          % first make a copy of the original fitparameters (as used to calculate yfit)
          par_store = par;
          % now determine, and add, the small step to add to fitparameter a(m) of the current scan - remember that a9
          % and a10 are the a7 and a8 from the trans part
          if m < 9
            if isequal(par(osa+m),0)
              % add a fixed number rather then taking a fraction of the parameter value where care needs to be taken
              % that this value is sensible (i.e. that it is indeed a small value for the fitparameter)
              hstep = Omin * INCR;
            else
              hstep = par(osa+m)*INCR;
            end
            par(osa+m) = par(osa+m) + hstep;
          else
            if isequal(par(8*scans+osa+m-2),0)
              hstep = Omin * INCR;
            else
              hstep = par(8*scans+osa+m-2)*INCR;
            end
            par(8*scans+osa+m-2) = par(8*scans+osa+m-2) + hstep;
          end
          
          % get the fitparameters for current scan (which now include the small change to parameter m)
          a = [par(osa+(1:8)) par(8*scans+osa+(7:8))];

          % calculate the fitfunction for the long part
          F1 = -LC_RP(a(4),a(6),CR,z(tL)-a(3)-LC_pos(1),a(1),Nseg);
          F2 = +LC_RP(a(4),a(6),CR,z(tL)-a(3)-LC_pos(2),a(1),Nseg);
          F3 = +LC_RP(a(4),a(6),CR,z(tL)-a(3)-LC_pos(3),a(1),Nseg);
          F4 = -LC_RP(a(4),a(6),CR,z(tL)-a(3)-LC_pos(4),a(1),Nseg);
          yfith(1,tL) = VLx*((mu0*CR)/(4*pi))*(F1+F2+F3+F4) + a(7) + a(8)*z(tL);
          % calculate the fitfunction for the trans part
          F1R = +TC_RP(a(4),a(5),a(6),CR,1*Ch,phi1R,phi2R,z(tT)-a(3)-TC_pos(1),a(1),(rot(k)-a(2))*pi/180,Nseg);
          F1L = -TC_RP(a(4),a(5),a(6),CR,1*Ch,phi1L,phi2L,z(tT)-a(3)-TC_pos(1),a(1),(rot(k)-a(2))*pi/180,Nseg);
          F2R = -TC_RP(a(4),a(5),a(6),CR,2*Ch,phi1R,phi2R,z(tT)-a(3)-TC_pos(2),a(1),(rot(k)-a(2))*pi/180,Nseg);
          F2L = +TC_RP(a(4),a(5),a(6),CR,2*Ch,phi1L,phi2L,z(tT)-a(3)-TC_pos(2),a(1),(rot(k)-a(2))*pi/180,Nseg);
          F3R = +TC_RP(a(4),a(5),a(6),CR,1*Ch,phi1R,phi2R,z(tT)-a(3)-TC_pos(3),a(1),(rot(k)-a(2))*pi/180,Nseg);
          F3L = -TC_RP(a(4),a(5),a(6),CR,1*Ch,phi1L,phi2L,z(tT)-a(3)-TC_pos(3),a(1),(rot(k)-a(2))*pi/180,Nseg);
          yfith(1,tT) = VTx*((mu0*CR)/(4*pi))*(F1R+F1L+F2R+F2L+F3R+F3L) + a(9) + a(10)*z(tT);
          
          % check for special cases (first reset flag)
          JzeroL = 0; % flag to set derivatives dy/da(m) to zero for long part
          JzeroT = 0; % flag to set derivatives dy/da(m) to zero for trans part
          if isequal(LMr,1) && isequal(m,4) % LMr flagged & m=4(Mr)
            JzeroL = 1;
          end
          if isequal(TMr,1) && isequal(m,4) % TMr flagged & m=4(Mr)
            JzeroT = 1;
          end
          if isequal(TMp,1) && isequal(m,5) % TMp flagged & m=5(Mp)
            JzeroT = 1;
          end
          if isequal(TMz,1) && isequal(m,6) % TMz flagged & m=6(Mz)
            JzeroT = 1;
          end
          
          % determine dyda for long part
          if isequal(JzeroL,0)
            dydaL = (yfith(1,tL) - yfit(1,tL)) / hstep;
          else
            dydaL = 0; % don't need to bother setting the correct length of vector here, can leave it single valued
          end
          % determine dyda for trans part
          if isequal(JzeroT,0)
            dydaT = (yfith(1,tT) - yfit(1,tT)) / hstep;
          else
            dydaT = 0; % don't need to bother setting the correct length of vector here, can leave it single valued
          end
          
          % Store in the Jaccobian under the correct scan. Derivatives to the fitparameters a1 to a8 are always placed
          % in the long block those to a9 and a10 are place in the trans block (unless common, the in first block).
          if isequal(nfo(2,m),1) % common: store under fitparameters of first scan
            J(tL,m) = dydaL;
            J(tT,m) = dydaT;
          elseif m > 8 % a9(a10): store under trans block a7(a8) if not common
            %J(tL,8*scans+osa+m-2) = dydaL; % dydaL for a9 and a10 will always be zero, so no need to set it
            J(tT,8*scans+osa+m-2) = dydaT; % dydaL for a9 and a10 will always be zero
          else % unique (a1 to a8): store under fitparameters of current scan
            J(tL,osa+m) = dydaL;
            J(tT,osa+m) = dydaT;
          end
          
          % restore the original parameter values
          par = par_store;
        end
      end % this closes the loop over m (the fitparameters)
      osa = osa + 8; % set offset to next scan
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