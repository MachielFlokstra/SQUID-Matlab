function [Vfit,Pchi,a,ZVA,ZVAi,VAfit,PAchi,aA] = mpmsfit(ZV,Wqnty,Wreps,Wstep,Wtype,fitfun,Pinit,Pnfo,FitAll,ROTr,...
  ROTp,ROTmr,ROTmp,h)
% INFO
%  This function makes the calls to the Levenberg-Marquardt fitting algorithm and returns all the results. Fits are made
%  to all single scans, but also to the averages of any repeat scans.
%
% INPUT
%  ZV     : 3-column matrix with all measured data (Z,V,Verr)
%  Wqnty  : 4-column matrix with the quantity (H,<T>,<T>-error,Rotation) per "new" scan
%  Wreps  : column vector with the number of scans per "new" scan
%  Wstep  : column vector with the number of steps per "new" scan
%  Wtype  : column vector with the coil-system used per "new" scan (0=longitudinal, 1=transverse)
%  fitfun : string with name of the fitfunction to use
%  Pinit  : row vector with the initial values for r,phi0,z0,Mr,Mp,Mz,c0,c1
%   > for z0,Mr,Mp,Mz the value can be set to 1i to indicate to get the value estimated from the raw data
%  Pnfo   : 4x8 matrix containing all required fitparameter infos [fit flag; common flag; lower lim; upper lim]
%  FitAll : flag for simultaneous fitting // typically only use for a single rotation measurement file
%  ROTr   : 2-column matrix with [rot r] // if not empty, used to overwrite Pinit(1)
%  ROTp   : 2-column matrix with [rot phi0] // if not empty, used to overwrite Pinit(2)
%  ROTmr  : 2-column matrix with [rot Mr] // if not empty, used to overwrite Pinit(4)
%  ROTmp  : 2-column matrix with [rot Mp] // if not empty, used to overwrite Pinit(5)
%   > the ROT* variables should all be logically ordered with rot ranging from [0 to 360]
%  h      : handle of the edit box to print fit progression into // set h=0 to print to command window 
%
% OUTPUT
%  Vfit   : column vector for the obtained fits to each single scan
%  Pchi   : row vector with the chisquared of each fit belonging to Vfit
%  a      : matrix with single fit results "[a aerr]" per scan
%  ZVA    : 3-column matrix with the averaged data (avg version of ZV)
%  ZVAi   : 2-column matrix with the first and last index of ZVA belonging to each "new" scan
%  VAfit  : column vector for the obtained fits to each scan
%  PAchi  : row vector with the chisquared of each fit
%  aA     : matrix with the averages of the single fit results "[<a> <aerr> std(a)]" (first three blocks of 8 columns)
%           and the average repeat scan fit results (last two blocks of 8 columns)
%
% CALLS
%  LMfit_mpms.m (the Levenberg-Marquardt fitting algorithm)
% ----------------------------------------------------------------------------------------------------------------------

% constants
tol = 5E-3; % tolerance of the fit
Npar = 8; % number of fit parameters per measurement // not meant as a variable, but for clarity of reading the code

% allocate memory for the fit results for each individual scan
Vfit = zeros(size(ZV,1),1);
Pchi = zeros(1,sum(Wreps));
Pfit = zeros(1,Npar*sum(Wreps)); % obtained fitparameters (stored in series for easier dy/da calculation)
Perr = zeros(1,Npar*sum(Wreps)); % errors of the fitparameters
% allocate memory for repeat scan averages
ZVA = zeros(sum(Wstep),3);
ZVAi = zeros(length(Wstep),2);
VAfit = zeros(sum(Wstep),1);
PAchi = zeros(1,length(Wreps));
PAfit = zeros(1,Npar*length(Wreps)); % obtained fitparameters (stored in series for easier dy/da calculation)
PAerr = zeros(1,Npar*length(Wreps)); % errors of the fitparameters

% get the initial values from Pinit
r = Pinit(1); % radial distance from coil center axis (m)
phi0 = Pinit(2); % phi-shift (deg)
z0 = Pinit(3); % z-shift (m) // set to i to estimate from raw data (at half the scanlength)
Mr = Pinit(4); % magnetic moment along r (Am^2 = 1E3 emu) // set to i to estimate from raw data
Mp = Pinit(5); % magnetic moment along phi (Am^2 = 1E3 emu) // set to i to estimate from raw data
Mz = Pinit(6); % magnetic moment along z (Am^2 = 1E3 emu) // set to i to estimate from raw data
c0 = Pinit(7); % background offset (V) // by default set to zero
c1 = Pinit(8); % background slope (V/m) // by default set to zero

% determine first and last index of ZVA belonging to each "new" scan
j0 = 1; % index indicating the start of the current "new" scan
for k=1:length(Wreps)
  ZVAi(k,1:2) = [j0 j0+Wstep(k)-1];
  j0 = j0+Wstep(k); % set to next scan
end

% determine the repeat scan averages (also if none were used, for plotting purposes)
RepScans = 0; % flag to indicate repeat scans were used
if ~isequal(length(Wreps),sum(Wreps))
  RepScans = 1; % repeat scans were used
end
j0 = 0; % index offset indicating the start of the current scan
for k=1:length(Wreps) % loop over each "new" scan
  for j=1:Wreps(k) % loop over all reps for the current "new" scan
    ZVA(ZVAi(k,1):ZVAi(k,2),:) = ZVA(ZVAi(k,1):ZVAi(k,2),:)+ZV(j0+(1:Wstep(k)),:); % add scan
    j0 = j0 + Wstep(k); % set to next scan
  end
  ZVA(ZVAi(k,1):ZVAi(k,2),:) = ZVA(sum(Wstep(1:k-1))+(1:Wstep(k)),:) / Wreps(k); % get average
end

% initialize all fitparameters
j0 = 0; % index offset indicating the start of the current scan
for k=1:length(Wreps) % loop over each "new" scan
  % Apply ROT* variables if supplied. These need to be applied per scan to effectively set the correct Pinit value for
  % the rotation used for the current scan. The ROT* are used to make a linear interpolation to find the correct value
  % for * at the current rotation (Wqnty(k,4))
  if ~isempty(ROTr)
    r = interp1(ROTr(:,1),ROTr(:,2),Wqnty(k,4)); % set r according to ROTr
  end
  if ~isempty(ROTp)
    phi0 = interp1(ROTp(:,1),ROTp(:,2),Wqnty(k,4)); % set phi0 according to ROTp
  end
  if ~isempty(ROTmr)
    Mr = interp1(ROTmr(:,1),ROTmr(:,2),Wqnty(k,4)); % set Mr according to ROTmr
  end
  if ~isempty(ROTmp)
    Mp = interp1(ROTmp(:,1),ROTmp(:,2),Wqnty(k,4)); % set Mp according to ROTmp
  end
  % Estimates for Mr,Mp,Mphi are based on the maximum measured scaled voltage. A Mz moment of 1 emu (1E-3 Am^2)
  % produces a measured maximum scaled voltage of about 1.68 V on the longitudinal pickup coils (at r=0) while a Mr
  % moment of 0.1 emu produces a measured maximum scaled voltage of about 2.355 V on the transverse pickup coils (at
  % r=phi=0). So a reasonable initial guess is the maximum scaled voltage times (1/1.68)*1E-3 for longitudinal data,
  % and times 1/23.55*1E-3 for transverse data.
  if isequal(Wtype(k),0)
    Vx = (1/1.68)*1E-3;
  else
    Vx = (1/23.55)*1E-3;
  end
  for j=1:Wreps(k) % loop over all reps for the current "new" scan
    ios = Npar*(sum(Wreps(1:k-1))+j-1); % set the index offset in Pfit where parameters of current scan start
    % Set initial values for fitparameters (r,phi0,z0,Mr,Mp,Mz,c0,c1) of current scan.
    Pfit(ios+1) = r;
    Pfit(ios+2) = phi0;
    if isequal(z0,1i)
      Pfit(ios+3) = mean(ZV(j0+(1:Wstep(k)),1)); % set value at half the scanlength
    else
      Pfit(ios+3) = z0;
    end
    % Mr
    if isequal(Mr,1i)
      Pfit(ios+4) = Vx*max(ZV(j0+(1:Wstep(k)),2)); % set intial guess for Mr at the estimated maximum
    else
      Pfit(ios+4) = Mr;
    end
    % Mp
    if isequal(Mp,1i)
      Pfit(ios+5) = Vx*max(ZV(j0+(1:Wstep(k)),2)); % set intial guess for Mp at the estimated maximum
    else
      Pfit(ios+5) = Mp;
    end
    % Mz
    if isequal(Mz,1i)
      Pfit(ios+6) = Vx*max(ZV(j0+(1:Wstep(k)),2)); % set intial guess for Mz at the estimated maximum
    else
      Pfit(ios+6) = Mz;
    end
    Pfit(ios+7) = c0;
    Pfit(ios+8) = c1;
    j0 = j0 + Wstep(k); % set index to start of next scan
  end
  % similar as for single scans, but adapted to take the averages of the repeat scans
  if isequal(RepScans,1)
    PAfit(Npar*(k-1)+1) = r;
    PAfit(Npar*(k-1)+2) = phi0;
    if isequal(z0,1i)
      PAfit(Npar*(k-1)+3) = mean(ZVA(sum(Wstep(1:k-1))+(1:Wstep(k)),1));
    else
      PAfit(Npar*(k-1)+3) = z0;
    end
    % Mr
    if isequal(Mr,1i)
      PAfit(Npar*(k-1)+4) = Vx*max(ZVA(sum(Wstep(1:k-1))+(1:Wstep(k)),2));
    else
      PAfit(Npar*(k-1)+4) = Mr;
    end
    % Mp
    if isequal(Mp,1i)
      PAfit(Npar*(k-1)+5) = Vx*max(ZVA(sum(Wstep(1:k-1))+(1:Wstep(k)),2));
    else
      PAfit(Npar*(k-1)+5) = Mp;
    end
    % Mz
    if isequal(Mz,1i)
      PAfit(Npar*(k-1)+6) = Vx*max(ZVA(sum(Wstep(1:k-1))+(1:Wstep(k)),2));
    else
      PAfit(Npar*(k-1)+6) = Mz;
    end
    PAfit(Npar*(k-1)+7) = c0;
    PAfit(Npar*(k-1)+8) = c1;
  end
end

% make the fits: all scans simultaneously or all separately
if isequal(FitAll,1)
  fprintf('\nfitting all datafiles simultaneously');
  % check if only dealing with longitudinal data (then phi0 and Mphi should be kept fixed)
  if isequal(sum(Wtype),0)
    Pnfo(1,[2 5]) = 0; % flag phi0 and Mphi as fixed fitparameters
  end
  % call the fit algorithm
  [Vfit,Pchi(:),Pfit,Perr] = LMfit_mpms(fitfun,ZV(:,1)',ZV(:,2)',ZV(:,3)',Wqnty(:,4),Wreps,Wstep,Wtype,Pfit,Pnfo,tol);
  %  > chisquared will be single valued and the value will be attached to the Pchi for each scan
  % if repeat scans were used, also fit the average scans
  if isequal(RepScans,1)
    [VAfit,PAchi(:),PAfit,PAerr] = LMfit_mpms(fitfun,ZVA(:,1)',ZVA(:,2)',ZVA(:,3)',Wqnty(:,4),ones(length(Wreps),1),...
      Wstep,Wtype,PAfit,Pnfo,tol);
  end
else
  fprintf(['fitting all datafiles (#scans = ' num2str(sum(Wreps)) ') : '])
  msgold = ''; % old output text
  j0 = 0; % index offset indicating the start of the current scan
  Pnfo_org = Pnfo(1,[2 5]); % save the original setting for fit flag of phi0 and Mp
  for k=1:length(Wreps) % loop over each "new" scan (belonging to range s1)
    % check if only dealing with longitudinal data (then phi0 and Mphi should be kept fixed)
    if isequal(Wtype(k),0)
      Pnfo(1,[2 5]) = [0 0]; % flag phi0 and Mphi as fixed fitparameters
    else
      Pnfo(1,[2 5]) = Pnfo_org;
    end
    % loop over all repeats of the current "new" scan
    for j=1:Wreps(k)
      t = j0 + (1:Wstep(k)); % select range for "new" scan k, repeat j
      j0 = j0 + Wstep(k); % set index to next repeat scan for next round
      n = sum(Wreps(1:k-1)) + j; % scan number
      ind = Npar*(n-1) + (1:Npar); % indices of the fitparameters for "new" scan k, repeat j
      % call the fit algorithm
      [Vfit(t,1),Pchi(n),Pfit(ind),Perr(ind)] = LMfit_mpms(fitfun,ZV(t,1)',ZV(t,2)',ZV(t,3)',Wqnty(k,4),1,Wstep(k),...
        Wtype(k),Pfit(ind),Pnfo,tol); % can't use common fitparameters here
      % update output to screen (or edit box)
      if isequal(h,0)
        msgclr = repmat('\b',1,length(msgold));
        msgnew = ['scans completed = ' num2str(n)];
        fprintf(msgclr);
        fprintf(msgnew);
        msgold = msgnew;
      else
        h.String = num2str(100*n/sum(Wreps));
        pause(0.001); % need a tiny pause to update the figure
      end
    end
    % if repeat scans were used, also fit the average scans
    if isequal(RepScans,1)
      tA = sum(Wstep(1:k-1)) + (1:Wstep(k));
      indA = Npar*(k-1) + (1:Npar);
      [VAfit(tA,1),PAchi(k),PAfit(indA),PAerr(indA)] = LMfit_mpms(fitfun,ZVA(tA,1)',ZVA(tA,2)',ZVA(tA,3)',Wqnty(k,4),...
        1,Wstep(k),Wtype(k),PAfit(indA),Pnfo,tol);
    end
  end
end
fprintf('\n');

% if no repeat scans were used, the average fits are the single fits (required for plotting)
if isequal(RepScans,0)
  VAfit = Vfit;
  PAchi = Pchi;
  PAfit = Pfit;
  PAerr = Perr;
end

% allocated memory to store fitparameters in matrix form in a per scan (or "new" scan) fashion
a = zeros(sum(Wreps),2*Npar); % single fit [a aerr] per scan
aA = zeros(length(Wreps),5*Npar); % single fit [<a> <aerr> std(a)] + <repeat scan> fit [a aerr], all per "new" scan
for n=1:Npar
  for k=1:length(Wreps)
    for j=1:Wreps(k)
      ind = sum(Wreps(1:k-1))+j;
      a(ind,n) = Pfit((ind-1)*Npar+n);
      a(ind,n+Npar) = Perr((ind-1)*Npar+n);
    end
    ios = sum(Wreps(1:k-1)); % index offset for a to mark the start of "new" scan k
    aA(k,n+0*Npar) = mean(a(ios+(1:Wreps(k)),n));        % single scan: <a> (of its repeats)
    aA(k,n+1*Npar) = mean(a(ios+(1:Wreps(k)),n+Npar));   % single scan: <aerr> (of its repeats)
    aA(k,n+2*Npar) = std(a(ios+(1:Wreps(k)),n));         % single scan: std(a) (of its repeats)
    aA(k,n+3*Npar) = PAfit((k-1)*Npar+n);                % averaged repeat scans: a
    aA(k,n+4*Npar) = PAerr((k-1)*Npar+n);                % averaged repeat scans: aerr
  end
end
end