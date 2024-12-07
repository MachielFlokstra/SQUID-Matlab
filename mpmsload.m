function [ZV,ZVi,Wqnty,Wreps,Wstep,Wtype,FNS,NSi,fname] = mpmsload(fname,fdir)
% INFO
%  Loads the raw mpms datafiles and stores all relevant data in the output "W" variables.
%
% INPUT
%  fname : column vector listing all datafiles to load (should all be in the same dir)
%  fdir  : directory containing the datafiles
%   > can also call mpmsload without input arguments, then a dialog box appears to select files
%
% OUTPUT
%  ZV    : 3-column matrix with all measured data (Z,V,Verr)
%  ZVi   : 2-column matrix with the first and last index of ZV belonging to each scan
%  Wqnty : 5-column matrix with the quantity (H,<T>,<T>-error,Rotation,time) per "new" scan
%  Wreps : column vector with the number of scans per "new" scan
%  Wstep : column vector with the number of steps per "new" scan
%  Wtype : column vector with the coil-system used per "new" scan (0=longitudinal, 1=transverse)
%  FNS   : column vector with the total number of "new" scans per datafile
%  NSi   : column vector with the corresponding "new" scan for each scan
%  fname : list with all datafiles
%
% NOTE
%  A scan is measured along z from z=0 to z=scanlength. A scan is either a "new" scan or a repeat of the previous scan
%  to allow averaging over mulitple repeat scans to lower noise.
% ----------------------------------------------------------------------------------------------------------------------

% datafile constants
DataHead = 31;  % number of header rows in the data file
DataTS = 1;     % time stamp of the collected data
DataH = 3;      % field (Oe)
DataT1 = 4;     % start temperature (K)
DataT2 = 5;     % end temperature (K)
DataRep = 6;    % (repeat) scan number
DataRej = 7;    %#ok<NASGU> % rejected flag
DataZ = 8;      % position (cm)
DataR = 9;      % rotation position (deg)
DataVlo = 17;   % longitudinal response (scaled voltage)
DataVtr = 31;   % transverse response (scaled voltage)

% MPMS sensitivity ranges in emu (1 emu = 1E-3 Am^2)
MPMSsens = [1.25E-4 2.5E-4 6.25E-4 1.25E-3 2.5E-3 6.25E-3 1.25E-2 2.5E-2 6.25E-2 1.25E-1 2.5E-1 6.25E-1 1.25];
RangeLim = 1.00; % threshold fraction of sensitivity ranges (for higher signals range is increased)
RangeErr = 2.3*1E-3; % fractional error of the range

usr = getenv('USERNAME'); % get username (windows)
addpath(['C:\Users\' usr '\OneDrive - University of St Andrews\work\data analysis\PSI LEM\lemcore_2018']);
if isequal(nargin,0)
  % select datafile(s) - can select multiple files
  raw = pwd; % store current work directory
  % set directory to where all MPMS data is stored (if dir exists)
  if isequal(exist(['C:\Users\' usr '\OneDrive - University of St Andrews\work\data\MPMS\'],'dir'),7)
    cd(['C:\Users\' usr '\OneDrive - University of St Andrews\work\data\MPMS\'])
  end
  [fname,fdir,flag] = uigetfile('*.raw','MPMS Data Analysis Program', '', 'MultiSelect', 'on');
  if isequal(flag,0)
    fprintf('ERROR: no file selected\n')
    return
  end
  cd(raw); % restore work directory (in case it was changed)
end

% get number of files (the stored filenames are in a cell structure for multiple files, but a string for one file)
fnum = 1;
if iscell(fname)
  fnum = size(fname,2);
end
FNS = zeros(fnum,1);

% load all datafiles
Wsize = [0 0 0]; % totals of measured data points (index 1) "new" scans (index 2), total scans (index 3)
for k=1:fnum
  if isequal(fnum,1) % only one file
    %fprintf(['loading file: ' fname '\n']) % print current datafile to screen
    disp(['loading file: ' fname]);
    raw = csvread([fdir fname],DataHead); % read the raw data
  else % multiple files in a cell structure
    %fprintf(['loading file: ' fname{k} '\n']) % print current datafile to screen
    disp(['loading file: ' fname{k}]) % print current datafile to screen
    raw = csvread([fdir fname{k}],DataHead); % read the raw data
  end
  % determine the number of "new" scans (i.e. excluding repeat scans)
  Nscan = 0;
  for j=1:size(raw,1)
    if isequal(raw(j,DataZ),0) && isequal(raw(j,DataRep),1) % each "new" scan starts at 0 position
      Nscan = Nscan + 1;
    end
  end
  % determine the number of repeat scans (Nreps), number of steps (Nstep) and type (long/trans) per "new" scan
  Nreps = ones(Nscan,1); % the first scan is the "new" scan (Nreps = 1)
  Nstep = zeros(Nscan,1);
  Ntype = zeros(Nscan,1); % 0 = longitudinal, 1 = transverse
  n = 0; % "new" scan index
  for j=1:size(raw,1)
    if isequal(raw(j,DataZ),0) % marks start of a scan: either a "new" scan (DataRep=1) or a repeat scan
      if isequal(raw(j,DataRep),1)
        n = n + 1; % update "new" scans index
        % the number of steps and type for the previous scan is now know (j0 still marks the previous scan)
        if n > 1
          Nstep(n-1) = (j - j0)/Nreps(n-1); % number of steps for "new" scan n-1
          % determine the type (if longitudinal voltages are zero, then transverse coil setup was used)
          if isequal(sum(raw(j0:(j-1),DataVlo)),0)
            Ntype(n-1) = 1; % change type to transverse
          end
        end
        j0 = j; % save the index which marks the start of "new" scan n
      else
        Nreps(n) = Nreps(n) + 1; % increase repeat counter of "new" scan n
      end
    end
  end
  % add the steps and type for the last "new" scan
  Nstep(n) = (size(raw,1)+1-j0)/Nreps(n); %#ok<NASGU>
  if isequal(sum(raw(j0:size(raw,1),DataVlo)),0)
    Ntype(n) = 1; %#ok<NASGU>
  end
  FNS(k) = Nscan; % total number of "new" scans in current datafile  
  Wsize = Wsize + [size(raw,1) Nscan sum(Nreps)]; % update totals
  eval(['temp' num2str(k) ' = {raw, Nreps, Nstep, Ntype};']); % temporarily store data into temp1, temp2, temp3, ...
  fprintf([' >> measurement scans = ' num2str(sum(Nreps)) '\n']) % print details to screen
end

% allocate memory to contain all the data
ZV = zeros(Wsize(1),3);
ZVi = zeros(Wsize(3),2);
Wqnty = zeros(Wsize(2),5);
Wreps = zeros(Wsize(2),1);
Wstep = zeros(Wsize(2),1);
Wtype = zeros(Wsize(2),1);
NSi = zeros(Wsize(2),1);

% merge all data of the W# variables into single variables (WZV, WHTR, Wreps, Wstep and Wtype) 
Noffset = 0; % index offset for "per measurement" variables
Woffset = 0; % index offset for the WZH matrix
for k=1:fnum
  raw = eval(['temp' num2str(k) '{1}']); % get raw from the current datafile
  Nreps = eval(['temp' num2str(k) '{2}']); % get Nreps from the current datafile
  Wreps(Noffset+(1:length(Nreps))) = Nreps; % merge Nreps into Wreps
  Wstep(Noffset+(1:length(Nreps))) = eval(['temp' num2str(k) '{3}']); % merge Nstep into Wstep
  Wtype(Noffset+(1:length(Nreps))) = eval(['temp' num2str(k) '{4}']); % merge Ntype into Wtype
  DataV = raw(:,DataVlo) + raw(:,DataVtr); % scaled voltage (info about trans/long is in Wtype)
  ZV(Woffset+(1:size(raw,1)),1:2) = [raw(:,DataZ) DataV]; % merge Z,V data into ZV
  % for each measurement, get field, temperature and rotation and merge into Wqnty
  j0 = 1; % starting index of the first measurement
  for j=1:length(Nreps)
    t = 1:Wreps(Noffset+j)*Wstep(Noffset+j);
    Wqnty(Noffset+j,1) = raw(j0,DataH); % field
    Wqnty(Noffset+j,2) = mean([raw(j0-1+t,DataT1); raw(j0-1+t,DataT2)]); % avg temperature
    Wqnty(Noffset+j,3) = std([raw(j0-1+t,DataT1); raw(j0-1+t,DataT2)]); % temperature error
    Wqnty(Noffset+j,4) = raw(j0,DataR); % rotation
    Wqnty(Noffset+j,5) = raw(j0,DataTS); % time
    j0 = j0 + Wreps(Noffset+j)*Wstep(Noffset+j); % update index
  end
  % update offsets
  Woffset = Woffset + size(raw,1);
  Noffset = Noffset + length(Nreps);
  % clear large variable
  eval(['clear temp' num2str(k)])
end

% rescale
ZV(:,1) = ZV(:,1)*1E-2; % convert cm to m

% determine the data ranges for each scan and the scan to "new" scan conversion
n = 1; % scan index counter
j0 = 1; % index indicating the start of the current scan
for k=1:length(Wreps)
  for j=1:Wreps(k)
    ZVi(n,1:2) = [j0 j0+Wstep(k)-1];
    NSi(n) = k;
    n = n + 1; % update counter
    j0 = j0+Wstep(k); % set to next scan
  end
end
  
% (try to) determine the measurement error
j0 = 1; % index indicating the start of the current scan
for k=1:length(Wreps)
  for j=1:Wreps(k)
    raw = ZV(j0:j0+Wstep(k)-1,2); % the scaled voltage of the current scan
    tempmax = max(abs(raw)); % the maximum absolute value
    if isequal(Wtype(k),0)
      % for longitudinal scans a max voltage of 1.68V corresponds to 1 emu (purely Mz at r=0)
      tempmax = tempmax/1.68; % the corresponding emu
    else
      % for longitudinal scans a max voltage of 2.355V corresponds to 0.1 emu (purely Mx at r=phi=0)
      tempmax = tempmax/23.55;
    end
    % find the range used and set error
    RangeObtained = 0;
    RangeIndex = 1;
    while isequal(RangeObtained,0)
      if tempmax > RangeLim*MPMSsens(RangeIndex)
        RangeIndex = RangeIndex + 1; % increase index if max voltage above the current range limit
        % error check
        if RangeIndex > length(MPMSsens)
          fprintf([' WARNING: estimated moment = ' num2str(tempmax,4) ' emu, outside detection limit of 1.25 emu\n']);
          RangeObtained = 1;
          RangeIndex = length(MPMSsens);
        end
      else
        RangeObtained = 1;
      end
    end    
    ZV(j0:j0+Wstep(k)-1,3) = MPMSsens(RangeIndex)*RangeErr;
    j0 = j0 + Wstep(k); % set to next scan
  end
end
end