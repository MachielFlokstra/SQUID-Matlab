% ======================================================================================================================
% FIT PROGRAM FOR MPMS DATA
%  This program loads ".raw" data from the MPMS-XL system to analyses the measured scans obtained using DC measurements.
%  The sample is approximated by a magnetic point dipole with coordinate (r,phi,z) and magnetic moment (Mr,Mphi,Mz). The
%  central axis of the pickup coil system is the z-axis and the z-axis is also the rotation axis(*). The center of the
%  pickup coil system is taken as the origin. When phi=0, a vector along positive r points to the center of 3 of the
%  transverse coils while a vector along negative r direction points to the center of the other 3 transverse coils.
%
%  The mpms measures for a given rotation angle (rot) along the z-axis starting at zm=0 up to zm=scanlength (zm being
%  the z value that the mpms is using). Both rot and and zm have an uncalibrated zero and need to be shifted such that
%  the zm=0 value coincides with z=0 and rot=0 coincides with phi=0. The calibration values are z0 and phi0 and should
%  be determined from dedicated scans in order to optimize fitting of the signals. 
%   > Thus: phi=rot-phi0 and z=zm-z0, with phi and z the azimuth and z-position in the coordinate framework of the model
%                                     funtions and rot and zm their respective values in the mpms coordinate framework.
%
%  For accurate measurements, the measured scanlength should be symmetric around the center of the pickup coil origin
%  (i.e. ideally z0 equals half the scanlength). When using the transverse coils to detect the Mr,Mphi components, the
%  rotation can be set such to maximize the signal (which is when the net moments points right into a saddle coil).
%
%  The fit to the measurement data is made using derived model functions for the longitudinal and transverse pickup coil
%  geometries. There are 8 fitparameters per scan: r,phi0,z0,Mr,Mphi,Mz,c0,c1 where co and c1 are the offset and slope
%  of the linear background signal. Once (r,phi0,z0) are determined from dedicated scans, all further scans on the
%  sample can be fit with r,phi0,z0 fixed. Only for large temperature changes z0 might change due to contraction of the
%  measurement rod.
%
% (*) It appears the system alignment is not perfect and the rotation of a sample is not symmetric around the coil
%     center axis. The coordinates r and phi thus depend on the rotation angle in a non-trivial way and dedicated
%     rotation scans are required to find r(rot) and phi0(rot).
%
% LONGITUDINAL COILS
%  When using the longitudinal pickup coils, any Mphi component always generates zero net flux so Mphi (and phi0) cannot
%  be determined. Any Mr component generates an odd signal (odd in z) and the amplitude of this signal scales with both
%  Mr and r. This should leave both undetermined (but related to one another). However, due to imperfect alignment, the
%  amplitude of the response from Mz also scales with r and making a rotation scan can result in determining r(rot).
%  This might work best if Mz>>Mr. The shape of the signal due to Mz depends on r as well, but the effect is not so big
%  and only for r~4-5mm does it lead to a visually more spiked signal.
%
% TRANSVERSE COILS
%  The transverse pickup-coil setup is a bit more complex and uses 6 saddle-type coils (6 half cylindrical coils wrapped
%  around the mounting cylinder (where also the longitudinal coils are wrapped around). However, the resulting response
%  function to a transverse moment is (often) similar to that of a longitudinal moment on the longitudinal pickup-coil
%  setup. For the transverse measurements, any transverse (Mr or Mphi) component generates an even signal. The angle phi
%  determines the "phase" of the signal both look similar, but are shifted in phase: Mphi(phi-90) ~ Mr(phi). The radius
%  r modifies the the amplitude a little, but mainly deforms the Mr signal (broadening of the central peak). Any Mz
%  component generates an odd signal. For a perfectly centered point dipole, any Mz component generates zero net flux.
%
% STRATEGY FOR OBTAINING THE SAMPLE COORDINATE
%  - Set a moderate field to achieve Mz >> Mr,Mphi and |Mz| clearly measurable.
%  - Make a rotation scan (from 0 to 360 in 5-10 deg steps) using both pickup coils.
%  - Fit the longitudinal scans with common z0,Mz (and Mr=0) to obtain r(rot).
%  - Fit the transverse scans, applying r(rot),z0,Mz (and Mr=Mp=0), to obtain phi(rot) and correct for degeneracy.
%  - Plot the rotation path to check if solution plausible.
%  - Refit the longitudinal scans, applying r(rot),z0,Mz, to obtain Mr(rot).
%  - Refit the transverse scans, applying r(rot),phi(rot),z0,Mr(rot),Mz, to obtain Mp(rot).
%  - Check if Mr(rot) and Mp(rot) are plausible solutions.
%     > above steps are integrated in mpmsROT.m
%  - Now fit all other data applying r(rot),phi(rot),z0 (and use Mr(rot),Mp(rot),Mz as initial guesses)
%  - Can chose to make pair-wise fits for longitudinal and transverse data sets linking Mr,Mp,Mz between the two
%     > requires fitfun = mpmsfun2 (and FitAll=1)
%
% CALLS
%  LMfit_mpms.m - the Levenberg-Marquardt algorithm
%  mpmsfit.m - to initial fitparameters and call the fit algorithm
%  mpmsfitlim.m - to set the default fitparameter limits
%  mpmsfun.m - fitfunctions for the measurement signal
%  mpmsfun2.m - fitfunctions for the measurement signal (linked version)
%  mpmsload.m - to select the .raw datafiles to analyse
%  mpmsROT.m - to extract the rotation path from rotation scans
%
%                                                                                           Machiel Flokstra, 23/08/2018
% ======================================================================================================================
function squidfit(fname,fdir,fitfun,Pinit,Pnfo,FitSim,ROTr,ROTp,ROTmr,ROTmp)
% INPUT
%  fname  : row cell vector listing all datafiles to load (should all be in the same dir)
%  fdir   : directory containing the datafiles
%  fitfun : string with name of the fitfunction to use
%  Pinit  : row cell containing the initial values for r,phi0,z0,Mr,Mp,Mz,c0,c1
%           use Pinit(3)=0 to set z0 at half the scanlength
%           use Pinit(4)='max' to set Mr at estimated maximum (similar for Mp and Mz)
%  Pnfo   : 4x8 matrix containing all required fitparameter infos [fit flag; common flag; lower lim; upper lim]
%  FitSim : flag for simultaneous fitting // typically only use for a single rotation measurement file
%  ROTr   : 2-column matrix with [rot r] // if not empty, used to overwrite Pinit(1)
%  ROTp   : 2-column matrix with [rot phi0] // if not empty, used to overwrite Pinit(2)
%  ROTmr  : 2-column matrix with [rot Mr] // if not empty, used to overwrite Pinit(4)
%  ROTmp  : 2-column matrix with [rot Mp] // if not empty, used to overwrite Pinit(5)
%   > the ROT* variables should all be logically ordered with rot ranging from [0 to 360]
%
% ALTERNATIVE CALL OPTIONS
%  squidfit()  : "stand-alone version" where all settings and data loading is accessible via the main figure.
%  squidfit(1) : Use to skip the fitting routine and instead load the results from disk (to go straight to plotting).
%
% GLOBAL/SHARED VARIABLES - FOR CALLBACK USAGE
%  fname, Pinit, Pnfo : see above
%  ZV      : 3-column matrix with all measured data (Z,V,Verr)
%  ZVi     : 2-column matrix with the first and last index of ZV belonging to each scan
%  Wqnty   : 4-column matrix with the quantity (H,<T>,<T>-error,Rotation) per "new" scan
%  Wreps   : column vector with the number of scans per "new" scan
%  Wstep   : column vector with the number of steps per "new" scan
%  Wtype   : column vector with the coil-system used per "new" scan (0=longitudinal, 1=transverse)
%  FNS     : column vector with the total number of "new" scans per datafile
%  NSi     : column vector with the corresponding "new" scan for each scan
%  Vfit    : column vector for the obtained fits to each single scan
%  Pchi    : row vector with the chisquared of each fit belonging to Vfit
%  a       : matrix with single fit results "[a aerr]" per scan
%  ZVA     : 3-column matrix with the averaged data (avg version of ZV)
%  ZVAi    : 2-column matrix with the first and last index of ZVA belonging to each "new" scan
%  VAfit   : column vector for the obtained fits to each scan
%  PAchi   : row vector with the chisquared of each fit
%  aA      : matrix with the averages of the single fit results "[<a> <aerr> std(a)]" (first three blocks of 8 columns)
%           and the average repeat scan fit results (last two blocks of 8 columns)
%
% GLOBAL/SHARED OBJECT HANDLES
%  ha1     : axes to plot M
%  ha2     : axes to plot V    
%  ha3     : axes to plot a/aA
%  pop1    : popupmenu for M selection
%  pop2    : popupmenu for the datafile selection
%  bg      : buttongroup for the M component selection
%  tb1     : toggle button for showing std
%  tb2     : button for showing avg(fit) or fit(avg)
%  sli1    : slider for the scan selection
%  sli2    : slider for the fitparameter selection
%  edt1    : edit box to show the scan index
%  edt2    : edit box to show the fitparameter
%  edt3    : edit box to input chisquared limit
%  edt4-11 : edit boxes for the initial values of the fitparameters
%  edt12   : edit box for the fit progression
%  ch1-8   : checkboxes for the fit flags of the fitparameters
%  pb1     : load data button
%  pb2     : fit data button
% ----------------------------------------------------------------------------------------------------------------------

if isequal(nargin,0)
  mode = 2; % stand-alone version: loading/fitting is operated through the figure
  % The stand-alone version (no input parameters) will have limited functionality:
  %  > doesn't allow simulaneous fits nor linked fits
  %  > doesn't allow setting ROT* variables
  %  > doesn't allow overwriting fitparameter limits
  initclear(); % clear all only if ran as a stand-alone code (i.e. without input arguments)
  FitSim = 0;
  fitfun = 'mpmsfun';
  ROTr = [];
  ROTp = [];
  ROTmr = [];
  ROTmp = [];
  fname = ''; % declare as empty (fname doesn't exist yet, but is required initializing the figure)
  Pinit = zeros(1,8);
  Pnfo = zeros(4,8);
  Pnfo(3:4,:) = mpmsfitlim();
elseif isequal(nargin,1) && isequal(fname,1)
  mode = 3; % when called from outside to load previous results and create the figure to plot results
  fprintf('loading results from disk\n')
  load('FitResults')
else
  mode = 1; % when called from outside to load the listed files, make fits and create the figure to plot results
  [ZV,ZVi,Wqnty,Wreps,Wstep,Wtype,FNS,NSi,fname] = mpmsload(fname,fdir);
  [Vfit,Pchi,a,ZVA,ZVAi,VAfit,PAchi,aA] = mpmsfit(ZV,Wqnty,Wreps,Wstep,Wtype,fitfun,Pinit,Pnfo,FitSim,ROTr,ROTp,...
    ROTmr,ROTmp,0);
  save('FitResults','fname','ZV','ZVi','Wqnty','Wreps','FNS','NSi','Vfit','Pchi','a','ZVA','ZVAi','VAfit','PAchi','aA')
end

% ----------------------------------------------------------------------------------------------------------------------

% create interactive figure
h = figure;
if isequal(mode,2)
  set(h,'NumberTitle','off','Name','squidfit','Position',[100 100 1400 420],'Resize','off');
  pw = 4/14;  % plot panel width
  p1 = 0;     % panel-1 x-offset
  p2 = 4/14;  % panel-2 x-offset
  p3 = 8/14;  % panel-3 x-offset
  p4 = 12/14; % panel-4 x-offset
else
  set(h,'NumberTitle','off','Name','squidfit','Position',[100 100 1200 420],'Resize','off');
  pw = 1/3; % plot panel width
  p1 = 0;   % panel-1 x-offset
  p2 = 1/3; % panel-2 x-offset
  p3 = 2/3; % panel-3 x-offset
end
BGC = [0.84 0.89 1]; % background color for panels with interactive bits
HLC = [0 0 0]; % highlight color for panel bounding boxes
TBC = [1 1 0.55]; % color for toggle buttons
PBC = [146 256 96]/256; % color for push buttons

% create the panels in which to place axes
pan = uipanel(h,'position',[p1 0.1428 pw 0.8572]); % position = [x0,y0,width,height] as a fraction of the parent size
ha1 = axes(pan); % axes to plot M
pan = uipanel(h,'position',[p2 0.1428 pw 0.8572]);
ha2 = axes(pan); % axes to plot V + fit
pan = uipanel(h,'position',[p3 0.1428 pw 0.8572]);
ha3 = axes(pan); % axes to plot the fitparameters

% create a panel below each axes to add info and interactive parts (buttongroups, sliders, etc)

% panel for the M plot
pan = uipanel(h,'position',[p1 0.0714 pw 0.0714]);
pan.BackgroundColor = BGC;
pan.HighlightColor = HLC;
tex = uicontrol(pan,'style','text','position',[5 2 70 20],'HorizontalAlignment','left','FontSize',9);
tex.String = 'Select Plot';
tex.BackgroundColor = BGC;
pop1 = uicontrol(pan,'style','popupmenu','position',[75 6 56 20],'FontSize',9);
pop1.String = {'M(T)','M(H)','M(rot)'}; % popup menu options
bg = uibuttongroup(pan,'position',[0.36 0 0.50 4]); % buttongroup to place mutually exclusive radiobuttons
bg.BorderType = 'none';
bg.BackgroundColor = BGC;
uicontrol(bg,'style','radiobutton','position',[0 5 45 20],'string','Mr','FontSize',9,'BackgroundColor',BGC);
uicontrol(bg,'style','radiobutton','position',[43 5 45 20],'string','Mphi','FontSize',9,'BackgroundColor',BGC);
uicontrol(bg,'style','radiobutton','position',[100 5 45 20],'string','Mz','FontSize',9,'BackgroundColor',BGC);
uicontrol(bg,'style','radiobutton','position',[144 5 45 20],'string','|M|','FontSize',9,'BackgroundColor',BGC);
tb1 = uicontrol(pan,'style','togglebutton','position',[335 5 52 18],'FontSize',9); % toggle button
tb1.String = 'add std'; % toggle to showing standard deviation of M (if repeat scans are used)
tb1.BackgroundColor = TBC;

% panel for the V plot
pan = uipanel(h,'position',[p2 0.0714 pw 0.0714]);
pan.BackgroundColor = BGC;
pan.HighlightColor = HLC;
tex = uicontrol(pan,'style','text','position',[5 2 70 20],'HorizontalAlignment','left','FontSize',9);
tex.String = 'Select Scan';
tex.BackgroundColor = BGC;
sli1 = uicontrol(pan,'style','slider','position',[90 4 220 20]);
sli1.Min = 1; % slider min value (max, value and sliderstep are correctly set during callback functions)
sli1.Max = 2; % slider max value: set correctly during callback functions
sli1.Value = 1; % slider initial value (set to first scan)
sli1.SliderStep = [1/(sli1.Max-1) 1]; % slider steps (arrow & trough): set correctly during callback functions
sli1.BackgroundColor = [1 1 1];
edt1 = uicontrol(pan,'style','edit','position',[330 5 50 18],'FontSize',9); % edit box to show scan number
edt1.Enable = 'off';
edt1.String = num2str(sli1.Min);

% panel for the fitparameters plot
pan = uipanel(h,'position',[p3 0.0714 pw 0.0714]);
pan.BackgroundColor = BGC;
pan.HighlightColor = HLC;
tex = uicontrol(pan,'style','text','position',[5 2 120 20],'HorizontalAlignment','left','FontSize',9);
tex.String = 'Select Fitparameter';
tex.BackgroundColor = BGC;
sli2 = uicontrol(pan,'style','slider','position',[130 4 180 20]);
sli2.Min = 1; % slider min value
sli2.Max = 9; % slider max value = total number of fitparameters (=Npar=8) + 1 (for chi2)
sli2.Value = 1; % slider initial value (set to first fitparameter)
sli2.SliderStep = [1/(sli2.Max-1) 1]; % slider steps (arrow & trough)
sli2.BackgroundColor = [1 1 1];
edt2 = uicontrol(pan,'style','edit','position',[330 5 50 18],'FontSize',9); % edit box to show fitparameter
edt2.Enable = 'off';
edt2.String = 'Rr';

% create a bottom panel for datafile select, chi2 threshold, and toggle between single/average scans
pan = uipanel(h,'position',[0 0 3*pw 0.0714]);
pan.BackgroundColor = BGC;
pan.HighlightColor = HLC;
tex = uicontrol(pan,'style','text','position',[5 2 70 20],'HorizontalAlignment','left','FontSize',9);
tex.String = 'Select File';
tex.BackgroundColor = BGC;
pop2 = uicontrol(pan,'style','popupmenu','position',[75 6 685 20],'FontSize',9);
pop2.String = fname; % set all the filenames as the popup menu options
tex = uicontrol(pan,'style','text','position',[765 2 50 20],'HorizontalAlignment','left','FontSize',9);
tex.String = 'chi2 lim';
tex.BackgroundColor = BGC;
edt3 = uicontrol(pan,'style','edit','position',[815 5 50 18],'FontSize',9); % edit box for chisquared lim input
edt3.String = '0';
tb2 = uicontrol(pan,'style','togglebutton','position',[880 4 300 20],'FontSize',9);
tb2.String = 'showing: single scan fits + averaged fitparameters';
tb2.BackgroundColor = TBC;
tb2.ForegroundColor = [1 0 0];

% ----------------------------------------------------------------------------------------------------------------------

% create the extra bits for mode 2
if isequal(mode,2)
  % top panel
  pan = uipanel(h,'position',[p4 0.1428 pw/2 0.8572],'BackgroundColor',BGC,'HighlightColor',HLC);
  % coordinate
  tex = uicontrol(pan,'style','text','position',[10 330 150 20],'HorizontalAlignment','center','FontSize',10);
  set(tex,'String','COORDINATE','FontWeight','bold','BackgroundColor',BGC)
  % r 
  tex = uicontrol(pan,'style','text','position',[10 300 35 20],'HorizontalAlignment','left','FontSize',10);
  set(tex,'String','r','BackgroundColor',BGC);
  tex = uicontrol(pan,'style','text','position',[120 300 35 20],'HorizontalAlignment','left','FontSize',10);
  set(tex,'String','mm','BackgroundColor',BGC);
  edt4 = uicontrol(pan,'style','edit','position',[45 302 70 20],'FontSize',10,'String',num2str(1));
  ch1 = uicontrol(pan,'style','checkbox','position',[170 302 20 20],'Value',0,'BackgroundColor',BGC);
  % phi0
  tex = uicontrol(pan,'style','text','position',[10 270 35 20],'HorizontalAlignment','left','FontSize',10);
  set(tex,'String','phi0','BackgroundColor',BGC);
  tex = uicontrol(pan,'style','text','position',[120 270 35 20],'HorizontalAlignment','left','FontSize',10);
  set(tex,'String','deg','BackgroundColor',BGC);
  edt5 = uicontrol(pan,'style','edit','position',[45 272 70 20],'FontSize',10,'String',num2str(0));
  ch2 = uicontrol(pan,'style','checkbox','position',[170 272 20 20],'Value',0,'BackgroundColor',BGC);
  % z0
  tex = uicontrol(pan,'style','text','position',[10 240 35 20],'HorizontalAlignment','left','FontSize',10);
  set(tex,'String','z0','BackgroundColor',BGC);
  tex = uicontrol(pan,'style','text','position',[120 240 35 20],'HorizontalAlignment','left','FontSize',10);
  set(tex,'String','cm','BackgroundColor',BGC);
  edt6 = uicontrol(pan,'style','edit','position',[45 242 70 20],'FontSize',10,'String','i');
  ch3 = uicontrol(pan,'style','checkbox','position',[170 242 20 20],'Value',0,'BackgroundColor',BGC);
  % moment
  tex = uicontrol(pan,'style','text','position',[10 200 150 20],'HorizontalAlignment','center','FontSize',10);
  set(tex,'String','MOMENT','FontWeight','bold','BackgroundColor',BGC)
  % Mr
  tex = uicontrol(pan,'style','text','position',[10 170 35 20],'HorizontalAlignment','left','FontSize',10);
  set(tex,'String','Mr','BackgroundColor',BGC);
  tex = uicontrol(pan,'style','text','position',[120 170 35 20],'HorizontalAlignment','left','FontSize',10);
  set(tex,'String','emu','BackgroundColor',BGC);
  edt7 = uicontrol(pan,'style','edit','position',[45 172 70 20],'FontSize',10,'String',num2str(0));
  ch4 = uicontrol(pan,'style','checkbox','position',[170 172 20 20],'Value',0,'BackgroundColor',BGC);
  % Mphi
  tex = uicontrol(pan,'style','text','position',[10 140 35 20],'HorizontalAlignment','left','FontSize',10);
  set(tex,'String','Mphi','BackgroundColor',BGC);
  tex = uicontrol(pan,'style','text','position',[120 140 35 20],'HorizontalAlignment','left','FontSize',10);
  set(tex,'String','emu','BackgroundColor',BGC);
  edt8 = uicontrol(pan,'style','edit','position',[45 142 70 20],'FontSize',10,'String',num2str(0));
  ch5 = uicontrol(pan,'style','checkbox','position',[170 142 20 20],'Value',0,'BackgroundColor',BGC);
  % Mz
  tex = uicontrol(pan,'style','text','position',[10 110 35 20],'HorizontalAlignment','left','FontSize',10);
  set(tex,'String','Mz','BackgroundColor',BGC);
  tex = uicontrol(pan,'style','text','position',[120 110 35 20],'HorizontalAlignment','left','FontSize',10);
  set(tex,'String','emu','BackgroundColor',BGC);
  edt9 = uicontrol(pan,'style','edit','position',[45 112 70 20],'FontSize',10,'String','i');
  ch6 = uicontrol(pan,'style','checkbox','position',[170 112 20 20],'Value',1,'BackgroundColor',BGC);
  % background
  tex = uicontrol(pan,'style','text','position',[10 70 150 20],'HorizontalAlignment','center','FontSize',10);
  set(tex,'String','BACKGROUND','FontWeight','bold','BackgroundColor',BGC)
  % c0
  tex = uicontrol(pan,'style','text','position',[10 40 35 20],'HorizontalAlignment','left','FontSize',10);
  set(tex,'String','c0','BackgroundColor',BGC);
  tex = uicontrol(pan,'style','text','position',[120 40 35 20],'HorizontalAlignment','left','FontSize',10);
  set(tex,'String','V','BackgroundColor',BGC);
  edt10 = uicontrol(pan,'style','edit','position',[45 42 70 20],'FontSize',10,'String',num2str(0));
  ch7 = uicontrol(pan,'style','checkbox','position',[170 42 20 20],'Value',1,'BackgroundColor',BGC);
  % c0
  tex = uicontrol(pan,'style','text','position',[10 10 35 20],'HorizontalAlignment','left','FontSize',10);
  set(tex,'String','c1','BackgroundColor',BGC);
  tex = uicontrol(pan,'style','text','position',[120 10 35 20],'HorizontalAlignment','left','FontSize',10);
  set(tex,'String','V/cm','BackgroundColor',BGC);
  edt11 = uicontrol(pan,'style','edit','position',[45 12 70 20],'FontSize',10,'String',num2str(0));
  ch8 = uicontrol(pan,'style','checkbox','position',[170 12 20 20],'Value',1,'BackgroundColor',BGC);
  % mid panel
  pan = uipanel(h,'position',[p4 0.0714 pw/2 0.0714],'BackgroundColor',BGC,'HighlightColor',HLC);
  pb1 = uicontrol(pan,'style','pushbutton','position',[10 4 80 20],'FontSize',10);
  set(pb1,'String','LOAD','FontWeight','bold','BackgroundColor',PBC);
  pb2 = uicontrol(pan,'style','pushbutton','position',[110 4 80 20],'FontSize',10);
  set(pb2,'String','FIT','FontWeight','bold','BackgroundColor',PBC,'Enable','off');
  % bot panel
  pan = uipanel(h,'position',[p4 0 pw/2 0.0714],'BackgroundColor',BGC,'HighlightColor',HLC);
  tex = uicontrol(pan,'style','text','position',[10 2 105 20],'HorizontalAlignment','left','FontSize',9);
  set(tex,'String','Scans completed:','BackgroundColor',BGC);
  edt12 = uicontrol(pan,'style','edit','position',[115 5 55 20],'FontSize',9,'String',num2str(0));
  tex = uicontrol(pan,'style','text','position',[172 2 20 20],'HorizontalAlignment','left','FontSize',9);
  set(tex,'String','%','BackgroundColor',BGC);
end

% ----------------------------------------------------------------------------------------------------------------------

% Attach callback functions to objects

pop1.Callback = {@Mselect};
bg.SelectionChangedFcn = {@Mselect};
tb1.Callback = {@Mselect};
sli1.Callback = {@Vselect};
sli2.Callback = {@Pselect};
pop2.Callback = {@Fselect};
tb2.Callback = {@Fselect};
edt3.Callback = {@Fselect};
if isequal(mode,2)
  pb1.Callback = {@loaddata};
  pb2.Callback = {@fitdata,fitfun,FitSim,ROTr,ROTp,ROTmr,ROTmp};
end

% disable plot manipulation objects if running as stand-alone, else make the plots
if isequal(mode,2)
  PlotManip([],[],'off');
else
  Fselect([],[]);
end

% ----------------------------------------------------------------------------------------------------------------------

% Declaration of callback functions.
%  > Callback function don't generate output. Any variables changes (or declarations) made during a callback that are
%    required to keep for future callbacks should be of "global" type, shared between all functions.
%  > The first two inputs are generated by default when the callback function is triggered by an object. The first input
%    is the responsible object handle and the second the event type.
%  > Callback functions can also be called "manually" (i.e. not by a trigger) in the standard way to call a function. In
%    that case the inputs 1 and 2 can be used to provide additional information.

% ----------------------------------------------------------------------------------------------------------------------

  % Callback function for plotting M (apply chisquared threshold to plot)
  function Mselect(h,~)
    % check if toggle button was pressed to get here
    if ~isempty(h)
      if isequal(h.Type,'uicontrol')
        if isequal(h.Style,'togglebutton')
          % change text
          if isequal(h.Value,1)
            h.String = 'remove';
          else
            h.String = 'add std';
          end
        end
      end
    end
    % determine which plot to make and which datafile to use
    Mtype = pop1.Value; % 1=M(T), 2=M(H), 3=M(rot)
    for MtypeInd=1:4
      if isequal(bg.Children(MtypeInd).Value,1)
        switch bg.Children(MtypeInd).String
          case 'Mr'
            Mcomp = 1;
          case 'Mphi'
            Mcomp = 2;
          case 'Mz'
            Mcomp = 3;
          case '|M|'
            Mcomp = 4;
        end
      end
    end
    dataf = pop2.Value; % index number of the datafile to use for plotting
    % select the row index
    switch Mtype
      case 1
        row = 2;
        xlab = 'temperature (K)';
      case 2
        row = 1;
        xlab = 'applied field (Oe)';
      case 3
        row = 4;
        xlab = 'rotation (deg)';
    end
    % determine the data range to use
    ros = sum(FNS(1:(dataf-1))); % index offset for Wreps,Wqnty and aA to mark the start of current datafile
    ri = ros+(1:FNS(dataf)); % corresponding range of index numbers for current datafile
    % ignore any index with chisquared above threshold
    chi2lim = str2double(edt3.String); % chisquared threshold
    % take out indices with chisquared above threshold
    if ~isequal(chi2lim,0)
      for k=length(ri):-1:1
        if PAchi(ri(k)) > chi2lim
          ri(k) = []; % take the index out
        end
      end
    end    
    % avg(fit) or fit(avg)
    aos = 0; % index offset for aA marking the start of [<a> <aerr> sta(a)] (aos=0) or [a aerr] (aos=3*Npar=24)
    lc = '.r-'; % line color
    if isequal(tb2.Value,1)
      aos = 24; % tb2 pressed: plot the fits of the repeat scan averages
      lc = '.b-'; % line color    
    end
    % make the V(z) plot and highlight the datapoint belonging to the scan in the M plot
    sind = round(sli1.Value); % index number of the scan (or "new" scan) of the current datafile to plot
    if isequal(tb2.Value,1)
      % tb2 pressed: sind is already the "new" scan index (but not yet the total/cummulative value)
      hp = ros + sind; % cummulative "new" scan index
    else
      % tb.2 nun-pressed: sind is the single scan index => need conversion to "new" scan index
      hp = sum(Wreps(1:ros)) + sind; % the cummulative scan index of the datapoint to highlight
      hp = NSi(hp); % the corresponding "new" scan index
    end
    axes(ha1); % make the axes active
    switch Mcomp
      case 1
        % plot standard deviation of repeats if button pressed
        if isequal(tb1.Value,1)
          x = Wqnty(ri,row)';
          y = 1E3*aA(ri,4)';
          ysig = 1E3*aA(ri,16+4)';
          plot([x; x],[y-ysig; y+ysig],'k'); % first plot in case of tb1 pressed
          hold on % set next plots to add (the first plot needs to erase the previous one)
        end
        % plot data on top of error bars (to be able to select them)
        plot(gca,Wqnty(ri,row),1E3*aA(ri,aos+4),lc); % first plot in case of tb1 unpressed
        hold on % set next plots to add (the first plot needs to erase the previous one)
        % add highlight only if below threshold (else the M point won't be in the plot)
        if (PAchi(hp) < chi2lim) || isequal(chi2lim,0)
          plot(gca,Wqnty(hp,row),1E3*aA(hp,aos+4),'s','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
        end
        % set Y-axes limits in case of constant, but non-zero y values without error bars
        if ~isequal(tb1.Value,1) && std(aA(ri,aos+4))<1E-10
          y = 1E3*mean(aA(ri,aos+4));
          if isequal(y,0)
            ha1.YLim = [-1 1];
          else
            ha1.YLim = y + [-0.1 0.1]*abs(y);
          end
        end
        hold off % set axes to erase when plotting next time
      case 2
        if isequal(tb1.Value,1)
          x = Wqnty(ri,row)';
          y = 1E3*aA(ri,5)';
          ysig = 1E3*aA(ri,16+5)';
          plot([x; x],[y-ysig; y+ysig],'k');
          hold on
        end
        plot(gca,Wqnty(ri,row),1E3*aA(ri,aos+5),lc);
        hold on
        if (PAchi(hp) < chi2lim) || isequal(chi2lim,0)
          plot(gca,Wqnty(hp,row),1E3*aA(hp,aos+5),'s','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
        end
        if ~isequal(tb1.Value,1) && std(aA(ri,aos+5))<1E-10
          y = 1E3*mean(aA(ri,aos+5));
          if isequal(y,0)
            ha1.YLim = [-1 1];
          else
            ha1.YLim = y + [-0.1 0.1]*abs(y);
          end
        end
        hold off
      case 3
        if isequal(tb1.Value,1)
          x = Wqnty(ri,row)';
          y = 1E3*aA(ri,6)';
          ysig = 1E3*aA(ri,16+6)';
          plot([x; x],[y-ysig; y+ysig],'k');
          hold on
        end
        plot(gca,Wqnty(ri,row),1E3*aA(ri,aos+6),lc);
        hold on
        if (PAchi(hp) < chi2lim) || isequal(chi2lim,0)
          plot(gca,Wqnty(hp,row),1E3*aA(hp,aos+6),'s','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
        end
        % set Y-axes limits in case of constant, but non-zero y values without error bars
        if ~isequal(tb1.Value,1) && std(aA(ri,aos+6))<1E-10
          y = 1E3*mean(aA(ri,aos+6));
          if isequal(y,0)
            ha1.YLim = [-1 1];
          else
            ha1.YLim = y + [-0.1 0.1]*abs(y);
          end
        end
        hold off
      case 4
        if isequal(tb1.Value,1)
          x = Wqnty(ri,row)';
          y = 1E3*sqrt(aA(ri,4).^2 + aA(ri,5).^2 + aA(ri,6).^2)';
          ysig = 1E3*sqrt(aA(ri,16+4).^2 + aA(ri,16+5).^2 + aA(ri,16+6).^2)';
          plot([x; x],[y-ysig; y+ysig],'k');
          hold on
        end
        plot(gca,Wqnty(ri,row),1E3*sqrt(aA(ri,aos+4).^2 + aA(ri,aos+5).^2 + aA(ri,aos+6).^2),lc);
        hold on
        if (PAchi(hp) < chi2lim) || isequal(chi2lim,0)
          plot(gca,Wqnty(hp,row),1E3*sqrt(aA(hp,aos+4).^2 + aA(hp,aos+5).^2 + aA(hp,aos+6).^2),'s',...
            'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
        end
        hold off
    end
    ha1.Title.String = 'Derived Moment';
    ha1.YLabel.String = 'moment (emu)';
    ha1.XLabel.String = xlab;
  end

% ----------------------------------------------------------------------------------------------------------------------
 
  % Callback function for plotting V (always plot - so one can see why chisquared was bad)
  function Vselect(~,~)
    dataf = pop2.Value; % index number of the datafile to use for plotting
    sind = round(sli1.Value); % index number of the scan (or "new" scan) of the current datafile to plot
    % > sind is scrolling over all scans or over "new" scans (which is setup in Fselect())
    edt1.String = num2str(sind); % udpate the edit box
    axes(ha2); % make the axes active
    % avg(fit) or fit(avg)
    if isequal(tb2.Value,1)
      % tb2 pressed: plot the fits of the repeat scan averages
      ros = sum(FNS(1:(dataf-1))); % index offset for VZAi to mark the start of current datafile
      plot(gca,1E2*ZVA(ZVAi(ros+sind,1):ZVAi(ros+sind,2),1),VAfit(ZVAi(ros+sind,1):ZVAi(ros+sind,2)),'k-');
      hold on
      plot(gca,1E2*ZVA(ZVAi(ros+sind,1):ZVAi(ros+sind,2),1),ZVA(ZVAi(ros+sind,1):ZVAi(ros+sind,2),2),'ob');
      hold off
    else
      % tb2 un-pressed: plot the single scan fits
      ros = sum(Wreps(1:sum(FNS(1:(dataf-1))))); % index offset for VZi to mark the start of current datafile
      plot(gca,1E2*ZV(ZVi(ros+sind,1):ZVi(ros+sind,2),1),Vfit(ZVi(ros+sind,1):ZVi(ros+sind,2)),'k-');
      hold on
      plot(gca,1E2*ZV(ZVi(ros+sind,1):ZVi(ros+sind,2),1),ZV(ZVi(ros+sind,1):ZVi(ros+sind,2),2),'or');
      hold off
    end
    ha2.Title.String = 'Measured Signal + Fit';
    ha2.XLabel.String = 'z (cm)'; ha2.YLabel.String = 'scaled voltage (V)';
    % replot M to update the highlighted datapoint
    Mselect([],[]);
  end

% ----------------------------------------------------------------------------------------------------------------------

  % Callback function for plotting fitparameters (apply chisquared threshold to plot)
  function Pselect(~,~)
    dataf = pop2.Value; % index number of the datafile to use for plotting
    pars = {'r' 'phi0' 'z0' 'Mr' 'Mphi' 'Mz' 'c0' 'c1' 'chi2'};
    ylab = {'r (mm)' 'phi0 (deg)' 'z0 (cm)' 'Mr (emu)' 'Mphi (emu)' 'Mz (emu)' 'c0 (V)' 'c1 (V/cm)' 'chisquared'};
    mult = [1E3 1 1E2 1E3 1E3 1E3 1 1E-2 1]; % multipliers required to plot in desired units
    pari = round(sli2.Value); % index number of the fitparameter to plot
    edt2.String = pars{pari}; % update the edit box
    chi2lim = str2double(edt3.String); % chisquared threshold
    axes(ha3); % make the axes active
    % avg(fit) or fit(avg)
    if isequal(tb2.Value,1)
      % tb2 pressed: plot the fits of the repeat scan averages
      ros = sum(FNS(1:(dataf-1))); % index offset for aA to mark the start of current datafile
      ri = ros+(1:FNS(dataf)); % corresponding range of index numbers for current datafile
      % take out indices with chisquared above threshold
      if ~isequal(chi2lim,0)
        for k=length(ri):-1:1
          if PAchi(ri(k)) > chi2lim
            ri(k) = []; % take the index out
          end
        end
      end
      if isequal(pari,9)
        plot(gca,ri-ros,mult(pari)*PAchi(ri),'.b-');
      elseif isequal(pari,2)
        % plot phi instead of phi0 (which is perhaps much more useful)
        y = conj(Wqnty(ri,4)' - aA(ri,24+pari)'); % a (of fits of repeat scan averages)
        yerr = conj(aA(ri,32+pari)'); % aerr (of fits of repeat scan averages)
        plot(gca,[ri-ros; ri-ros],mult(pari)*[y-yerr; y+yerr],'-k');
        hold on;
        plot(gca,ri-ros,mult(pari)*(Wqnty(ri,4)-aA(ri,24+pari)),'.b-'); % a (of fits of repeat scan averages)
        hold off;
        % change ylabel
        ylab{2} = 'phi = rot - phi0 (deg)';
      else
        y = conj(aA(ri,24+pari)'); % a (of fits of repeat scan averages)
        yerr = conj(aA(ri,32+pari)'); % aerr (of fits of repeat scan averages)
        plot(gca,[ri-ros; ri-ros],mult(pari)*[y-yerr; y+yerr],'-k');
        hold on;
        plot(gca,ri-ros,mult(pari)*aA(ri,24+pari),'.b-'); % a (of fits of repeat scan averages)
        hold off;
      end
    else
      % tb2 un-pressed: plot the single scan fits
      ros = sum(Wreps(1:sum(FNS(1:(dataf-1))))); % index offset for a to mark the start of current datafile
      ri = ros+(1:sum(Wreps( sum(FNS(1:(dataf-1))) + (1:FNS(dataf)) )));  % corresponding range of index
      % take out indices with chisquared above threshold
      if ~isequal(chi2lim,0)
        for k=length(ri):-1:1
          if Pchi(ri(k)) > chi2lim
            ri(k) = []; % take the index out
          end
        end
      end
      if isequal(pari,9)
        plot(gca,ri-ros,mult(pari)*Pchi(ri),'.r-');
      else
        y = conj(a(ri,pari)');
        yerr = conj(a(ri,8+pari)'); % a (of single scan fits)
        plot(gca,[ri-ros; ri-ros],mult(pari)*[y-yerr; y+yerr],'-k'); % aerr (of single scan fits)
        hold on;
        plot(gca,ri-ros,mult(pari)*a(ri,pari),'.r-'); % a (of single scan fits)
        hold off;
      end
    end
    % use full x-axis (unless only one datapoint to show)
    if ~isequal(length(ri),1)
      set(gca,'XLim',[ri(1) ri(length(ri))]-ros);
    end
    % set labels
    ha3.Title.String = 'Fitparameter Values';
    if isequal(tb2.Value,1)
      ha3.XLabel.String = 'new scan number';
    else
      ha3.XLabel.String = 'scan number';
    end
    ha3.YLabel.String = ylab(pari);
  end

% ----------------------------------------------------------------------------------------------------------------------

  % Callback function to be used when selecting a new datafile of toggling the avg(fit)/fit(avg) button
  function Fselect(h,~)
    if ~isempty(h)
      if isequal(h.Style,'togglebutton')
        % tb2 presses/un-pressed
        if isequal(h.Value,1)
          tb2.String = 'showing: results of fitting repeat scan averages';
          tb2.ForegroundColor = [0 0 1];
          % tb1 (for adding/removing the std of a) has no function now, set it to not add std and disable
          tb1.Value = 0;
          tb1.Enable = 'off';
          tb1.String = 'add std';
        else
          tb2.String = 'showing: single scan fits + averaged fitparameters';
          tb2.ForegroundColor = [1 0 0];
          tb1.Enable = 'on'; % tb1 can be used again
        end
      end
    end
    % setting to plot a of averaged scans, can not plot std, so grey-out the button and set to off in that case!
    dataf = pop2.Value; % index number of the datafile to use for plotting
    % reset the select scan slider (needs to be done for both togglebutton change and file change)
    sli1.Value = 1; % set to first scan
    if isequal(tb2.Value,1)
      % tb2 pressed: scroll over "new" scans
      sli1.Max = FNS(dataf); % total number of "new" scans in currently selected datafile
    else
      % tb2 un-pressed: scroll over single scans
      roff = sum(FNS(1:(dataf-1))); % index offset for Wreps to mark the start of current datafile
      sli1.Max = sum(Wreps(roff+(1:FNS(dataf)))); % total number of scans in currently selected datafile
    end
    if isequal(sli1.Max,1)
      % if only 1 scan to show, disable slider
      sli1.Enable = 'off';
    else
      sli1.Enable = 'on'; % enable slider (in case it was disabled)
      sli1.SliderStep = [1/(sli1.Max-1) 1]; % set slider steps
    end
    % replot
    Vselect([],[]);
    Pselect([],[]);
  end

% ----------------------------------------------------------------------------------------------------------------------

  % Callback function to load data
  function loaddata(~,~)
    [ZV,ZVi,Wqnty,Wreps,Wstep,Wtype,FNS,NSi,fname] = mpmsload();
    PlotManip([],[],'off'); % disable plot manipulations
    InputManip([],[],'on'); % enable fit input manipulation
    pop2.String = fname; % update the menu options of the select file dropbox
    pop2.Value = 1; % set to first file in the list
  end

% ----------------------------------------------------------------------------------------------------------------------

  % Callback function to fit data
  function fitdata(~,~,fitfun,FitAll,ROTr,ROTp,ROTmr,ROTmp)
    % disable all inputs during calculation
    InputManip([],[],'off');
    PlotManip([],[],'off');
    % load initial fitparameter values from figure
    Pinit(1) = 1E-3*str2double(edt4.String);
    Pinit(2) = str2double(edt5.String);
    Pinit(3) = str2double(edt6.String);
    if ~isequal(Pinit(3),1i)
      Pinit(3) = 1E-2*Pinit(3);
    end
    Pinit(4) = str2double(edt7.String);
    Pinit(5) = str2double(edt8.String);
    Pinit(6) = str2double(edt9.String);
    for k=[4 5 6]
      if ~isequal(Pinit(k),1i)
        Pinit(k) = 1E-3*Pinit(k);
      end
    end
    Pinit(7) = str2double(edt10.String);
    Pinit(8) = 1E2*str2double(edt11.String);
    % load fit flags from figure
    Pnfo(1,1) = ch1.Value;
    Pnfo(1,2) = ch2.Value;
    Pnfo(1,3) = ch3.Value;
    Pnfo(1,4) = ch4.Value;
    Pnfo(1,5) = ch5.Value;
    Pnfo(1,6) = ch6.Value;
    Pnfo(1,7) = ch7.Value;
    Pnfo(1,8) = ch8.Value;
    % make the fits
    [Vfit,Pchi,a,ZVA,ZVAi,VAfit,PAchi,aA] = mpmsfit(ZV,Wqnty,Wreps,Wstep,Wtype,fitfun,Pinit,Pnfo,FitAll,ROTr,ROTp,...
      ROTmr,ROTmp,edt12);
    % ensable all inputs after calculation
    InputManip([],[],'on');
    PlotManip([],[],'on');
    % plot
    Fselect([],[]);
  end

% ----------------------------------------------------------------------------------------------------------------------

  % Callback function to set enable state of input objects for plotting
  function PlotManip(~,~,state)
    % state = 'on' to enable objects, state = 'off' to disable objects
    pop1.Enable = state;
    for k=1:4
      bg.Children(k).Enable = state;
    end
    tb1.Enable = state;
    sli1.Enable = state;
    sli2.Enable = state;
    pop2.Enable = state;
    tb2.Enable = state;
    edt3.Enable = state;
  end

% ----------------------------------------------------------------------------------------------------------------------

  % Callback function to set enable state of input objects for fitting
  function InputManip(~,~,state)
    % state = 'on' to enable objects, state = 'off' to disable objects
    edt4.Enable = state;
    edt5.Enable = state;
    edt6.Enable = state;
    edt7.Enable = state;
    edt8.Enable = state;
    edt9.Enable = state;
    edt10.Enable = state;
    edt11.Enable = state;
    ch1.Enable = state;
    ch2.Enable = state;
    ch3.Enable = state;
    ch4.Enable = state;
    ch5.Enable = state;
    ch6.Enable = state;
    ch7.Enable = state;
    ch8.Enable = state;
    pb1.Enable = state;
    pb2.Enable = state;
  end
end

% ======================================================================================================================
function initclear()
close all; % close all figures
clearvars; % clear variable and globals from workspace, clear functions and mex links
clc; % clear command window
warning('off', 'all'); % disable all warnings
format long g % set display precision
tic; % start timer
end