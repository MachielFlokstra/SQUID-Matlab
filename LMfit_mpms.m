function [yfit,chisq,par,err] = LMfit_mpms(fitfun,x,y,ysig,rot,reps,step,type,par,nfo,tol) %#ok<INUSL>
% INFO
%  This function contains the Levenberg-Marquardt algorithm which tries to reduce the chisquared value of a fit between
%  a set of data points y(x), with standard deviation ysig(x), and a nonlinear function dependent on coefficients a. The
%  Levenberg-Marquardt algorithm needs to be accomanied by a fitfunction that calculates yfit(x,a) and dyda(x,a) (which
%  is the Jacobian, J). The program returns the obtained best fit (+ chisquared, fitparameters, fitparameter errors).
%  This version is adapted to fit MPMS data (fitting multiple scans simulateounsly).
%
% DETAILS OF THE FIT ALGORITHM
%  Each iteration step the algorithm tries to improve the values of the fitparameters by a small amount of da which is
%  calculated by solving the equation: (J'*J + lamda*diag(J'*J))*da = J'*(y-yfit), where ' denotes the transpose and
%  lamda is the damping factor which is tuned during the iteration process. The step da is thus calculated based on the
%  Jacobian and yfit that belong to the current best found values for a. When adding da to a would results in a
%  lowering of chisquared, a is updated. Otherwise, the damping factor is decreased and a new (smaller) da is
%  calculated. The iteration stops when either the calculated step length da, or the reduction of sum of squares from
%  the latest parameter vector (a + da) fall below predefined limits.
%
%  The equation is derived from setting d(chi^2)/d(a) = 0, but it no longer contains explicitly the standard deviations.
%  Once convergence is reached (chisquared is minimized), the errors of the fitparameters can be determined using:
%  asig = sqrt(sum(dy.*dy) * diag(inv(H))/(N-m)), where H = J'*J, N the total number of fitted points, m the number of
%  fit parameters and dy = y-yfit.
%
% INPUT
%  fitfun : name (m-file) of the fitfunction
%  x      : row vector with x data of all scans in series
%  y      : row vector with y data of all scans in series
%  ysig   : row vector with y-error of all scans in series
%  rot    : column vector with rotation per "new" scan // rot values are always in range [0 to 360]
%  reps   : column vector with number of scans per "new" scan
%  step   : column vector with number of steps per "new" scan
%  type   : column vector with the coil-system used per "new" scan (0=longitudinal, 1=transverse)
%  par    : row vector with fitparameter values of all scans in series
%  nfo    : 4-row matrix with for each fitparameter: fit flag / common flag / lower limit / upper limit
%  tol    : tolerance of the fit
%
% OUTPUT
%  yfit   : fit curve for the given set of fitparameters (par)
%  chisq  : reduced chisquared of the fit
%  par    : (optimized) fitparameters (for each scan in series)
%  err    : errors of the fitparameters
% ----------------------------------------------------------------------------------------------------------------------

lamda = 0.01; % set an initial value of the damping factor for the LM
% note: a lower value means a larger value in dpar and leads to faster reaching convergence, but can also lead to ending
%       up in a local minimum with fitparameter values far from their optimum values. 
m = (sum(nfo(1,:))-sum(nfo(2,:)))*sum(reps) + sum(nfo(2,:)); % fitparameters to fit: sum(reps)*uniques + commons
W = 1./(ysig.^2); % diagonal of the weighting matrix

% calculate fitfunction (yfit) and jacobian (J) for initial value of fitparameters (par)
[yfit,J,par] = eval([fitfun '(x,rot,reps,step,type,par,nfo,1)']);
dy = y - yfit; % the difference between fit and data
chisq = sum((dy.*dy)./(ysig.*ysig)); % chisquared of the fit

% calculate covariance matrix H = J' * (diag(W) * J), can't use diag(W) for long W vectors
WJ = zeros(length(W),size(J,2));
for k=1:length(W)
  WJ(k,:) = W(k)*J(k,:);
end
H = J'*WJ;

% iterative fitting routine
chisq_old = chisq; % the chisquared value of the previous iteration
while ((lamda < 1e10) && ( (abs(1-chisq/chisq_old) > tol) || isequal(chisq_old,chisq)))
  chisq_old = chisq; % update value if chisq was improved (need to update after the while() test)
  % calculate new fitparameters using current damping and H matrix by solving: [H + lambda*diag(H)]*dpar = [J'W]*dy
  H_lm = H + lamda*diag(diag(H)); % the new H matrix
  % set the diagonal elements corresponding to the fixed fitparameters to 1
  for k=1:length(par)
    if isequal(H_lm(k,k), 0) % (these elements are zero due zeros in the Jacobian)
      H_lm(k,k) = 1;
    end
  end
  dpar = H_lm \ (J'*(W.*dy)'); % the change needed in a to optimize
  par_lm = par + dpar'; % the optimized fitparameter values
  par_lm = applylim(par_lm,nfo); %#ok<NASGU> % keep fitparameters within limits
  
  % calculate goodness-of-fit (chisquared) when using the new fitparameters
  [yfit_lm,~,par_lm] = eval([fitfun '(x,rot,reps,step,type,par_lm,nfo,0)']); % no need to calculate J here
  dy_lm = y - yfit_lm;
  chisq_lm = sum((dy_lm.*dy_lm)./(ysig.*ysig));
  
  % did the fit improve? (if so, keep the updated fitparameters)
  if (chisq_lm < chisq)
    lamda = lamda/10; % reduce damping for next iteration round
    % keep the obtained fitparameters, yfit, dy and chisq
    yfit = yfit_lm;
    par = par_lm;
    dy = dy_lm;
    chisq = chisq_lm;

    % calculate Jaccobian of the new fitparameters
    [~,J,~] = eval([fitfun '(x,rot,reps,step,type,par,nfo,1)']);
    WJ = zeros(length(W),size(J,2));
    for k=1:length(W)
      WJ(k,:) = W(k)*J(k,:);
    end
    H = J'*WJ;
  else % fit did not improve: increases the damping factor
    lamda = lamda*10;
  end
end
          
% determine the errors in the optimized parameter
dy = y - yfit;
wi2 = sum(dy.*dy)/(length(y)-m+1); % all the w^2(i) are now equal to the mean square measurement error
% note: W thus becomes the identity matrix multiplied by 1/wi2
H = (J'*J)/wi2;
for k=1:length(par)
  if isequal(H(k,k), 0)
    H(k,k) = 1;
  end
end
err = sqrt(diag(inv(H)))';

% get errors for common fitparameters from the first scan
for k=1:length(par)
  if isequal(nfo(2,1+rem(k-1,size(nfo,2))),1)
    err(k) = err(1+rem(k-1,size(nfo,2)));
  end
end

% set errors to zero for fixed fitparameters
for k=1:length(par)
  if isequal(nfo(1,1+rem(k-1,size(nfo,2))),0)
    err(k) = 0;
  end
end

% calculate the reduced chisquared of the fit
chisq = sum((dy.*dy)./(ysig.*ysig));
chisq = chisq/(length(y)-m); % reduced chi^2
end

% ======================================================================================================================
function parout = applylim(parin,nfo)
% keeps parameters within limits
parout = parin;
for k=1:length(parin)
  ai = 1+rem(k-1,size(nfo,2));
  % only deal with parameter that are kept free
  if isequal(nfo(1,ai),1)
    % for phi0 subtract/add integers of 360 deg to stay in [0-360] range
    if isequal(ai,2)
      if parin(k) > 360
        parout(k) = rem(parin(k),360);
      end
      if parin(k) < 0
        parout(k) = 360 + rem(parin(k),360);
      end
    else
      % check lower limit
      if parin(k) < nfo(3,ai)
        parout(k) = nfo(3,ai);
      end
      % check upper limit
      if parin(k) > nfo(4,ai)
        parout(k) = nfo(4,ai);
      end
    end
  end
end
end