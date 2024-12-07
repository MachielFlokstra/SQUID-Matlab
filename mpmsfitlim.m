function lim = mpmsfitlim()
% INFO
%  Sets the default fitparameter limits for r, phi0, z0, Mr, Mp, Mz, c0 and c1
% ----------------------------------------------------------------------------------------------------------------------
lim = zeros(2,8); % min/max for the 8 fitparameters
lim(:,1) = [0 10]*1E-3; % r (m) // don't allow negative r
lim(:,2) = [0 360]; % phi0 (deg) // keep phi0 within [0 to +360] by taking off integers of 360
lim(:,3) = [-10 10]*1E-2; % z0 (m)
lim(:,4) = [-1 1]*1E10; % Mr (Am^2)
lim(:,5) = [-1 1]*1E10; % Mp (Am^2)
lim(:,6) = [-1 1]*1E10; % Mz (Am^2)
lim(:,7) = [-1 1]*1E10; % c0 (V) // the background offset
lim(:,8) = [-1 1]*1E10; % c1 (V/m) // the background slope
end