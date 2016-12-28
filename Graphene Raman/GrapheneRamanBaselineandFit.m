% This is a program to baseline and fit Raman spectra of graphene
% (single- and multi-layer). Currently suitable with D, G, and 2D peaks.
%
% 
% Instructions:
% 1. Place 'GrapheneRamanBaselineandFit.m' and your data file into 
%       your current Matlab directory.
% 2. Import your data file as 2 separate 1xn doubles (1 for x values and 1
%       for y values) called 'xdata' and 'ydata'.
% 3. Run "GrapheneRamanBaselineandFit(xdata,ydata)" from the Matlab 
%       command window (for more advanced options, see below).
% 4. The command window will output the intensities (a.u.), location
%       (cm-1), and half-widths at half-maximum (HWHM) for each peak in
%       addition to a plot and data of the baselined spectrum.
%
%
% Baselining is performed via asymmetric reweighted penalized least 
% square (arPLS) baseline removal, which was reported by Baek et al.
% in their 2015 paper entitled "Baseline correction using asymmetrically 
% reweighted penalized least squares smoothing." Part of the baseline code
% is also taken from CRIkit by Charles H. Camp Jr (charles.camp@nist.gov, 
% ccampjr@gmail.com) and can be found at https://github.com/CoherentRaman
% NIST/CRIkit. If you do use this code, please do cite them.
% 
%
% Peak fitting is performed using Matlab's built-in nonlinear least-
% squares solver where each Raman mode is a baselined Gaussian curve.
%
% Input:    xdata (1xn double containing Raman spectrum x-values in cm-1)
%           ydata (1xn double containing Raman spectrum y-values in any units)
%           smoothness_param (*optional*) (default is 1e3)
%           min_diff (*optional*) (default is 1e-6)
%
% Output:   Baselined data (nx2 double labeled as 'parseddata')
%           Peak (D, G, 2D) intensities
%           Peak locations
%           Peak half-widths at half-maximum (HWHM)


function [parseddata,xD,xG,x2D] = GrapheneRamanBaselineandFit(xdata,ydata,smoothness_param,min_diff)

ORDER = 2; % Difference filter order
MAX_ITER = 100; % Maximum iterations

[m,n] = size(ydata);
if (m ~= 1 && n ~= 1)
    error('This function only accepts 1D (effective) signals');
end
ydata = ydata(:);

if nargin == 2
    smoothness_param = 1e3;
    min_diff = 1e-6;
end

% Initialize baselining parameters
signal_length = length(ydata);
difference_matrix = diff(speye(signal_length), ORDER);
minimization_matrix = (smoothness_param*difference_matrix')*difference_matrix;
penalty_vector = ones(signal_length,1);

% Execute baselining
for count = 1:MAX_ITER
    penalty_matrix = spdiags(penalty_vector, 0, signal_length, signal_length);
    % Cholesky decomposition
    C = chol(penalty_matrix + minimization_matrix);
    baseline = C \ (C'\(penalty_vector.*ydata));
    d = ydata - baseline;
    % make d-, and get penalty_vector^t with m and s
    dn = d(d<0);
    m = mean(dn);
    s = std(dn);
    penalty_vector_temp = 1./(1+exp(2*(d-(2*s-m))/s));
    % check exit condition and backup
    if norm(penalty_vector-penalty_vector_temp)/norm(penalty_vector) < min_diff
%         count
        break;
    end
    penalty_vector = penalty_vector_temp;
end

% Parse baselined data
baselineddata = [transpose(xdata),ydata - baseline];
output = baselineddata;
[m,n] = size(baselineddata);
parseddata = [];
for i = 1:m
    if baselineddata(i,1) > 600
        parseddata = [parseddata;baselineddata(i,1:2)];
    else
        break
    end
end
[o,p] = size(parseddata);
xvals = parseddata(1:o,1);
yvals = parseddata(1:o,2);

% Plot baselined data
plot(xvals,yvals);
title('Please click on the peaks (D, G, and 2D) from left to right')
xlabel('Raman Shift (cm-1)')
ylabel('Intensity (a.u.)')

% Use user input as initial guesses (i.e. x0 matrix) for lsqcurvefit
[xguess,yguess] = ginput(3);

% Initialize fitting parameters
x0 = [yguess(1) yguess(2) yguess(3) xguess(1) xguess(2) xguess(3) 5 5 5];
fun1 = @(xi,xvals)(xi(1)./(1+(((xvals-xi(2))./xi(3)).^2)));
fun3 = @(x,xvals)(x(1)./(1+(((xvals-x(4))./x(7)).^2)))+(x(2)./(1+(((xvals-x(5))./x(8)).^2)))+(x(3)./(1+(((xvals-x(6))./x(9)).^2)));

% Perform fitting and isolate values for D, G, and 2D peaks
x = lsqcurvefit(fun3,x0,xvals,yvals);
xD = [x(1) x(4) x(7)];
xG = [x(2) x(5) x(8)];
x2D = [x(3) x(6) x(9)];

% Plot the baselined data and three corresponding fits
plot(xvals,yvals,xvals,fun1(xD,xvals),xvals,fun1(xG,xvals),xvals,fun1(x2D,xvals))
legend('Baselined Data', 'D fit', 'G fit', '2D fit','Location','northwest')
title('Baselined Data and Final Fits')
xlabel('Raman Shift (cm-1)')
ylabel('Intensity (a.u.)')

D_intensity = xD(1)
D_location = xD(2)
D_HWHM = xD(3)

G_intensity = xG(1)
G_location = xG(2)
G_HWHM = xG(3)

intensity_2D = x2D(1)
location_2D = x2D(2)
HWHM_2D = x2D(3)