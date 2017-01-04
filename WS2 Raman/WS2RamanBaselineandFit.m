% This is a program to baseline and fit Raman spectra of two-dimensional
% tungsten disulfide (WS2). Currently suitable with the following peaks:
%
% LA(M) ~ 175 cm^(-1)
% A1g-LA(M) ~ 230 cm^(-1)
% 2LA(M)-3E2g2 ~ 265 cm^(-1)
% 2LA(M)-2E2g2 ~ 296 cm^(-1)
% 2LA(M)-E2g2 ~ 323 cm^(-1)
% 2LA(M) ~ 350 cm^(-1)
% E2g1 ~ 355 cm^(-1)
% A1g ~ 418 cm^(-1)
%
% 
% Instructions:
% 1. Place 'WS2RamanBaselineandFit.m' and your data file into your current
%       Matlab directory.
% 2. Import your data file as 2 separate 1xn doubles (1 for x values and 1
%       for y values) called xdata and ydata.
% 3. Run "WS2RamanBaselineandFit(xdata,ydata)" from the Matlab command
%       window (for more advanced options, see below).
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
%           Peak (LA(M), A1g-LA(M), 2LA(M)-3E2g2, 2LA(M)-2E2g2, 2LA(M)-E2g2, 
%                   2LA(M), E2g1, A1g) intensities
%           Peak locations
%           Peak half-widths at half-maximum (HWHM)


function [parseddata,xLAM,xA1g_LAM,x2LAM_3E2g2,x2LAM_2E2g2,x2LAM_E2g2,x2LAM,xE2g1,xA1g] = WS2RamanBaselineandFit(xdata,ydata,smoothness_param,min_diff)

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
    if baselineddata(i,1) < 500
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
title('Please click on the peaks (LA(M),A1g-LA(M),2LA(M)-3E2g2,2LA(M)-2E2g2,2LA(M)-E2g2,2LA(M),E2g1,A1g) from left to right')
xlabel('Raman Shift (cm-1)')
ylabel('Intensity (a.u.)')

% Use user input as initial guesses (i.e. x0 matrix) for lsqcurvefit
[xguess,yguess] = ginput(8);

% Initialize fitting parameters
x0 = [yguess(1) yguess(2) yguess(3) yguess(4) yguess(5) yguess(6) yguess(7) yguess(8) xguess(1) xguess(2) xguess(3) xguess(4) xguess(5) xguess(6) xguess(7) xguess(8) 5 5 5 5 5 5 5 5];
fun1 = @(xi,xvals)(xi(1)./(1+(((xvals-xi(2))./xi(3)).^2)));
fun3 = @(x,xvals)(x(1)./(1+(((xvals-x(9))./x(17)).^2)))+(x(2)./(1+(((xvals-x(10))./x(18)).^2)))+(x(3)./(1+(((xvals-x(11))./x(19)).^2)))+(x(4)./(1+(((xvals-x(12))./x(20)).^2)))+(x(5)./(1+(((xvals-x(13))./x(21)).^2)))+(x(6)./(1+(((xvals-x(14))./x(22)).^2)))+(x(7)./(1+(((xvals-x(15))./x(23)).^2)))+(x(8)./(1+(((xvals-x(16))./x(24)).^2)));

% Perform fitting and isolate values for LA(M), A1g-LA(M), 2LA(M)-3E2g2, 2LA(M)-2E2g2, 2LA(M)-E2g2, 2LA(M), E2g1, A1g peaks
x = lsqcurvefit(fun3,x0,xvals,yvals);
xLAM = [x(1) x(9) x(17)];
xA1g_LAM = [x(2) x(10) x(18)];
x2LAM_3E2g2 = [x(3) x(11) x(19)];
x2LAM_2E2g2 = [x(4) x(12) x(20)];
x2LAM_E2g2 = [x(5) x(13) x(21)];
x2LAM = [x(6) x(14) x(22)];
xE2g1 = [x(7) x(15) x(23)];
xA1g = [x(8) x(16) x(24)];

% Plot the baselined data and three corresponding fits
plot(xvals,yvals,xvals,fun1(xLAM,xvals),xvals,fun1(xA1g_LAM,xvals),xvals,fun1(x2LAM_3E2g2,xvals),xvals,fun1(x2LAM_2E2g2,xvals),xvals,fun1(x2LAM_E2g2,xvals),xvals,fun1(x2LAM,xvals),xvals,fun1(xE2g1,xvals),xvals,fun1(xA1g,xvals))
legend('Baselined Data', 'LA(M) fit', 'A1g-LA(M) fit', '2LA(M)-3E2g2 fit','2LA(M)-2E2g2 fit', '2LA(M)-E2g2 fit', '2LA(M) fit','E2g1 fit', 'A1g fit', 'Location','northwest')
title('Baselined Data and Final Fits')
xlabel('Raman Shift (cm-1)')
ylabel('Intensity (a.u.)')

LAM_intensity = xLAM(1)
LAM_location = xLAM(2)
LAM_HWHM = xLAM(3)

A1g_LAM_intensity = xA1g_LAM(1)
A1g_LAM_location = xA1g_LAM(2)
A1g_LAM_HWHM = xA1g_LAM(3)

intensity_2LAM_3E2g2 = x2LAM_3E2g2(1)
location_2LAM_3E2g2 = x2LAM_3E2g2(2)
HWHM_2LAM_3E2g2 = x2LAM_3E2g2(3)

intensity_2LAM_2E2g2y = x2LAM_2E2g2(1)
location_2LAM_2E2g2 = x2LAM_2E2g2(2)
HWHM_2LAM_2E2g2 = x2LAM_2E2g2(3)

intensity_2LAM_E2g2 = x2LAM_E2g2(1)
location_2LAM_E2g2 = x2LAM_E2g2(2)
HWHM_2LAM_E2g2 = x2LAM_E2g2(3)

intensity_2LAM = x2LAM(1)
location_2LAM = x2LAM(2)
HWHM_2LAM = x2LAM(3)

E2g1_intensity = xE2g1(1)
E2g1_location = xE2g1(2)
E2g1_HWHM = xE2g1(3)

A1g_intensity = xA1g(1)
A1g_location = xA1g(2)
A1g_HWHM = xA1g(3)