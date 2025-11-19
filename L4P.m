function [cf, G] = L4P(x, y, varargin)
%L4P Four-parameter logistic regression (4PL).
%
%   [cf, G] = L4P(x, y)
%   [cf, G] = L4P(x, y, st)
%   [cf, G] = L4P(x, y, st, L, U)
%
%   Description:
%       L4P fits a four-parameter logistic (4PL) model to experimental
%       data, commonly used in immunoassays and dose-response curves
%       (e.g., ELISA, RIA, IRMA). The 4PL model captures the sigmoidal
%       shape of standard curves, including minimum asymptote (A), maximum
%       asymptote (D), slope (B), and inflection point (C).
%
%       The 4PL equation is:
%
%           F(x) = D + (A - D) / (1 + (x / C)^B)
%
%       where:
%           A : Minimum asymptote (response at zero concentration)
%           B : Hill slope (steepness and direction of curve)
%           C : Inflection point (EC50)
%           D : Maximum asymptote (response at infinite concentration)
%
%   Inputs:
%       x   - Column vector (N x 1) of x-values (e.g., concentrations).
%       y   - Column vector (N x 1) of responses, or matrix (N x M) of
%             replicates. When y is a matrix, row means are used as
%             responses and row-wise standard deviations as weights.
%
%       st  - (Optional) Starting values [A0, B0, C0, D0].
%             If empty or not provided, starting points are estimated.
%
%       L   - (Optional) Lower bounds for parameters [Amin, Bmin, Cmin, Dmin].
%             If empty or not provided, bounds are inferred from data.
%
%       U   - (Optional) Upper bounds for parameters [Amax, Bmax, Cmax, Dmax].
%             If empty or not provided, bounds are inferred from data.
%
%   Outputs:
%       cf  - A cfit object representing the fitted 4PL curve.
%       G   - Structure with goodness-of-fit measures:
%                 G.sse, G.rsquare, G.adjrsquare, G.rmse, G.dfe
%
%   Example:
%       x = [0; 4.5; 10.6; 19.7; 40; 84; 210];
%       y = [0.0089; 0.0419; 0.0873; 0.2599; 0.7074; 1.528; 2.7739];
%
%       [cf, G] = L4P(x, y);
%
%       plot(x, y, 'ro'); hold on; plot(cf, 'r'); hold off;
%
%   GitHub repository:
%       https://github.com/dnafinder/logistic4
%
%   Citation:
%       Cardillo G. (2025). logistic4: Four-parameter logistic regression
%       tools in MATLAB (L4P and L4Pinv). Available at:
%       https://github.com/dnafinder/logistic4
%
%   License:
%       Distributed under the terms specified in the LICENSE file
%       in the logistic4 repository.
%
%   Author:
%       Giuseppe Cardillo
%       giuseppe.cardillo.75@gmail.com
%
%   Created:
%       2012-01-01 (original concept)
%
%   Updated:
%       2025-11-19 (refactored and documented version)
%
%   Version:
%       1.1.0
%
%   See also:
%       L4Pinv, L3P, L3Pinv

% -----------------------------
% Input parsing and validation
% -----------------------------
p = inputParser;

addRequired(p, 'x', @(v) validateattributes(v, ...
    {'numeric'}, {'column', 'real', 'finite', 'nonnan', 'nonempty'}));

addRequired(p, 'y', @(v) validateattributes(v, ...
    {'numeric'}, {'2d', 'real', 'finite', 'nonnan', 'nonempty'}));

addOptional(p, 'st', [], @(v) isempty(v) || ...
    validateattributes(v, {'numeric'}, ...
        {'row', 'real', 'finite', 'nonnan', 'numel', 4}));

addOptional(p, 'L', [], @(v) isempty(v) || ...
    validateattributes(v, {'numeric'}, ...
        {'row', 'real', 'nonnan', 'numel', 4}));

addOptional(p, 'U', [], @(v) isempty(v) || ...
    validateattributes(v, {'numeric'}, ...
        {'row', 'real', 'nonnan', 'numel', 4}));

parse(p, x, y, varargin{:});

x  = p.Results.x;
y  = p.Results.y;
st = p.Results.st;
L  = p.Results.L;
U  = p.Results.U;

clear p;

% Ensure x and y have same number of rows
assert(size(x,1) == size(y,1), ...
    'L4P:SizeMismatch', 'x and y must have the same number of rows.');
assert(size(x,1) >= 4, ...
    'L4P:NotEnoughData', 'At least 4 data points are required.');

% -----------------------------
% Compute means & weights if replicates are provided
% -----------------------------
if ~isvector(y)
    we = std(y, 0, 2);
    y  = mean(y, 2);
else
    y  = y(:);
    we = zeros(size(x));
end

% -----------------------------
% Starting values estimation
% -----------------------------
slope = (y(end) - y(1)) / (x(end) - x(1));

if isempty(st)
    yRange = max(y) - min(y);
    target = yRange / 2;
    [~, Idx] = min(abs(y - target));
    
    st = [
        min(y)           ... % A0
        sign(slope)      ... % B0
        x(Idx)           ... % C0
        max(y)             ... % D0
    ];
end

% -----------------------------
% Parameter bounds
% -----------------------------
if isempty(L)
    L = zeros(1,4);
    if slope < 0
        L(2) = -Inf;
    end
end

if isempty(U)
    U = Inf(1,4);
    if slope < 0
        U(2) = 0;
    end
end

% -----------------------------
% Define model and fit options
% -----------------------------
fo = fitoptions('Method', 'NonlinearLeastSquares', ...
    'StartPoint', st, ...
    'Lower', L, ...
    'Upper', U);

if all(we)
    set(fo, 'Weights', we);
end

ft = fittype('D + (A - D) / (1 + (x / C)^B)', ...
    'dependent',  {'y'}, ...
    'independent',{'x'}, ...
    'coefficients',{'A', 'B', 'C', 'D'});

% -----------------------------
% Fit the 4PL model
% -----------------------------
[cf, G] = fit(x, y, ft, fo);

end
