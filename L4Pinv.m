function x = L4Pinv(cf, y)
%L4PINV Inverse of the four-parameter logistic (4PL) equation.
%
%   x = L4Pinv(cf, y)
%
%   Description:
%       L4Pinv computes the inverse of the four-parameter logistic (4PL)
%       model used in L4P. Given a fitted 4PL curve (or a numeric parameter
%       vector) and a set of response values y, this function returns the
%       corresponding x values that satisfy:
%
%           y = D + (A - D) / (1 + (x / C)^B)
%
%       where:
%           A : minimum asymptote (response at zero concentration)
%           B : Hill slope (steepness and direction of curve)
%           C : inflection point (EC50)
%           D : maximum asymptote (response at infinite concentration)
%
%   Syntax:
%       x = L4Pinv(cf, y)
%
%   Inputs:
%       cf  - Either:
%              • A cfit object returned by L4P, containing parameters
%                A, B, C, and D of the 4PL model, or
%              • A numeric vector [A, B, C, D] (1x4 or 4x1).
%
%       y   - Numeric array of response values for which you want the
%             corresponding x values. y must be real, finite, non-empty.
%             It may be a scalar, vector, or matrix; x will match its size.
%
%   Outputs:
%       x   - Numeric array of x values such that the 4PL model evaluated
%             at x returns (approximately) the values in y. The size of x
%             is the same as the size of y.
%
%   Example:
%       % Fit a 4PL curve using L4P:
%       xdata = [0; 4.5; 10.6; 19.7; 40; 84; 210];
%       ydata = [0.0089; 0.0419; 0.0873; 0.2599; 0.7074; 1.528; 2.7739];
%
%       [cf, G] = L4P(xdata, ydata);
%
%       % Invert the fitted curve at a specific response:
%       yq = 1.782315;
%       xq = L4Pinv(cf, yq);
%
%       % Alternatively, using the explicit parameter vector:
%       params = coeffvalues(cf);      % [A, B, C, D]
%       xq2   = L4Pinv(params, yq);    % should match xq
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
%       This function is distributed under the terms specified in the
%       LICENSE file of the logistic4 repository.
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
%       L4P

% -----------------------------
% Input parsing and validation
% -----------------------------

p = inputParser;

% cf must be either a fit-like object or a numeric vector of length 4.
addRequired(p, 'cf', @(v) ...
    isobject(v) || ...
    (isnumeric(v) && isvector(v) && numel(v) == 4 && ...
     all(isfinite(v(:)) & isreal(v(:)) & ~isnan(v(:)))));

% y must be a numeric, real, finite, non-empty array (any shape allowed).
addRequired(p, 'y', @(v) ...
    validateattributes(v, {'numeric'}, ...
                       {'real', 'finite', 'nonnan', 'nonempty'}));

parse(p, cf, y);
cf = p.Results.cf;
y  = p.Results.y;

clear p;

% -----------------------------
% Extract 4PL parameters (A, B, C, D)
% -----------------------------
% If cf is a fit object, extract its coefficients; otherwise use the
% numeric vector directly.

if isobject(cf)
    params = coeffvalues(cf); % expects [A, B, C, D]
else
    % Ensure params is a row vector [A, B, C, D].
    params = cf(:).';
end

A = params(1);
B = params(2);
C = params(3);
D = params(4);

% -----------------------------
% Basic sanity checks on y
% -----------------------------
% For a standard 4PL curve, the response typically lies between A and D.
% Values outside the open interval (min(A, D), max(A, D)) may still be
% computed, but can lead to extrapolation or non-physical x values.
yMinModel = min(A, D);
yMaxModel = max(A, D);

if any(y(:) <= yMinModel)
    warning('L4Pinv:BelowModelRange', ...
        ['Some response values are <= the lower asymptote. ' ...
         'Inversion may be unreliable or extrapolative.']);
end

if any(y(:) >= yMaxModel)
    warning('L4Pinv:AboveModelRange', ...
        ['Some response values are >= the upper asymptote. ' ...
         'Inversion may be unreliable or extrapolative.']);
end

% -----------------------------
% Invert the 4PL equation
% -----------------------------
% Starting from:
%   y = D + (A - D) / (1 + (x / C)^B)
%
% Rearranging:
%   y - D = (A - D) / (1 + (x / C)^B)
%   (A - D) / (y - D) = 1 + (x / C)^B
%   (x / C)^B = (A - D) / (y - D) - 1
%   x / C     = ((A - D) / (y - D) - 1)^(1 / B)
%   x         = C * ((A - D) / (y - D) - 1)^(1 / B)
%
% The computation is performed element-wise over y.

ratio = ((A - D) ./ (y - D)) - 1;
x     = C .* (ratio .^ (1 ./ B));

end
