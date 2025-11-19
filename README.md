[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/logistic4)

📘 Overview
logistic4 is a MATLAB mini-toolbox implementing the Four-Parameter Logistic (4PL) regression model, widely used in immunoassays and dose-response analyses such as ELISA, RIA, IRMA, and general bioassay standard curves. The 4PL curve describes the characteristic sigmoidal (“S-shaped”) relationship between analyte concentration and measured signal while modelling both the lower asymptote (A) and upper asymptote (D).

The mathematical form of the 4PL equation is:
    F(x) = D + (A - D) / (1 + (x / C)^B)

where:
- A : Minimum asymptote (response at zero concentration)
- B : Hill slope (steepness and direction)
- C : Inflection point (EC50)
- D : Maximum asymptote (response at infinite concentration)

The repository provides two main functions:
- L4P    : Fits a 4PL model to experimental data
- L4Pinv : Computes the inverse of the 4PL model, estimating x from a given response y

✨ Features
- Robust 4PL nonlinear regression with support for replicate data
- Parameter starting points automatically estimated when not provided
- Lower and upper bounds for parameters inferred from data (or user-defined)
- Weighted regression when replicates are provided (row-wise standard deviations)
- Returns a cfit model object and full goodness-of-fit statistics
- Forward evaluation (model prediction) and inverse evaluation (interpolation)
- Suitable for bioassay calibration curves, dose-response studies, and QC applications

📥 Installation
1. Download the logistic4 repository:
   https://github.com/dnafinder/logistic4

2. Place the files (L4P.m and L4Pinv.m) in any directory.

3. Add the folder to your MATLAB path:
      addpath('path_to_logistic4')

4. Verify that MATLAB can find the functions:
      which L4P
      which L4Pinv

⚙️ Requirements
- MATLAB (any recent version)
- Curve Fitting Toolbox (for fit, fittype, cfit objects)

📈 Usage Examples
Fitting a 4PL curve:

    x = [0; 4.5; 10.6; 19.7; 40; 84; 210];
    y = [0.0089; 0.0419; 0.0873; 0.2599; 0.7074; 1.528; 2.7739];

    [cf, G] = L4P(x, y);

Plotting the fitted curve:

    plot(x, y, 'ro');
    hold on;
    plot(cf, 'r');
    hold off;

Interpolating unknown samples:

    yq = 1.8;
    x_est = L4Pinv(cf, yq);

Using a numeric parameter vector:

    params = [A, B, C, D];
    x_est = L4Pinv(params, yq);

🔢 Inputs
L4P(x, y)
L4P(x, y, st)
L4P(x, y, st, L, U)

- x : Column vector of concentrations (N×1)
- y : Column vector of responses OR matrix of replicates (N×M)
- st: Optional starting values [A0 B0 C0 D0]
- L : Optional lower bounds [Amin Bmin Cmin Dmin]
- U : Optional upper bounds [Amax Bmax Cmax Dmax]

L4Pinv(cf, y)
- cf : cfit object from L4P OR numeric vector [A B C D]
- y  : Query response values (any size)

📤 Outputs
L4P returns:
- cf : cfit object representing the fitted 4PL curve
- G  : Goodness-of-fit structure with fields:
       sse, rsquare, adjrsquare, dfe, rmse

L4Pinv returns:
- x : Array of interpolated x-values matching the size of y

🧠 Interpretation
- A and D represent the lower and upper plateaus of the curve
- B controls steepness and direction (positive = increasing curve)
- C represents the EC50 or midpoint of the transition zone
- A good fit typically shows:
  • R² and adjusted R² close to 1  
  • Low SSE and RMSE  
  • A smooth curve matching the sigmoidal trend of the data  

📌 Notes
- If y is provided as an N×M matrix, L4P automatically averages replicates and uses row-wise standard deviations as weights.
- Nonlinear regression may require reasonable starting values and bounds for best convergence.
- Response values equal to or outside the asymptotes A and D may lead to extrapolation when inverted with L4Pinv.

🧾 Citation
If you use logistic4 in publications or analytical workflows, please cite:

Cardillo G. (2025). logistic4: Four-parameter logistic regression tools in MATLAB (L4P and L4Pinv).  
Available at: https://github.com/dnafinder/logistic4

👤 Author
Giuseppe Cardillo  
Email: giuseppe.cardillo.75@gmail.com  
GitHub: https://github.com/dnafinder

📄 License
logistic4 is distributed under the terms specified in the LICENSE file available at:
https://github.com/dnafinder/logistic4
