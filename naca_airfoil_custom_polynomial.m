% =========================================================================
% Custom Polynomial Airfoil Surface Points with Chebyshev Spacing
% Written by Claude Sonnet 4.6, edited by Matthew Brody. Made from multiple scripts
% Compatible with MATLAB and GNU Octave 
%
% Defines upper and lower surfaces via independent polynomial equations of the form:
%
%   y(x) = a0*sqrt(x) + a1*x + a2*x^2 + a3*x^3 + a4*x^4
%   For polynomials multiplied by the chord length and a constant, distribute the constant to
%   all coefficients first, or else your airfoil will be deformed.
%
% where x is the normalized chord coordinate in [0, 1].
%
% Output: a pair of three-column .txt files (0, “y”, -“x”) (<-- code edited, was (x,y) before), 
% tracing the surfaces in the order expected by most CAD curve importers (x, y, z)
%
% Reference for NACA-equivalent coefficients:
%   Abbott & von Doenhoff, "Theory of Wing Sections", 1959
%   NACA thickness:  a0= 0.2969, a1=-0.1260, a2=-0.3516, a3= 0.2843, a4=-0.1015
%   (multiply each by t/0.2 to scale to a desired thickness fraction t)
% =========================================================================

clear; clc; close all;

% =========================================================================
%% CONFIGURATION
% =========================================================================

chord    = 235;    % Chord length [any consistent unit]
n_points = 23;    % Chebyshev points per surface (loop will have 2*n+1 pts after start, end points)

output1 = 'naca_airfoil_top.txt';  % Output CSV filename
output2 = 'naca_airfoil_bottom.txt';
filepath = 'C:/Users/Default/Documents/MATLAB/';  % Please specify your desired directory

% -------------------------------------------------------------------------
% Upper surface polynomial coefficients
%   y_upper(x) = u0*sqrt(x) + u1*x + u2*x^2 + u3*x^3 + u4*x^4
%   x is normalized chord in [0, 1]; output is scaled by chord.
%
% Example below reproduces self-closing NACA 0012 upper surface:
% -------------------------------------------------------------------------
u0 =  0.1773484518;   % coefficient of sqrt(x)
u1 = -0.0756029397;   % coefficient of x
u2 = -0.2128454979;   % coefficient of x^2
u3 =  0.1736375862;   % coefficient of x^3
u4 = -0.0625435473;   % coefficient of x^4

% -------------------------------------------------------------------------
% Lower surface polynomial coefficients
%   y_lower(x) = l0*sqrt(x) + l1*x + l2*x^2 + l3*x^3 + l4*x^4

% -------------------------------------------------------------------------
l0 = -0.1773484518;   % coefficient of sqrt(x)
l1 =  0.0756029397;   % coefficient of x
l2 =  0.2128454979;   % coefficient of x^2  
l3 = -0.1736375862;   % coefficient of x^3
l4 =  0.0625435473;   % coefficient of x^4

% =========================================================================
% CHEBYSHEV SPACING  (normalized chord, ascending from 0 to 1)
% =========================================================================
k = (1 : n_points)';
cos_k  = cos((2*k - 1) * pi / (2 * n_points));
x_norm = sort(0.5 * (1 - cos_k));   % [0, 1], denser near 0 and 1

x = chord * x_norm;   % Physical chord coordinates

% =========================================================================
%% EVALUATE SURFACE POLYNOMIALS
% =========================================================================

xn = x_norm;   % shorthand for normalized coordinate

y_upper = chord .* ( u0*sqrt(xn) + u1*xn + u2*xn.^2 + u3*xn.^3 + u4*xn.^4 );
y_lower = chord .* ( l0*sqrt(xn) + l1*xn + l2*xn.^2 + l3*xn.^3 + l4*xn.^4 );

% =========================================================================
%% ASSEMBLE CONTINUOUS LOOP
%
% CAD importers build a spline through points in order. The loop must be
% non-self-intersecting and close exactly, so the first and last point are
% identical (trailing edge). The leading edge (x=0) is shared between the
% two surfaces — it appears once in the middle of the sequence.
%
% Prepend exact leading edge (0,0) to ensure closure at the LE; similar (0, chord) for TE
x = [0; x; chord];
y_upper = [0; y_upper; 0];
y_lower = [0; y_lower; 0];
x_norm  = [0; x_norm; chord];
%

% =========================================================================
%% WRITE SPACE-DELIMITED .TXT FILES
% =========================================================================

% 61mm +z offset for spar in twist morphing wing simulation, chord reaches in -z.

T1 = [zeros(n_points+2,1), y_upper, 61-x]; 
T2 = [zeros(n_points+2,1), y_lower, 61-x];

csvwrite(fullfile(filepath,output1),T1);
csvwrite(fullfile(filepath,output2),T2);

% =========================================================================
%% PLOT
% =========================================================================

figure('Name', 'Custom Polynomial Airfoil', 'Color', 'white');
hold on; grid on; axis equal;

% Individual surfaces
plot(x, y_upper, 'b-o', 'MarkerSize', 3, 'DisplayName', 'Upper surface');
plot(x, y_lower, 'r-o', 'MarkerSize', 3, 'DisplayName', 'Lower surface');

% Chord and reference points
plot([0 chord], [0 0], 'k:', 'LineWidth', 0.8, 'DisplayName', 'Chord line');
plot(x(end),   y_upper(end), 'ks', 'MarkerSize', 7, 'DisplayName', 'TE (loop start/end)');
plot(x(1),     y_upper(1),   'k^', 'MarkerSize', 7, 'DisplayName', 'LE (loop midpoint)');

% Rug plot of Chebyshev nodes
rug_y = min(y_lower) - 0.04*chord;
plot(x, repmat(rug_y, size(x)), 'k|', 'MarkerSize', 5);
text(chord/2, rug_y - 0.025*chord, 'Chebyshev nodes', ...
    'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.4 0.4 0.4]);

title(sprintf('Custom polynomial airfoil  |  chord = %.3g  |  %d pts/surface  |  %d loop pts', ...
    chord, n_points));
xlabel('x');
ylabel('y');
legend('Location', 'northeast');
xlim([-0.05*chord, 1.1*chord]);
