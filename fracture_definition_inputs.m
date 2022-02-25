function [L, N, dx, H, roughness, mismatch, aniso, seed, lambda_0, model] = fracture_definition_inputs()
% This function generates all the necesary inputs for rough fracture
% generator
L         = 25.4*3;     % Fracture length                                                     [mm]
N         = 2^10;       % Number of divisions along fracture side                             [power of 2]
H         = 0.8;        % Hurst exponent. Determines fractal dimension D=3-H                  [0.45<H<0.85]
roughness = 0.005*L;    % root-mean-square roughness (~ standard deviation of heights)        [same units as L]
mismatch  = 0.05;       % Mismatch length scale (wavelength) as a fraction of fracture size   [0<Mismatch<1]
aniso     = 1.0;        % anisotropy ratio                                                    [0<Aniso<1]
seed      = 2;          % Random number seed (phase generation)                               [+ve integer]
dx        = L/N;        % Grid spacing                                                        [same units as L]
lambda_0  = 0.25;       % Roll-ff wavelenth                                                   [0<lambda_0<1]
model     = 'bilinear'; % Power spectal density model                                         ['linear', 'linear', 'bilinear', 'smooth']
end

