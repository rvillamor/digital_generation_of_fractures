% Rafael Villamor Lora
% February 11th, 2020
% GENERATE 3D MODELS (STL) OF A ROUGH FRACTURE
%
% This script shows, step-by-step, how to generate the STL files of a rough
% fractures using teh RSG_brown1995() and fields2bodies() functions

clc, close all, clear all

%% --- FIRST STEP: GEOMETRY GENERATION
% One can define the 3D geometry of a rough fracture using different 
% methods, for instance:
% (1) reconstruction from real rocks using 3D scanning 
% (2) generation using traditional Computer-Aided Design software
% (3) from  mathematical models <-- USED HERE

% Generate 3D model of a fracture using Brown's (1995) model:
% Brown, S. R. (1995). Simple Mathematical-Model of a Rough Fracture. 
% Journal of Geophysical Research-Solid Earth, 100(B4), 5941–5952.
% https://doi.org/10.1029/94JB03262
% Model inputs
%{
L         : Fracture length                                                     [mm]
N         : Number of divisions along fracture side                             [power of 2]
dx        : Grid spacing                                                        [same units as L]
H         : Hurst exponent. Determines fractal dimension D=3-H                  [1<H<1]
roughness : root-mean-square roughness (~ standard deviation of heights)        [same units as L]
mismatch  : Mismatch length scale (wavelength) as a fraction of fracture size   [0<Mismatch<1]
aniso     : anisotropy ratio                                                    [0<Aniso<1]
seed      : Random number seed (phase generation)                               [+ve integer]
lambda_0  : Roll-ff wavelenth                                                   [0<lambda_0<1]
model     : Power spectal density model                                         ['linear', 'linear', 'bilinear', 'smooth']

%}
[L, N, dx, H, roughness, mismatch, aniso, seed, lambda_0, model] = fracture_definition_inputs();
% Model outputs
%{
upperSurface  : Upper surface of the fracture
lowerSurface  : Lower surface of the fracture 
apertureField : Aerture field map
%}
[upperSurface, lowerSurface, apertureField] = RSG_brown1995(H, roughness, mismatch, N, aniso, seed, lambda_0, model);

% Translate the fracture surfaces until they are in contact
lowerSurface  = lowerSurface - min(lowerSurface(:));
shift_upper   = upperSurface - lowerSurface;
upperSurface  = upperSurface - min(shift_upper(:));
% Translate the fractyre surfaces so the aperture volume is located at 
% z = 0 (origin of heights)
shift_both   = (mean(lowerSurface(:)) + mean(upperSurface(:))) / 2;
lowerSurface = lowerSurface - shift_both;
upperSurface = upperSurface - shift_both;

% Composite apeture specimen: since larger aperture = larger gaps = lower
% topography level --> the composite aperture is the opposite of the
% apertureField. Note that now, zero apertures (i.e., contacting points) 
% are located at z = 0
compositeAperture = - apertureField;

%% --- SECOND STEP: GENERATE 3D BODIES FROM 2D FIELDS
% Now, we will convert these 2D matrices which define the topography of the
% fracture surfaces into 3D bodies, and wil; output the results as STL 
% files

% Specify the settings for the generation of the 3D bodies
%{
geometry   : Specimen geometry --> 'prism' or 'cylinder'
outputFile : Output file address
bodyDim    : Body dimensions as a fraction of the length [length, width, height OR diameter, height]
%}
settings_lowerSurface                 = settings_3Dbodies();
settings_lowerSurface.outputFile      = "Results/lowerSurface.stl";
settings_upperSurface                 = settings_3Dbodies();
settings_upperSurface.outputFile      = "Results/upperSurface.stl";
settings_compositeAperture            = settings_3Dbodies();
settings_compositeAperture.outputFile = "Results/compositeAperture.stl";
% Generate the bodies and STL files
%{
FV : triangulated patch defined by a structure with fields 'vertices' and 
     'faces'
%}
[FV_lowerSurface]      = fields2bodies(lowerSurface, dx, settings_lowerSurface, "lowerSurface");
[FV_upperSurface]      = fields2bodies(upperSurface, dx, settings_upperSurface, "upperSurface");
[FV_compositeAperture] = fields2bodies(compositeAperture, dx, settings_compositeAperture); 

%% --- THIRD STEP (optional): PLOT RESULTS

% DRAW UPPER AND LOWER FARCTURE BODIES
figure
% Lower volume
FV = FV_lowerSurface;
trisurf(FV.faces, FV.vertices(:,1), FV.vertices(:,2), FV.vertices(:,3)); hold on  
% Upper volume
FV        = FV_upperSurface;
rot_mat   = [-1 0 0; 0 1 0; 0 0 -1];                                              % Rotation matrix arounf y-axis (180 degrees)
R         = rot_mat * [FV.vertices(:,1) FV.vertices(:,2) FV.vertices(:,3)]';      % Flip the upper body so it's aligned with the lower one
topHeight = 2 * max(FV.vertices(:,3));
trisurf(FV.faces, R(1,:)', R(2,:)', R(3,:)' + topHeight);                
axis equal vis3d, shading interp, camlight left, lighting gouraud, colormap white % Make everything look pretty 

% DRAW THE COMPOSITE APERTURE SPECIMEN
figure
FV = FV_compositeAperture;
trisurf(FV.faces, FV.vertices(:,1), FV.vertices(:,2), FV.vertices(:,3));            % Plot aperture volume
axis equal vis3d, shading interp, camlight left, lighting gouraud, colormap white   % Make everything look pretty 
