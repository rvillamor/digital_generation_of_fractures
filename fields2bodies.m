function [FV_fracture] = fields2bodies(surface, dx, settings, position)
% GENERATE 3D MODELS FROM 2D FIELDS
%
% [outputArg1,outputArg2] = fields2bodies(upperSurface, lowerSurface, dx, settings)
%
% This function takes 2D matrices that define the z-coordinates of the
% farcture surface, generates 3D bodies, and output the results as a 
% STL file
%
% Inputs
% -------------------------------------------------------------------------
% surface  : 2D matrix defining the topograpgy of the field   [length]
% dx       : Pixel size                                       [length]
% settings : Structure with fields 
%              'geometry'   - Specimen geometry --> 'prism' or 'cylinder'
%              'bodyDim'    - Body dimensions as a fraction of the length 
%                             [length, width, height OR diameter, height]
%              'outputFile' - Output file address 
% position : Position of the surface relative to the fracture plane -->
%            "upperSurface" or "LowerSurface"
%
% Outputs
% -------------------------------------------------------------------------
% FV       : triangulated patch defined by a structure with fields 
%            'vertices' and 'faces'
%
% References
% -------------------------------------------------------------------------
% This function uses the surf2solid() and stlwrite() functions from Sven, 
% which are available from:
% https://www.mathworks.com/matlabcentral/fileexchange/42876-surf2solid-make-a-solid-volume-from-a-surface-for-3d-printing
% https://www.mathworks.com/matlabcentral/fileexchange/20922-stlwrite-write-ascii-or-binary-stl-files
%
% See some examples of MATLAB 3D printing @
% Reference: https://efcms.engr.utk.edu/ef230-2020-06/modules/3dprinting/matlab3d.php
%
% Author
% -------------------------------------------------------------------------
% Rafael Villamor Lora
% February 11th, 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            NOW THE CODE HERE ^^                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin <3 | nargin > 4
    error('Oooopsss... incorrect number of inputs. Check the function definition');
end
%--------------------------------------------------------------------------
% READ SETTINGS
outputFile   = settings.outputFile;   % Address of output files
geometry     = settings.geometry;     % Prismatic ('prism') or Cylindrical ('cylinder) fracture
bodyDim      = settings.bodyDim;      % Body dimensions as a fraction of the length [length, width, height OR diameter, height]
if nargin == 3                        % Position of the surface relative to the fracture plane
    position = "lowerSurface";        % If not specified: apply default value
elseif nargin == 4                    % If specified: check the spelling  
    if strcmp(position, "lowerSurface")
        position = "lowerSurface";
    elseif strcmp(position, "upperSurface")
        position = "upperSurface";
    else
        error('Ooooopsss ... surface position not recognized')
    end
end

%--------------------------------------------------------------------------
% SPECIMEN DIMENSIONS
[Nx, Ny] = size(surface); % Dimensions of the 2D matrix
if Ny > Nx                % We assume the length (longest dimension) is X-dir
    surface = surface';   % ... if not, transpose the matrix
end
N = Nx;  
switch geometry
    case 'prism'
        spcm_length   = dx * N * bodyDim.length;                           % Specimen length  [same units as dx]
        spcm_width    = dx * N * bodyDim.width;                            % Specimen width   [same units as dx]
        spcm_height   = dx * N * bodyDim.height;                           % Specimen height  [same units as dx]
    case 'cylinder'
        spcm_diameter = dx * N * bodyDim.diameter;                         % Specimen diameter [same units as dx]
        spcm_height   = dx * N * bodyDim.height;                           % Specimen height   [same units as dx]
    otherwise
        error('Geometry not supported');
end

%--------------------------------------------------------------------------
% CROP SURFACE TO DESIRED DIMENSIONS
switch geometry
    case 'prism'
        % Region of interest
        surface_cropped = surface(1:spcm_length/dx, 1:spcm_width/dx);
        % X-Y coordinates of surface
        [XX, YY]        = meshgrid(linspace(0, spcm_width, size(surface_cropped,2)), linspace(0, spcm_length, size(surface_cropped,1))); 
    case 'cylinder'
        surface_cropped = surface(1:spcm_diameter/dx, 1:spcm_diameter/dx);
        [XX, YY]     = meshgrid(linspace(0, spcm_diameter, size(surface_cropped,1)), linspace(0, spcm_diameter, size(surface_cropped,1)));       
    otherwise
        error('Geometry not supported');
end

%--------------------------------------------------------------------------
% SHIFT THE SURFACE TO THE APROPRIATE ELEVATION
if position == "lowerSurface"
    ZZ = surface_cropped + spcm_height;       % Translate
elseif position == "upperSurface"
    ZZ      = surface_cropped;
    rot_mat = [-1 0 0; 0 1 0; 0 0 -1];                                     % Rotation matrix along Y-direction (180 degrees)
    R       = rot_mat * [XX(:) YY(:) ZZ(:)]';                              % Flip the top surface (the printed base will be at the bottom)
    XX      = reshape(R(1,:), size(XX,1), size(XX,2));                     % New x-coordinate matrix
    YY      = reshape(R(2,:), size(YY,1), size(YY,2));                     % New y-coordinate matrix
    ZZ      = reshape(R(3,:), size(ZZ,1), size(ZZ,2)) + spcm_height;       % New z-coordinate matrix + Translate
end

%--------------------------------------------------------------------------
% FOR CYLINDRICAL SPECIMENS: GET A CIRCULAR SECTION
if geometry == "cylinder"
    [~,rho]     = cart2pol(XX - mean(XX(1,:)),YY - mean(YY(:,1)));         % Convert to polar coordinates
    outside     = find(rho > Cpcm_diameter/2);                             % Points outside the cylinder
    ZZ(outside) = 0;                                                       % Set points outside = 0
end

%--------------------------------------------------------------------------
% CONVERT SURFACE TO SOLID AND SAVE STL FILE
FV_fracture = surf2solid(XX, YY, ZZ, 'elevation', 0);                      % Convert to a closed solid
stlwrite(outputFile, FV_fracture);                                         % Save as .STL
end