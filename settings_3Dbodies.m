function [settings] = settings_3Dbodies()
% This function generates the necessary settings for the fields2bodies
% function
outputFile = "Results/surface.stl"; % Output file address
geometry   = "prism";               % Specimen geometry: 'prism' or 'cylinder'
switch geometry
    case "prism"
        bodyDim.length = 1;    % Body dimensions as a fraction of the length (i.e., values 0-1) [Length]
        bodyDim.width  = 0.5;  % Body dimensions as a fraction of the length (i.e., values 0-1) [Width]
        bodyDim.height = 0.10; % Body dimensions as a fraction of the length (i.e., values 0-1) [Height]
    case "cylinder"
        bodyDim.diameter = 1;  % Upper body dimensions as a fraction of the length (i.e., values 0-1) [Diameter]
        bodyDim.height = 1;    % Upper body dimensions as a fraction of the length (i.e., values 0-1) [Height]

    otherwise
        error('Geometry not supported');
end
% Outputs
settings.geometry    = geometry;
settings.bodyDim     = bodyDim;
settings.outputFile  = outputFile;
end

