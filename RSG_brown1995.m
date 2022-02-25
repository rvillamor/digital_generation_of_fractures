function [Zupper, Zlower, Zaperture, PSD] = RSG_brown1995(H, roughness, mismatch, N, aniso, seed, lambda_0, model)
% SIMPLE MATHEMATICAL MODEL OF A ROUGH FRACTURE (Brown, 1995)
%
% [Zupper, Zlower, Zaperture, PSD] = RSG_brown1995(H, roughness, mismatch, N, aniso, seed, lambda_0, model)
%
% This function implements the model described in Brown (1995):
% "A simple mathematical model of rough-walled fractures in rock is 
% described which requires the specification of only three main parameters:
% the fractal dimension, the rms roughness at a reference length scale, and
% a length scale describing the degree of mismatch between the two fracture
% surfaces"
%
% Inputs
% -------------------------------------------------------------------------
% H              : Hurst exponent. Determines fractal dimension D=3-H                                  [0 < H < 1]
% roughness      : root-mean-square roughness (i.e.standard deviation of heights)                      [m2]
% mismatch       : Mismatch length scale (wavelength) as a fraction of fracture size                   [0 < Mismatch < 1]                      
% N              : Number of divisions along fracture side                                             [power of 2]               
% aniso          : anisotropy ratio                                                                    [0 < Aniso < 1]
% seed           : random number seed                                                                  [+ve integer]
% lambda_0       : (optional) "roll-off" length scale as a fraction of fracture size                   [0 < lambda_0 < 1 (default)]
% model          : (optional) Power spectral density model: 'linear' (default),'bilinear', 'smooth'
%
% Outputs
% -------------------------------------------------------------------------
% Zupper, Zlower : Upper and lower rough surfaces
% Zaperture      : Aperture field map
% PSD            : Structure containing the PSD matrices (qxx, qyy, Cq)
% 
% References:
% -------------------------------------------------------------------------
% Brown, S.R. (1995). Simple mathematical model of a rough farcture. JGR
%    100(B4) 5941-5952
%
% Author:
% -------------------------------------------------------------------------
% Rafael Villamor Lora
% May 5, 2020 (last modified: 06/18/2020)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            NOW THE CODE HERE ^^                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% LOCAL VARIABLES
MaxSeed  = 10; % Maximum random number seed value
L        = 1;  % Fracture length (for now this code is scale ignorant --> 
               % the user will define the RMS at 'current' scale)

% CHECK THE FUNCTION INPUTS
% Check number of inputs
if nargin <6 | nargin > 8
    error('Upsss... incorrect number of inputs. Check the function definition');
end
% Apply default values
if nargin == 6
    lambda_0 = 1;        % No "roll-off"
end
if nargin < 8
    model    = 'linear'; % Linear model
end
% Check mismatch length scale
if mismatch < 0 | mismatch > 1
    error('Upsss... Check your inputs! [0 < mismatch < 1] ');
end
% Check anisotropy ratio
if aniso < 0 | aniso > 1
    error('Upsss... Check your inputs! [0 < aniso < 1] ');
end
% Check the "roll-off" length scale
if lambda_0 < 0 | lambda_0 > 1
    error('Upsss... Check your inputs! [0 < lambda_0 < 1]');
end
% Check the PSD model
if strcmp(model, 'linear') + strcmp(model,'bilinear') + strcmp(model,'smooth') < 1
    error('Upsss... Check your input! model = "linear", "bilinear" or "smooth"');
end

%-------------------------------------------------------------------------%
% GENERAL APPROACH
%
% We are going to generate two surface topographies using the inverse fast 
% Fourier transform, ifft(). And for that we need two basic ingredients 
% (for each surface)
%
%         The modulus of the Fourier transform, |F| --> Fmodulus
%                            +
%             The phase of the Fourier transform    --> Fphase
%
% Fmodulus is going to be the same for both surfaces, and the Fphase will
% be too except for the 'small' wavelengths (i.e. all those wavelengths <
% the mismatch length scale).
%
% Then, we will rescale the surfaces in the z direction so we can impose
% the input rms roughness. In other words, we will generate surfaces with
% a generic scale, and then we will specify the scale at that particular 
% scale. 

%-------------------------------------------------------------------------%
% INITIATE RANDOM NUMBER GENERATORS
% These are two independent random number streams that use the uniform 
% pseudorandom algorithm specified by 'gentype' -> 'mt19937ar'
% Here myStream1 and myStream2 use the same uniform pseudorandom number 
% generator algorithm (i.e. 'mt19937ar'), but use different seeds
myStream1 = RandStream('mt19937ar', 'Seed', seed);
myStream2 = RandStream('mt19937ar', 'Seed', MaxSeed * rand(myStream1));

%-------------------------------------------------------------------------%
% BUILD POWER SPECTRAL DENSITY MATRIX ==> GET DFT MODULUS, |F|
% Wavenumbers
dx          = L/N;                                % Sample interval
Nyq         = 1/(2*dx);                           % Nyquist of data
dq          = 1/(N*dx);                           % Wavenumber increment
wavenumbers = 2 * pi * [-Nyq : dq : Nyq-dq];      % 1D Wavenumbers     
[qx, qy]    = meshgrid(wavenumbers);              % 2D Wavenumbers [m^-1]
q           = (qx.^2 + (qy/aniso).^2).^0.5;       % Wavenumber modulus (Apply here the anisotropy conditions)
lambda      = 2*pi ./ q;                          % Wavelengths modulus

% Find wavelengths above the "roll-off" lengthscale
rollOff     = find(lambda > lambda_0 * L);

% An interesting characteristic of self-affine surfaces is that they have 
% decreasing power-law PSD across most scales. For the case of the C^{iso}, 
% we have 
%
%                  C^(iso) = Cq(q) ~ q^{-2-2H)
%
% where H is the Hurst exponent [0<H<1], which is related to the fractal
% dimension by D = 3 - H. And q = sqrt(qx.^2 + qy.^2).
%
% Now, we know that the Cq is proportional to the square of the amplitude
% spetrum [A(q)], i.e.
%
%                  Cq(q) ~ [A(q)].^2
%
% And we also know that the amplitude spectrum is proportional to the
% modulus of the Fourier transform of the signal, |F|,
%
%                  |F| ~ A(q)
%
% which yields to
%
%                  Cq(q) ~ |F|.^2 --> |F| ~ [Cq(q)]^0.5
%
% and finally we have
%
%                  |F| ~ [q^{-2-2H)]^0.5 = q^{-1-H}
%
% Strictly speaking the above relationships are of the form of
%
%                  Cq(q) = (dx / N)^2 .* |F|.^2
%                  Cq(q) = C0 .* q^{-2-2H)
%
% and therefore we should have
%
%                 |F| = sqrt(C0) * N / dx .* q^{-1-H} 
%
% where C0 is a proprotionality constant that defines the roughness scale.
% See Parsson et al., 2005 or Appendix X for more details.

% However, in this implementation we will apply the scaling later using the
% parameter 'roughess'. And therefore, we use the form: |F| ~ q^{-1-H}

% For the anisotropy case, the modulus takes the form of
if     strcmp(model, 'linear')
       % Compute the FFT modulus
       Fmodulus          = q .^ (-1-H);                                    % |F| ~ q^{-1-H}
       %Fmodulus         = sqrt(C0) * N / dx .* q .^ (-1-H);               % |F| = sqrt(C0) * N / dx .* q^{-1-H}
       % Apply "roll-off" conditions
       Fmodulus(rollOff) = 0;                        
elseif strcmp(model, 'bilinear')
       % Compute the FFT modulus
       Fmodulus          = q .^ (-1-H);                                    % |F| ~ q^{-1-H}
       % Apply "roll-off" conditions
       Fmodulus(rollOff) = (2*pi / (lambda_0 * L)) .^ (-1-H);
elseif strcmp(model, 'smooth')
       % Compute the FFT modulus + Apply "roll-off" conditions
       Fmodulus          = (2*pi / (lambda_0 * L) + q ) .^ (-1-H);         % |F| ~ (q0 + q)^{-1-H}
       warning('here')
end

% Set the power of short wavenumbers (i.e. long wavelengths) to zero
zAmp                                                   = 0; % This number should be >= 0 (use zAmp = 0 --> for k0)
Fmodulus(N/2+1-zAmp:N/2+1+zAmp, N/2+1-zAmp:N/2+1+zAmp) = 0;
% Also along qx,qy =0 --> not sure if I need this
% Fmodulus(N/2+1, :) = 0;
% Fmodulus(:, N/2+1) = 0;


%-------------------------------------------------------------------------%
% BUILD PHASE SPECTRAL DENSITY MATRIX==> GET DFT angle, ang(F)
% Use the two independent random number streams to build two independent
% phase spectra
phi1 = 2 * pi * rand(myStream1,N);     
phi2 = 2 * pi * rand(myStream2,N);
% Make phase spectra to be symmetric with respect the origin --> Not sure if I need this
%{
phi1_up = rand(myStream1,N/2,N);
phi1    = 2 * pi * [phi1_up; -flipud(fliplr(phi1_up))];
phi2_up = rand(myStream2,N/2,N);
phi2    = 2 * pi * [phi2_up; -flipud(fliplr(phi2_up))];
%}

% Find mismatching wavelengths
mismatching = find(lambda < mismatch*L);  % Find mismatching wavelengths

% Build the phase spectral matrix for each surface
Phase1              = phi1;               % First surface
Phase2              = Phase1;             % Second surface
Phase2(mismatching) = phi2(mismatching);  % The phase spectrum of the second surface is different for short wavelengths

%-------------------------------------------------------------------------%
% BUILD A COMPLEX ARRAY TO INPUT INTO THE FFT2
% Note both surfaces share the same amplitude spectrum, they only differ in
% the phase spectrum @ small wavelengths
A1 = Fmodulus .* complex(cos(Phase1), sin(Phase1)) ./ abs(complex(cos(Phase1), sin(Phase1))); % First surface
A2 = Fmodulus .* complex(cos(Phase2), sin(Phase2)) ./ abs(complex(cos(Phase2), sin(Phase2))); % Second surface

% COMPUTE THE SURFACES USING FFT2
% Note 1: X = ifftshift(Y) rearranges a zero-frequency-shifted Fourier 
%         transform Y back to the original transform output. I.e. it shifts 
%         the zero-frequency component from the center of the array to the
%         corners.
% Note 2: Since A is a  nearly conjugate symmetric matrices, you can 
%         compute the inverse Fourier transform faster by specifying the 
%         'symmetric' option, which also ensures that the output is real.
Z1  = ifft2(ifftshift(A1),'symmetric');  % Upper surface
Z2  = ifft2(ifftshift(A2),'symmetric');  % Lower surface

%-------------------------------------------------------------------------%
% SET THE STANDARD DEVIATION OF THE SURFACE HEIGHTS
% Recall that so far we have not specify the length scale of the surface...
% and actually we will not :)
%
% What we are going to do is the following:
% 1) Assume that the generated surfaces have a length of L_generated, which
%    surprisingly matches the surface length the user had in mind, i.e.
%    L_generated = L_user = L
% 2) But the standard deviation of the generated surface, std(Z) does not
%    match the std (i.e. the rms roughness) the user wanted... :(
% 3) We will fix this by saying Z = Z / std(Z) * std_user

% Note that std(Z1(:)) = std(Z2(:)) since the std(Z) (or the rms roughness)
% depends only of the Cq(q), which is the same for both surfaces,
roughnessScale = roughness / std(Z1(:)); 

% And now apply it...
Z1             = Z1 * roughnessScale;
Z2             = Z2 * roughnessScale;

%-------------------------------------------------------------------------%
% COMPUTE APERTURE
% Our upper and lower surfaces are perfectly aligned, so the composite
% aperture is given by
Zap = (Z1 - Z2);                   
Zap = Zap - min(Zap(:));

%-------------------------------------------------------------------------%
% COMPUTE THE 2DPOWER SPECTRUM DENSITY OF THE RESULTING FRACTURE
% Check Parsson et al., 2005 – Appendices A and D if you don't know
% where the following is coming from. It is also in the appendices of my 
% thesis. But basically,
%
% Cq(q) = (Nxdx * Nydy) * A(q).^2
% Cq(q) = (Ndx)^2 * [|F| / N^2].^2
% Cq(q) = dx^2 / N^2 * |F|.^2
%
% And recall to scale the Fmodulus we used before  with the
% 'roughnessScale' factor
Fmodulus = Fmodulus * roughnessScale;
Cq       = Fmodulus.^2 * (dx / N)^2;

%-------------------------------------------------------------------------%
% OUTPUTS

% Surfaces
Zupper    = Z1;                       % Upper surface
Zlower    = Z2;                       % Lower surface
Zaperture = Zap;                      % Aperture field

% Power Spectral Density
PSD.qxx = qx;                         % x-wavevectors
PSD.qyy = qy;                         % y-wavevectors
PSD.Cq  = Cq;                         % 2D PSD

end
