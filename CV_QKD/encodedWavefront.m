% Create coherent states with Gaussian modulation as wavefronts
% Simulate Gaussian modulated coherent quantum states by displacing
% vaccum states [Weedbrook (2004), Laudenbach (2018)]
% 
% [Weedbrook, C. et al. (2004) ‘Quantum Cryptography Without
% Switching’, Physical Review Letters, 93(17), p. 170504. Available at: 
% https://doi.org/10.1103/PhysRevLett.93.170504.]
% 
% Laudenbach, F. et al. (2018) ‘Continuous‐Variable Quantum Key 
% Distribution with Gaussian Modulation—The Theory of Practical 
% Implementations’, Advanced Quantum Technologies, 1(1), p. 1800011. 
% Available at: https://doi.org/10.1002/qute.201800011.
% 
% this code needs to be modified to use SI units to make it usable for
% atmpshperic propagation, etc.
% 
% date created: 11/01/25
% date modified: 27/01/25
% 
% author: Pranay Khosla
% 

% standard preface
close all; clear; clc;
% clear all; % run when facing problems or running first time

% rng("default")

%% Parameters
% set these parameters such that numerical error is minimised
% N = 1000, paths = 100 seems to work
% since the code is vectorised, this may not be a lot of computaitonal load
paths = 100; % trajectories of each state (quantumness)
N = 1000; % # Gaussian modulated states (alpha_i)
res = 16; % use a small resolution for testing, might need 128 px for simulations

%% Quantum operators for a wavefront
% realisations of vacuum field: (1/2) * (q + ip)
vac = @(paths,res) (1/2) * (randn(res,res,paths) + 1i*randn(res,res,paths));

% quadrature operators: (a + a*), -i(a - a*)
X = @(a) (a + conj(a));
P = @(a) -1i*(a - conj(a));

%% Wavefront with vacuum in each pixel

var_mod = 5;
% displacement parameters
x = normrnd(0,sqrt(var_mod),1);
p = normrnd(0,sqrt(var_mod),1);

%% Telescope pupil
pupil = @(res) sqrt(meshgrid(linspace(-1,1,res)).^2 + meshgrid(linspace(-1,1,res))'.^2) <= 1;

% imagesc(x*pupil(res));
% axis square;

%% Quadratures method
vac_wavefront = vac(paths,res);
pupilX = x*pupil(res) + X(vac_wavefront);
pupilP = p*pupil(res) + P(vac_wavefront);
%%
meanX = (1/paths)*sum(pupilX,3);
meanP = (1/paths)*sum(pupilP,3);

varX = (1/paths)*sum(pupilX.*pupilX,3) - meanX.^2;
varP = (1/paths)*sum(pupilP.*pupilP,3) - meanP.^2;

pupilAmpl = mean(sqrt(pupilX.^2+pupilP.^2),3);
% which one to use?
% |alpha> = |q +ip>
% a = 1/2 * (q +ip)
pupilAlpha = (1/2) * (pupilX + 1i*pupilP);
pupilAmpl2 = mean(abs(pupilAlpha),3);
pupilAmpl3 = mean(sqrt(pupilAlpha .* conj(pupilAlpha)),3);


%% State method 
% (mean is ambiguous, there is a factor of 2)
% (var is fine)
% if alpha0 = 0, then also there is a mean of 0.6
% alpha0 = sqrt(x^2+p^2);
% phi = atan(p/x);
% 
% pupilAlpha = (alpha0*pupil(res) + vac_wavefront) .* exp(1i*phi*pupil(res));
% 
% alpha = pupilAlpha;
% pupilX2 = (alpha + conj(alpha));
% pupilP2 = -1i*(alpha - conj(alpha));
% 
% meanX2 = mean(pupilX2,3);
% meanP2 = mean(pupilP2,3);
% 
% varX2 = (1/paths)*sum(pupilX2.*pupilX2,3) - meanX2.^2;
% varP2 = (1/paths)*sum(pupilP2.*pupilP2,3) - meanP2.^2;
% 
% meanAlpha0 = mean(abs(pupilAlpha),3);
% meanPhi = mean(angle(pupilAlpha),3);