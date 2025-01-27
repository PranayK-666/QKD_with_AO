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
% date created: 18/12/24
% date modified: 27/01/25
% 
% author: Pranay Khosla
% 

% standard preface
close all; clear; clc;
% clear all; % run when facing problems or running first time

%% Parameters
% set these parameters such that numerical error is minimised
% N = 1000, paths = 100 seems to work
% since the code is vectorised, this may not be a lot of computational load
paths = 1000; % trajectories of each state (quantumness)
N = 100; % # Gaussian modulated states (alpha_i)

%% Quantum operators
% N realisations of vacuum field: (1/2) * (q + ip)
vac = @(paths,N) (1/2) * (randn(N,paths) + 1i*randn(N,paths));

% quadrature operators: (a + a*), -i(a - a*)
X = @(a) (a + conj(a));
P = @(a) -1i*(a - conj(a));

%% Displaced coherent states with Gaussian modulation
% Samples N points from zero-centred Gaussian distribution with given variance
Gaussian = @(V,N) normrnd(0,sqrt(V),1,N);

% draw real numbers from Gaussian distribution
var_mod = 5;
x = Gaussian(var_mod,N);
p = Gaussian(var_mod,N);

% x=0;p=0;

vac_noise = vac(paths,N);
% displace quadratures of vacuum state by the Gaussian sampled real numbers
X_Gaussian = x' + X(vac_noise);
P_Gaussian = p' + P(vac_noise);

meanX_Gaussian = mean(X_Gaussian,2);
meanP_Gaussian = mean(P_Gaussian,2);

% take care of statisitcs here, calculate variance (scale by N or N-1)
varX_Gaussian = (1/paths)*sum(X_Gaussian.*X_Gaussian,2) - meanX_Gaussian.^2;
varP_Gaussian = (1/paths)*sum(P_Gaussian.*P_Gaussian,2) - meanP_Gaussian.^2;

% vac_noise = vac(paths,N);
% X_vac = X(vac_noise);
% P_vac = P(vac_noise);
% meanX = mean(X_vac,2);
% meanP = mean(P_vac,2);
% varX = (1/paths)*sum(X_vac.*X_vac,2)-meanX.^2;
% varP = (1/paths)*sum(P_vac.*P_vac,2)-meanP.^2;

%% Plot the states in phase space

% using for loop for plotting is slower but looks nicer
figure(1); clf
for i = 1:N
    plot(X_Gaussian(i,:),P_Gaussian(i,:),'.')
    hold on
    plot(meanX_Gaussian(i),meanP_Gaussian(i),'kx','MarkerSize',8,'LineWidth',2)   
end

% % faster plotting using vectorised code (kind of)
% figure(2); clf
% cmap = lines(N);
% scatter(X_Gaussian(:),P_Gaussian(:),10,repmat(cmap, paths,1),'filled',MarkerFaceAlpha=0.4);
% hold on
% scatter(meanX_Gaussian,meanP_Gaussian,40,'kx','LineWidth',2,'MarkerEdgeAlpha',0.4);

xmax = 6*sqrt(var_mod); % std = sqrt(var), -3*std to + 3*std must contain 99.7%
axis equal
axis([-xmax xmax -xmax xmax])

% Plot x- and y-axes as hard lines
line([-xmax, xmax], [0, 0], 'Color', 'k'); % x-axis
line([0, 0], [-xmax, xmax], 'Color', 'k'); % y-axis

xlabel('X - Quadrature')
ylabel('P - Quadrature')
title('Gaussian modulated coherent states')

%% Verifying the simulation

% 1) V(x,p) = V_mod (used for Gaussian distribution)

fprintf('Variance of the quadrature components x and p should match the variance of the Gaussian distribution\n V_mod: %0.3f\n',var_mod)
fprintf('V(x): %.3f, V(p): %.3f\n',var(x),var(p));

% 2) <n> = 2V_mod (for the ensemble); numerically this should equal var(x) +
% var(p)

fprintf('\nMean photon number of the ensemble should be twice the modulation variance, 2V_mod: %.3f\n',2*var_mod)
alpha0 = sqrt(meanX_Gaussian.^2+meanP_Gaussian.^2);
meanPhotonNumber = mean(alpha0.*alpha0);
fprintf('Mean of alpha mod squared is: %.3f\n',meanPhotonNumber)

% 3) V(X,P) = V(x,p) + 1 (variance of the ensemble)
% here, we assume that the whole ensemble is 1 state with 0 mean photon
% number and calculate the variance (looks like huge blob at the origin,
% similar to vaccum state with higher variance

X1 = reshape(X_Gaussian,1,[]);
var_X1 = 1/(paths*N)*sum(X1.*X1)-((1/(paths*N))*sum(X1))^2;
disp(var_X1 - var(x))

P1 = reshape(P_Gaussian,1,[]);
var_P1 = 1/(paths*N)*sum(P1.*P1)-((1/(paths*N))*sum(P1))^2;
disp(var_P1 - var(p))

%% Save the Gaussian states to a file
% save('X_Gaussian.mat','X_Gaussian')
% save('P_gaussian.mat','P_Gaussian')
% save('GaussianQuadratures.mat','X_Gaussian','P_Gaussian')

%% Plot one alpha
figure(3); clf
subplot(2,2,1)
i = 55; % 55
plot(X_Gaussian(i,:),P_Gaussian(i,:),'.')
hold on
plot(meanX_Gaussian(i),meanP_Gaussian(i),'rx','MarkerSize',8,'LineWidth',2)

xmax = 4*sqrt(var_mod);
axis equal
axis([-xmax xmax -xmax xmax])

% Plot x- and y-axes as hard lines
line([-xmax, xmax], [0, 0], 'Color', 'k'); % x-axis
line([0, 0], [-xmax, xmax], 'Color', 'k'); % y-axis

xlabel('X - Quadrature')
ylabel('P - Quadrature')
title('Coherent state')

subplot(2,2,2)
hP = histogram(meanP_Gaussian,'Orientation','horizontal');
xlabel('PDF (P)')

subplot(2,2,3)
hX = histogram(meanX_Gaussian);
ylabel('PDF (X)')

subplot(2,2,4)
for i = 1:N
    plot(X_Gaussian(i,:),P_Gaussian(i,:),'.')
    hold on
    plot(meanX_Gaussian(i),meanP_Gaussian(i),'kx','MarkerSize',4,'LineWidth',2)   
end

xmax = 4*sqrt(var_mod);
axis equal
axis([-xmax xmax -xmax xmax])

% Plot x- and y-axes as hard lines
line([-xmax, xmax], [0, 0], 'Color', 'k'); % x-axis
line([0, 0], [-xmax, xmax], 'Color', 'k'); % y-axis

xlabel('X - Quadrature')
ylabel('P - Quadrature')
title('Gaussian modulated coherent states')