% Formula found in paper "A study of degraded light coupling into
% single-mode fibers" - C. Ruilier - 1998
clearvars;

%%
lambda=1550e-9;
D=0.16;
w0=5e-6;

% Range of focal values
fini=D/0.6;
fend=D/0.01;
fpoints = 1000;
f = linspace(fini,fend,fpoints);

beta = pi/2*D./f*w0/lambda;
alpha = linspace(0,1,fpoints);

rho=@(beta,alpha) 2*(exp(-beta.^2).*(1-exp(beta.^2*(1-alpha^2)))./beta/sqrt(1-alpha^2)).^2;

figure()
for i=1:100:fpoints
    hold on;
    plot(beta,rho(beta,alpha(i)));
end
hold off;

%% NGS SOURCE
% Specific wavelengths defined in photometry.m
% Added bands for FSOC use

ngs = source('wavelength',photometry.FSOC_L2,'height',38000e3);
%% TELESCOPE
nPx = 280;
tel = telescope(D,... % Diameter
    'fieldOfViewInArcMin',0.0,...
    'obstructionRatio',0.0,... % Obstruction ratio
    'resolution',nPx,... % pupil sampling resolution
    'samplingTime',1/250);

%% Atmosphere
% paramters of the atmosphere, characterised by the site
f0 = 5e-2; % Fried parameter describes atmospheric turbulence strength, [m]
L0 = 30; % Outer scale length of turbulent structures, [m]

% use Cn2 profile for the fractioanlR0
fractionalR0 = [0.5, 0.3, 0.2]; % Contribution of each layer to the trubulence
% lower layers have higher turbulence and decreases with altitude

atm = atmosphere(photometry.V, f0, L0,...
    'altitude', [0, 4, 10]*1e3,'fractionnalR0', fractionalR0,...
    'windSpeed', [21, 21, 21], 'windDirection', [0, pi/3, pi]);

%% MAIN
% Fried parameter is given for 0.55um wavelength. Atmosphere class scales
% to appropriate wavelength.

% fried = [15e-2,6e-2,3e-2]; % Array containing Fried parameters
obsRatio = [0,0.2,0.34,0.4];
lambdaArray = [lambda];
% fried = 1e-2; % Single Fried value
Natm = 10; % Number of independent atmospheres to simulate for each Fried value

% Parameters
lambda = ngs.wavelength;
% w0=5e-6; % Fiber mode radius
D = tel.D; % Telescope diameter
res = tel.resolution; % Sampling of the telescope pupil

% Range of focal values
fini=tel.D/0.6;
fend=tel.D/0.01;
fpoints = 1000;
f = linspace(fini,fend,fpoints);

% Uncomment if only one focal to evaluate
% fopt=tel.D/0.22; % Infri case 

beta = pi/2*D./f*w0/lambda;
etaF_stored = zeros(numel(f),numel(obsRatio));
beta_opt = zeros(numel(obsRatio),1);
eta_opt = zeros(numel(obsRatio),1);
% etaF_fopt = zeros(Natm,numel(fried));
tic

for i = 1:numel(obsRatio)


tel = telescope(D,... % Diameter
    'fieldOfViewInArcMin',0.0,...
    'obstructionRatio',obsRatio(i),... % Obstruction ratio
    'resolution',nPx,... % pupil sampling resolution
    'samplingTime',1/250);

tel = tel + atm;

        ngs=ngs.*tel; % Propagation
        EAS = ngs.amplitude.*exp(1i.*ngs.phase); % Construction of the wave
        disp(ngs.phase); imagesc(ngs.phase);
        for k = 1:numel(f)
            % Efficiency
            [etaF_stored(k,i),~] = eta_Chen(EAS,w0,lambda,f(k),tel.D);
        end
        [pks, ind] = findpeaks(etaF_stored(:,i));
        beta_opt(i) = beta(ind);
        eta_opt(i) = etaF_stored(ind,i);
%         TF = islocalmax(etaF_stored(:,i));
%         beta_opt(i) = beta(TF);
%         eta_opt(i) = etaF_stored(TF,i);

end

toc

%%
figure
for i=1:numel(obsRatio)
    hold on;
    txt = ['\alpha = ',num2str(obsRatio(i)),sprintf('\n'),'(',num2str(beta_opt(i)),'  ',num2str(eta_opt(i)),')'];
    plot(beta,etaF_stored(:,i),'DisplayName',txt);
end
hold off;
ylabel('\eta');
xlabel('\beta');
legend show;

D_f = beta.*(2/pi).*(lambda/w0);
D_f_opt = beta_opt.*(2/pi).*(lambda/w0);

figure;
for i=1:numel(obsRatio)
    hold on;
    txt = ['\alpha = ', num2str(obsRatio(i)),sprintf('\n'),'(',num2str(D_f_opt(i)),'  ',num2str(eta_opt(i)),')'];
    plot(D_f,etaF_stored(:,i),'DisplayName',txt);
end
hold off;
ylabel('\eta')
xlabel('Relative aperture D/f')
legend show;
    
