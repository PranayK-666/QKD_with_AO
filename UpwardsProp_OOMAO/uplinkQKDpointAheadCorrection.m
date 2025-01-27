function [EE, LE_psf, pxSize_m] = uplinkQKDpointAheadCorrection(D,resolution,profile,sample_delay,correction,samplingTime,iterations,pointAhead,zen_scaleFactor)

% INPUTS:
%           D:                  telescope diameter [m]
%           resolution:         telescope resolution [px]
%           profile:            {'gaussian'};
%           sample_delay:       AO delay samples
%           correction:         1 - with AO; 0 - no correction
%           samplingTime:       1/1000 for 1KHz; 1/500 for 500 Hz
%           iterations:         number of iterations per simulation case
%           pointAhead:         uplink point ahead [arcsec]
%           zen_scaleFactor:    zen_scaleFactor takes into account
%                               elevation change. Use 1 for zenith. Use: 
%                               zen_scaleFactor = (1./cos(deg2rad(90-beta))).^(-3/5)
%                               for beta elevation 
%
% OUTPUTS:
%           EE:                 encircled energy
%           LE_psf:             long exposure PSF
%           pxSize_m:           pixel size in m at the satellite

r0 = zen_scaleFactor * 3e-2;             % r0 = 3e-2 m [m] 
% pointAhead = 10;                        % arcsec
nLenslet = 19; %16

s = RandStream('mt19937ar', 'seed', 5);

atmRx = atmosphere(photometry.V,r0,30,...
    'altitude',[0,4,10]*1e3,...
    'fractionnalR0',[0.7,0.05,0.02],...
    'windSpeed',[5,20,15],...
    'windDirection',[0,pi/4,-pi/4],...
    'randStream',s);
atmLaunch = atmosphere(photometry.V,r0,30,...
    'altitude',[0,4,10]*1e3,...
    'fractionnalR0',[0.7,0.05,0.02],...
    'windSpeed',[5,20,15],...
    'windDirection',[0,pi/4,-pi/4],...
    'randStream',s);

% atmRx = atmosphere(photometry.V,r0,30,...
%     'altitude',[0,4,10]*1e3,...
%     'fractionnalR0',[0.7,0.05,0.02],...
%     'windSpeed',[0,0,0],...
%     'windDirection',[0,pi/4,-pi/4],...
%     'randStream',s);
% atmLaunch = atmosphere(photometry.V,r0,30,...
%     'altitude',[0,4,10]*1e3,...
%     'fractionnalR0',[0.7,0.05,0.02],...
%     'windSpeed',[0,0,0],...
%     'windDirection',[0,pi/4,-pi/4],...
%     'randStream',s);

telRx = telescope(D,'resolution',resolution,'samplingTime',samplingTime,'fieldOfViewInArcsec',40);
telLaunch = telescope(D,'resolution',resolution,'samplingTime',samplingTime,'fieldOfViewInArcsec',40);


% deformable mirror influence function class
bifa = influenceFunction('monotonic',0.75);
% deformable mirror class
dm = deformableMirror(nLenslet+1,'modes',bifa,'resolution',telRx.resolution);

telRx = telRx + atmRx;
telLaunch = telLaunch + atmLaunch;

% SATELLITE PARAMETERS
orbit_height = 550e3;                   % 550 km
D_sat = 0.16;                           % satellite receiver D 16cm

%ENCIRCLED ENERGY ASSESSMENT
comms = source('wavelength',photometry.FSOC_Q2,  'zenith', arcsec(0), 'azimuth', deg2rad(0),'height',orbit_height);
% footprint = 2.44*(comms.wavelength/telLaunch.D).*orbit_height;         % Laser footprint @ sat [m] - full angle

w0 = (1-0.15)* D/2;                     % Assume 15% smaller than launch radius
zR = (pi*w0^2)/comms.wavelength;        % Rayleigh length (zR); n=1
% 
footprint = 2*w0*sqrt(1+(orbit_height/zR).^2);

% Assume no atm to calculate the number of px per footprint
atmIdeal = atmosphere(photometry.V,1000,30,...
    'altitude',[0,4,10]*1e3,...
    'fractionnalR0',[0.7,0.05,0.02],...
    'windSpeed',[5,20,15],...
    'windDirection',[0,pi/4,-pi/4],...
    'randStream',s);
telIdeal = telescope(D,'resolution',resolution,'samplingTime',samplingTime);
telIdeal = telIdeal + atmIdeal;

if strcmp(profile,'gaussian') 
    amplitude = tools.gaussian(telIdeal.resolution,0.8*telIdeal.resolution);
    amplitude = amplitude./max(amplitude(:));
else
    amplitude = 1;
end
comms.mask = telIdeal.pupil;
launchAtm = getLaunchAtm(telIdeal,telIdeal);
telIdeal.opticalAberration = launchAtm;

waveUpIdeal = upwardsPropagation_FSOC_QKD(comms,telIdeal,dm,amplitude);
figure(1);
imagesc(angle(waveUpIdeal));
title("phase")
colorbar; 

figure(2);
imagesc(abs(waveUpIdeal));
title("amplitude")
colorbar; pause(2);

% Calculate the FWHM of the footprint [px] in the diffractive-limit image at
% 550km to determine the px size.
frame_DL = abs(waveUpIdeal).^2;

% figure; imagesc(frame_DL); title('DL');

fwhm_DL = fwhm([1:telIdeal.resolution],frame_DL(size(frame_DL,1)/2,:));
w_1e2 = 2*0.8493218*fwhm_DL;              % relation between size at 1/e2 points and FWHM of a gaussian

pxSize_m = footprint/w_1e2;
% pxSize_m = footprint/fwhm_DL;

% Satellite receiver size in px
D_sat_px = D_sat/pxSize_m;
disp(D_sat_px);
% Sat area is represented by 'sat_rx'.
[sat_x,sat_y] = circle_coord(telLaunch.resolution/2, telLaunch.resolution/2, D_sat_px/2);

%%
circleSize_p = [1:1:telLaunch.resolution];

EE = zeros(iterations, size(circleSize_p,2));
TT = zeros(iterations,2);
LE_psf = zeros(iterations, telLaunch.resolution, telLaunch.resolution);


%NO CORRECTION
% figure;

if correction ==0

    w=waitbar(0,'No correction INITIALIZING...','Name','Processing Status');

    % Q) should each of this iteration corresponds to a different state?
    % I think these are propagation iterations and not AO correction
    % iterations
    for j=1:iterations

        waitbar(j/iterations,w,sprintf('D=%d m; point ahead = %d - No correction processing...%d of %d',[D,pointAhead,j,iterations]));
        if ~ishandle(w)
            break
        end
    
        dm.coefs = 0;

        % comms QKD uplink 810 nm (FSOC_Q2); downlink beacon 770 nm
        % (FSOC_Q3)
        comms = source('wavelength',photometry.FSOC_Q2,  'zenith', arcsec(0), 'azimuth', deg2rad(0),'height',orbit_height);
        
        if strcmp(profile,'gaussian') 
            amplitude = tools.gaussian(telLaunch.resolution,0.8*telLaunch.resolution);
            amplitude = amplitude./max(amplitude(:));
        else
            amplitude = 1;
        end
        comms.mask = telLaunch.pupil;
        
        % Upwards propagation
        +telRx;
        launchAtm = getLaunchAtm_FOV(telLaunch,telRx);
        telLaunch.opticalAberration = launchAtm;

        % add loop for each coherent state

        waveUp = upwardsPropagation_FSOC_QKD(comms,telLaunch,dm,amplitude);
        figure(3)
        imagesc(angle(waveUp)); 
        title("Phase"); colorbar;
        figure(4)
        imagesc(abs(waveUp));
        title("Amplitude"); colorbar;

        I2 = abs(waveUp).^2;
        
        % Downwards propagation
        down_beacon = source('wavelength',photometry.FSOC_Q3,  'zenith', arcsec(pointAhead), 'azimuth', deg2rad(0),'height',orbit_height);
        % lgsDown.extent = abs(waveUp);
        down_beacon = down_beacon.*telRx;

        % loop around for each coherent state
        % each coherent state will give one X and P

        LE_psf(j,:,:) = I2;

        % here goes the measurement of T using X and P distribution
        
        % EE @ sat (comms uplink)
        for n = 1:length(circleSize_p)
            EE(j,n) = encircledEnergy(I2, circleSize_p(n),'centralPx');
        end
        % % TT @ sat (comms uplink)
        % [reconstructedZernike, residual] = lf_valdivia_zerfit(comms.phase,3,telRx.pupil);
        % TT(j,:) = reconstructedZernike(2:3);
        % % 
        % imagesc(I2);
        % drawnow;

        
    end

    delete(w)

%CORRECTION
else

    w=waitbar(0,'AO correction INITIALIZING...','Name','Processing Status');
    
    % comms QKD uplink 810 nm (FSOC_Q2); downlink beacon 770 nm
    % (FSOC_Q3)

    %Propagation
    comms = source('wavelength',photometry.FSOC_Q2,  'zenith', arcsec(0), 'azimuth', deg2rad(0),'height',orbit_height);
    
    if strcmp(profile,'gaussian')
        amplitude = tools.gaussian(telLaunch.resolution,0.8*telLaunch.resolution);
        amplitude = amplitude./max(amplitude(:));
        threshold = 1;
    else
        amplitude = 1;
        threshold = 5;
    end
    comms.mask = telLaunch.pupil;
    
    % Upwards propagation
    if sample_delay==1
        +telRx;
    % else
    %     +telRx;
    %     +telRx;
    end
    
    % Upwards propagation (QKD uplink @ 810 nm)
    launchAtm = getLaunchAtm_FOV(telLaunch,telRx);
    telLaunch.opticalAberration = launchAtm;
    waveUp = upwardsPropagation_FSOC_QKD(comms,telLaunch,dm,amplitude);
    

    % Downwards propagation (downlink beacon @ 770nm)
    down_beacon = source('wavelength',photometry.FSOC_Q3,  'zenith', arcsec(pointAhead), 'azimuth', deg2rad(0),'height',orbit_height);
    % lgsDown.extent = abs(waveUp);
    down_beacon = down_beacon.*telRx;
    
    for j=1:iterations

        waitbar(j/iterations,w,sprintf('D=%d m; point ahead = %d - AO correction processing...%d of %d',[D,pointAhead,j,iterations]));
        if ~ishandle(w)
            break
        end
    
        dm.coefs = 0;

        %Perfect correction - perfect measurement in the downlink beacon
        perfectMeasurement = getPhaseDown(telLaunch,telRx,down_beacon.phase); % down_beacon.phase
        modes=dm.modes.modes;
        coefs=modes\(perfectMeasurement(:)*(down_beacon.wavelength/(2*pi)));
        
        dm.coefs = 0.5*coefs;

        %Propagation
        comms = source('wavelength',photometry.FSOC_Q2,  'zenith', arcsec(0), 'azimuth', deg2rad(0),'height',orbit_height);
        
        if strcmp(profile,'gaussian')
            amplitude = tools.gaussian(telLaunch.resolution,0.8*telLaunch.resolution);
            amplitude = amplitude./max(amplitude(:));
            threshold = 1;
        else
            amplitude = 1;
            threshold = 5;
        end
        comms.mask = telLaunch.pupil;
           
         % Upwards propagation
        if sample_delay==1
            +telRx;
        else
            +telRx;
            +telRx;
        end

        % Upwards propagation
        launchAtm = getLaunchAtm_FOV(telLaunch,telRx);
        telLaunch.opticalAberration = launchAtm;
        waveUp = upwardsPropagation_FSOC_QKD(comms,telLaunch,dm,amplitude);
        % imagesc(abs(waveUp));
        figure(3)
        imagesc(angle(waveUp)); 
        title("Phase"); colorbar;
        figure(4)
        imagesc(abs(waveUp));
        title("Amplitude"); colorbar;
        % Downwards propagation
        down_beacon = source('wavelength',photometry.FSOC_Q3,  'zenith', arcsec(pointAhead), 'azimuth', deg2rad(0),'height',orbit_height);
        % lgsDown.extent = abs(waveUp);
        down_beacon = down_beacon.*telRx;

        frame = abs(waveUp).^2;
        % frame(frame<threshold) = 0; %1 as a threshold for the EE computation

        LE_psf(j,:,:) = frame;

        % EE @ sat (comms uplink)
        for n = 1:length(circleSize_p)
            EE(j,n) = encircledEnergy(frame, circleSize_p(n), 'centralPx');
        end
        % % TT % sat (comms uplink)
        % [reconstructedZernike, residual] = lf_valdivia_zerfit(comms.phase,3,telRx.pupil);
        % TT(j,:) = reconstructedZernike(2:3);

        % imagesc(frame);
        % drawnow;
        % 
    end
    
end

delete(w)

EE = mean(EE,1);
LE_psf = sum(LE_psf,1);

disp(D_sat_px);

end