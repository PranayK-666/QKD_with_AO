% standard preface
close all; clc;
clear;
% clear ll; % use when facing problems

%% Uplink Simulation

% simulation parameters
D = 0.7; % transmitter telescope (ground station) diameter [m]
res = 128; % telescope resolution [px]
profile = 'gaussian';
sample_delay = 1;
samplingTime = 1/500; % 500 Hz
iter = 10; % number of iterations per simulation case
pointAhead = 2; % uplink point ahead [arcsec]
zen_scaleFactor = 5/3; % zenith

%% without AO correction
correction  = 0; % 1 - with AO; 0 - no correction
[EE, LE_psf, pxSize_m] = uplinkQKDpointAheadCorrection(D, res, profile, sample_delay, correction, samplingTime, iter, pointAhead, zen_scaleFactor);
%% good atmosphere
correction  = 0;
zen_scaleFactor = 10; % used this value to simulate ideal atmosphere (very large r0)
[EE2, LE_psf_2, pxSize2] = uplinkQKDpointAheadCorrection(D, res, profile, sample_delay, correction, samplingTime, iter, pointAhead, zen_scaleFactor);
%% with AO correction
zen_scaleFactor = 5/3;
correction  = 1; % 1 - with AO; 0 - no correction
[EE_AO, LE_psf_AO, pxSize_m_AO] = uplinkQKDpointAheadCorrection(D, res, profile, sample_delay, correction, samplingTime, iter, pointAhead, zen_scaleFactor);

%% Encircled energy
figure;
plot(EE, 'r--', 'LineWidth', 2)
hold on
plot(EE_AO, 'b', 'LineWidth', 2)
hold on
plot(EE2, 'g-.', 'LineWidth', 2)
legend('No Correction', 'Corrected with AO', 'Good Atmosphere')

%% Save PSF

% without AO correction
filename = 'LE_psf_animation.gif'; 

figure; colormap('pink');

for i = 1:size(LE_psf, 1)
    % Extract the i-th image from LE_psf
    image_data = squeeze(LE_psf(i, :, :));

    imagesc(image_data); colorbar;
    title(['Image from LE\_psf - Frame ' num2str(i)]);
    axis image; 

    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    % Write to the GIF File
    if i == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end

    % Pause to create the animation effect
    pause(0.05);
end

% good atmosphere
filename = 'LE_psf_animation2.gif'; 

figure; colormap('pink');

for i = 1:size(LE_psf_2, 1)
    % Extract the i-th image from LE_psf
    image_data = squeeze(LE_psf_2(i, :, :));

    imagesc(image_data); colorbar;
    title(['Image from LE\_psf - Frame ' num2str(i)]);
    axis image;

    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    % Write to the GIF File
    if i == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end

    % Pause to create the animation effect
    pause(0.05);
end


% wih AO correction
filename = 'LE_psf_AO_animation.gif';

figure; colormap('pink');

for i = 1:size(LE_psf_AO, 1)
    % Extract the i-th image from LE_psf
    image_data = squeeze(LE_psf_AO(i, :, :));

    imagesc(image_data); colorbar;
    title(['Image from LE\_psf - Frame ' num2str(i)]);
    axis image;

    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    % Write to the GIF File
    if i == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end

    % Pause to create the animation effect
    pause(0.05);
end
