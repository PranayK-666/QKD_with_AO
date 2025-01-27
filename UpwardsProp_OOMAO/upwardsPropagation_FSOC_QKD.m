function out = upwardsPropagation_FSOC_QKD(src,tel,dm,amplitudeIn)
        
            atm_m           = tel.opticalAberration;
            nLayer          = atm_m.nLayer;
            layers          = atm_m.layer;
            altitude_m      = [layers.altitude];
            phase_m         = [ layers.phase ];
            layersNPixel    = [ layers.nPixel ];
            srcHeight = src.height;
            wavenumber      = 2*pi/src.wavelength;
            out_ = zeros(tel.resolution,tel.resolution,nLayer);
            out = zeros(tel.resolution,tel.resolution,nLayer);
            nOut = size(out,1);

            %%%%%%%%%%% Create the QKD Phase %%%%%%%%%%%%
            res = tel.resolution;
            pupil = @(res) sqrt(meshgrid(linspace(-1,1,res)).^2 + meshgrid(linspace(-1,1,res))'.^2) <= 1;
            QKDPhase = (zeros(res)+101*pi).*pupil(res);
            %%%%%%%%%%% Create the QKD Phase %%%%%%%%%%%%
        
            for kLayer = 1:nLayer

                height = altitude_m(kLayer);
%                 [xs,ys] = meshgrid(layerSampling_m{kLayer});
                if height==0 && layersNPixel(kLayer)==nOut
                    % FIRST LAYER
                    %Add atm layer
                    out_(:,:,kLayer) = phase_m{kLayer};
                    out(:,:,kLayer) = out_(:,:,kLayer)*(tel.opticalAberration.wavelength/src(1).wavelength);
                    
                    %Add mirror focus to 90km (as a fixed value in order to check with other srcHeight the spot will be larger)
                    ngsTemp = source('wavelength', src(1).photometry);
                    telTemp = telescope(tel.D,'resolution',tel.resolution);
                    zern = zernike(4, tel.D, 'resolution', tel.resolution);
                    zern.c = 1.635*(-ngsTemp.wavelength/(2*pi))*(90e3*(1-cos(tel.R/90e3))/ngsTemp.wavelength); %2.35 /0.94//1.635%value??? defocus = z*(1-cos(tel.R/90e3)) [m] 
%                     zern.c = 1.635*(-ngsTemp.wavelength/(2*pi))*((90e3*(1-cos(tel.R/90e3)))/(ngsTemp.wavelength/2*pi)); %value??? defocus = z*(1-cos(tel.R/90e3)) [m] 


                    ngsTemp = ngsTemp*zern*telTemp;  

                    % out(:,:,kLayer) = out(:,:,kLayer) + ngsTemp.meanRmPhase;    %phase

                    % Add QKD phase
                    ngsTemp2 = source('wavelength', src(1).photometry);
                    telTemp2 = telescope(tel.D,'resolution',tel.resolution);
                    zern2 = zernike(1, tel.D, 'resolution', tel.resolution);
                    

                    ngsTemp2 = ngsTemp2.*zern2*telTemp2; 

                    

                    % out(:,:,kLayer) = out(:,:,kLayer) + ngsTemp.meanRmPhase + ngsTemp2.meanRmPhase;    %ngsTemp2phase (inmput as arg, remains same for all propagation itr, atm is random)
                    
                    %%%%%% Added Phase %%%%%%
                    out(:,:,kLayer) = out(:,:,kLayer) + ngsTemp.meanRmPhase + QKDPhase;
                    %%%%%%%%%%%%%%%%%%%%%%%%%

                    figure(5)
                    imagesc(angle(out(:,:,kLayer)));
                    title("Phase Value"); colorbar;
                    
                    disp(max(max(ngsTemp.meanRmPhase)));
                    %DEFORMABLE MIRROR ---------------------------------
                    dmPhase = dm.surface*wavenumber;
                    
                    out(:,:,kLayer) = out(:,:,kLayer) - 2*dmPhase;                %phase; 2*dmPhase
                    %----------------------------------------------------
                    
                    %Complex field
                    out(:,:,kLayer) = amplitudeIn.*exp(1i.*out(:,:,kLayer));  %amplitude+phase
                else
%                     % Higher atm layers (interpolation depending on height)
%                     layerR = R_;
%                     u = sampler_m*layerR;
%                     xc = height.*srcDirectionVector1;
%                     yc = height.*srcDirectionVector2;
%                     [xi,yi] = meshgrid(u+xc+m_origin(1),u+yc+m_origin(2));
% 
%                     out_(:,:,kLayer) = linear(xs,ys,phase_m{kLayer},xi,yi);     %phase
                    out_(:,:,kLayer) = phase_m{kLayer};% + QKDPhase; % QKD phase
                    
                    %Fresnel propagation from kLayer-1 to kLayer
                    Height = height - altitude_m(kLayer-1);
                    wave = propTF(src,Height,tel,out(:,:,kLayer-1));         %propTF     %amplitude+phase
                    amplitude = abs(wave);
                    phase = (angle(wave));    %*1/wavenumber???
                    
                    %Add corresponding kLayer phase
                    out(:,:,kLayer) = amplitude.*exp(1i.*(phase+(out_(:,:,kLayer)*(tel.opticalAberration.wavelength/src(1).wavelength))));
                    
                    
                end
                
                
            end
            
            %Add propagation from last atm layer to source height
            %FRESNEL
            outLast = propTF(src,90e3-altitude_m(end),tel,out(:,:,end));    %amplitude+phase
            
            %FRAUNHOFFER
%             outLast = fProp(src,srcHeight-altitude_m(end),tel,out(:,:,end));    %amplitude+phase
            
% % %             outLast = propTF(src,srcHeight-45e3,tel,outLast);    %amplitude+phase
% % % 
% % %             %DEFORMABLE MIRROR ---------------------------------           
% % %             tmp1 = abs(outLast);
% % %             tmp2 = (angle(outLast));
% % % 
% % % 
% % %             dmPhase = (-2*dm.surface*wavenumber);
% % % 
% % % %                     out(:,:,kLayer) = out(:,:,kLayer) + dmPhase;                %phase
% % %        
% % %             %Add mirror focus to 90km (as a fixed value in order to check with other srcHeight the spot will be larger)
% % %             ngsTemp = source;
% % %             telTemp = telescope(tel.D,'resolution',tel.resolution);
% % %             zern = zernike(4, tel.D, 'resolution', tel.resolution);
% % %             zern.c = (-ngsTemp.wavelength/(2*pi))*((90e3*(1-cos(tel.R/90e3)))/(ngsTemp.wavelength/2*pi)); %value??? defocus = z*(1-cos(tel.R/90e3)) [m]
% % %             ngsTemp = ngsTemp*zern*telTemp;
% % % 
% % % %             out(:,:,kLayer) = out(:,:,kLayer) + ngsTemp.meanRmPhase;    %phase
% % % 
% % % 
% % %             outLast =tmp1.*exp(1i.*(tmp2+dmPhase+ngsTemp.meanRmPhase));
% % %             
% %             %----------------------------------------------------
% %             
% % %             out = sum(out,3) + outLast;
% % %             
            out = outLast;

            

            
            
            
            
            
            
        end
        
%%

function[u2,L2] = fProp(src,height,tel,u1)
    % FRAUNHOFFER propagation
    % assumes uniform sampling
    % u1 - source plane field
    % L1 - source plane side length
    % lambda - wavelength
    % z - propagation distance
    % L2 - observation plane side length
    % u2 - observation plane field
    %

    L1 = tel.D;
    lambda = src(1).photometry.wavelength;
    z = height;

    [M,~]=size(u1); %get input field array size
    dx1=L1/M; %source sample interval
    k=2*pi/lambda; %wavenumber
    %
    L2=lambda*z/dx1; %obs sidelength
    dx2=lambda*z/L1; %obs sample interval
    x2=-L2/2:dx2:L2/2-dx2; %obs coords
    [X2,Y2]=meshgrid(x2,x2);
    %
    
    % cc=1/(1i*lambda*z)*exp(1i*k/(2*z)*(X2.^2+Y2.^2));
    % H=fft2(fftshift(cc))*dx1^2; %shift trans func
    % U1=fft2(fftshift(u1)); %shift, fft src field
    % U2=H.*U1; %multiply
    % u2=ifftshift(ifft2(U2)); %inv fft, center obs field
    
    
    cc=1/(1i*lambda*z)*exp(1i*k/(2*z)*(X2.^2+Y2.^2));
    u2=cc.*ifftshift(fft2(fftshift(u1)))*dx1^2;

end

% Q) how is the height of each layer selected?
% does the input side length change for every iteration due to divergence?
% Q) why is Fresnel prop used for propagation when the atmosphere has only
% 3 layers; this makes distance much larger than aperture and hence Fresnel
% approximation des not hold?
 function[u2]=propTF(src,height,tel,u1)
% FRESNEL propagation - transfer function approach
% assumes same x and y side lengths and (i.e. the matrix is square)
% uniform sampling
% u1 - source plane field
% L - source and observation plane side length
% lambda - wavelength
% z - propagation distance
% u2 - observation plane field

L = tel.D;
% L = 4;  % NM - THIS IS WHERE THE PROBLEM IS
lambda = src.photometry.wavelength;
z = height;

[M,N]=size(u1); %get input field array size
dx=L/M; %sample interval

% Q) what is this?
% dx = 0.0078; % NM sampling that works best
% L = dx*M;

k=2*pi/lambda; %wavenumber
fx=-1/(2*dx):1/L:1/(2*dx)-1/L; %freq coords - by doing 1/2dx we are sampling at Nyquist - Nyquist: dx = 1/2f
[FX,FY]=meshgrid(fx,fx);
% H=exp(j*k*z)*exp(-j*pi*lambda*z*(FX.^2+FY.^2)); %trans func
H=exp(-j*pi*lambda*z*(FX.^2+FY.^2)); %trans func  NM this is original

% H =
% (exp(1i*2*pi*(z/lambda))/(1i*lambda*z))*exp(1i*pi*((FX^2+FY^2)/(lambda*z)));
% %impulse response approach
H = fftshift(H); %shift trans func
U1=fft2(fftshift(u1)); %shift, fft src field
U2=H.*U1; %multiply
u2=ifftshift(ifft2(U2)); %inv fft, center obs field


end


   
function F = linear(arg1,arg2,arg3,arg4,arg5)
%LINEAR 2-D bilinear data interpolation.
%   ZI = LINEAR(EXTRAPVAL,X,Y,Z,XI,YI) uses bilinear interpolation to
%   find ZI, the values of the underlying 2-D function in Z at the points
%   in matrices XI and YI.  Matrices X and Y specify the points at which
%   the data Z is given.  X and Y can also be vectors specifying the
%   abscissae for the matrix Z as for MESHGRID. In both cases, X
%   and Y must be equally spaced and monotonic.
%
%   Values of EXTRAPVAL are returned in ZI for values of XI and YI that are
%   outside of the range of X and Y.
%
%   If XI and YI are vectors, LINEAR returns vector ZI containing
%   the interpolated values at the corresponding points (XI,YI).
%
%   ZI = LINEAR(EXTRAPVAL,Z,XI,YI) assumes X = 1:N and Y = 1:M, where
%   [M,N] = SIZE(Z).
%
%   ZI = LINEAR(EXTRAPVAL,Z,NTIMES) returns the matrix Z expanded by
%   interleaving bilinear interpolates between every element, working
%   recursively for NTIMES. LINEAR(EXTRAPVAL,Z) is the same as
%   LINEAR(EXTRAPVAL,Z,1).
%
%   See also INTERP2, CUBIC.

[nrows,ncols] = size(arg3);
%     mx = numel(arg1); my = numel(arg2);
s = 1 + (arg4-arg1(1))/(arg1(end)-arg1(1))*(ncols-1);
t = 1 + (arg5-arg2(1))/(arg2(end)-arg2(1))*(nrows-1);


% Matrix element indexing
ndx = floor(t)+floor(s-1)*nrows;

% Compute intepolation parameters, check for boundary value.
if isempty(s), d = s; else d = find(s==ncols); end
s(:) = (s - floor(s));
if ~isempty(d), s(d) = s(d)+1; ndx(d) = ndx(d)-nrows; end

% Compute intepolation parameters, check for boundary value.
if isempty(t), d = t; else d = find(t==nrows); end
t(:) = (t - floor(t));
if ~isempty(d), t(d) = t(d)+1; ndx(d) = ndx(d)-1; end

% Now interpolate.
onemt = 1-t;
F =  ( arg3(ndx).*(onemt) + arg3(ndx+1).*t ).*(1-s) + ...
    ( arg3(ndx+nrows).*(onemt) + arg3(ndx+(nrows+1)).*t ).*s;


end