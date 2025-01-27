function atmOut = getLaunchAtm(launch,telRx)
    %Get launch atm by selecting central area of telRx atm in order to be
    %sure of using same atmosphere when launching and receiving.
    
    %telRx atm
    atmRx                = telRx.opticalAberration;
    nLayerRx             = atmRx.nLayer;
    layersRx             = atmRx.layer;
    altitude_m           = [layersRx.altitude];
    sampler_m            = linspace(-1,1,telRx.resolution);
    phaseRx              = { layersRx(:).phase }; %{}
    R_Rx                 = telRx.R;
%     layerSamplingRx      = { layersRx.sampling };
    layersNPixelRx       = [ layersRx.nPixel ] ;
    
    %launch atm   
    atmLaunch             = launch.opticalAberration;
    nLayerLaunch          = atmRx.nLayer;
    layersLaunch          = atmRx.layer;
    altitude_m            = [layersLaunch.altitude];
    R_Launch              = launch.R;
    layersNPixelLaunch    = [ layersLaunch.nPixel ] ;
    
    phaseOut = cell(1,nLayerLaunch);
    for l=1:nLayerRx
        
        phaseIn = phaseRx{l};
        
        [mIn,nIn] = size(phaseIn);
        ratio = round(R_Rx/R_Launch);
        mOut = mIn/ratio;
        nOut = nIn/ratio;
        if ratio == 1
            phaseOutTmp = phaseIn;
        else
            phaseOutTmp = phaseIn(mIn/2-mOut/2:mIn/2+mOut/2,nIn/2-nOut/2:nIn/2+nOut/2);
        end
        
        phaseOut{l} = imresize(phaseOutTmp,[layersNPixelLaunch(l) layersNPixelLaunch(l)]);
        
    end
    
    atmOut.nLayer = nLayerLaunch;
    atmOut.layer.phase = phaseOut;
    atmOut.layer.altitude = altitude_m;
    atmOut.layer.nPixel = layersNPixelLaunch;
    atmOut.wavelength = atmLaunch.wavelength;
    atmOut.R = R_Launch;



end