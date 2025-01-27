function atmOut = getLaunchAtm_FOV(launch,telRx)
    %Get launch atm by selecting central area of telRx atm in order to be
    %sure of using same atmosphere when launching and receiving.
    %Use this function when telRx has FOV
    
    %telRx atm
    atmRx                = telRx.opticalAberration;
    nLayerRx             = atmRx.nLayer;
    layersRx             = atmRx.layer;
    altitude_m           = [layersRx.altitude];
    sampler_m            = linspace(-1,1,telRx.resolution);
    phaseRx              = { layersRx.phase };
    R_Rx                 = telRx.R;
    layerSamplingRx      = { layersRx.sampling };
    layersNPixelRx       = [ layersRx.nPixel ] ;
    
    %launch atm   
    atmLaunch             = launch.opticalAberration;
    nLayerLaunch          = atmLaunch.nLayer;
    layersLaunch          = atmLaunch.layer;
    altitude_m            = [layersLaunch.altitude];
    R_Launch              = launch.R;
    layersNPixelLaunch    = [ layersLaunch.nPixel ] ;
    
    pxLayer = layersNPixelLaunch(1);
    
    phaseOut = cell(1,nLayerLaunch);
    for l=1:nLayerRx
        
        phaseIn = phaseRx{l};
        pxLayerIn = layersNPixelRx(l);
        phaseOut{l} = phaseIn(pxLayerIn/2+1-pxLayer/2:pxLayerIn/2+pxLayer/2,pxLayerIn/2+1-pxLayer/2:pxLayerIn/2+pxLayer/2);

               
       
    end
    
    atmOut.nLayer = nLayerLaunch;
    atmOut.layer.phase = phaseOut;
    atmOut.layer.altitude = altitude_m;
    atmOut.layer.nPixel = layersNPixelLaunch;
    atmOut.wavelength = atmLaunch.wavelength;
    atmOut.R = R_Launch;



end