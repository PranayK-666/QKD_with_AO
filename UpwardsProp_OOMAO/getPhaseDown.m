function phaseDown = getPhaseDown(launch,telRx,phase)
   
    %Select central area of lgsDown.phase as a perfect measurement for the
    %DM coefs calculation.

    %launch atm   
    atmLaunch             = launch.opticalAberration;
    layersLaunch          = atmLaunch.layer;
    R_Launch              = launch.R;
    R_Rx                  = telRx.R;
    layersNPixelLaunch    = [ layersLaunch.nPixel ] ;
    
    phaseIn = phase;

    [mIn,nIn] = size(phaseIn);
    ratio = round(R_Rx/R_Launch);
    mOut = mIn/ratio;
    nOut = nIn/ratio;
    
    if ratio == 1
        phaseOutTmp = phaseIn;
    else
        phaseOutTmp = phaseIn(mIn/2-mOut/2:mIn/2+mOut/2,nIn/2-nOut/2:nIn/2+nOut/2);
    end
    
    % phaseOutTmp = phaseIn;
        
    phaseDown = imresize(phaseOutTmp,[layersNPixelLaunch(1) layersNPixelLaunch(1)]);
        
 
end