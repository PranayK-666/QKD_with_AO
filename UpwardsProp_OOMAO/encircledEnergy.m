function EE = encircledEnergy(frame, radius, method)
% Encircled energy on a cirlce of size "circle_size" in terms of % of
% total energy
% "radius" in pixel
% D_sat_px = number of px that represent the satellite receiver
% Method = 'CoG' uses the centre of gravity of the frame as reference
% point for the EE calculation;
% Method = 'centralPx' uses the frame central pixel as reference point
% for the EE calculation.

total_energy = sum(frame(:));

if strcmp(method, 'CoG') == 1
  

    %Center of mass
    %These next 4 lines produce a matrix C whose rows are
    % the pixel coordinates of the image.

    C = cellfun(@(n) 1:n, num2cell(size(frame)),'uniformoutput',0);
    [C{:}] = ndgrid(C{:});
    C = cellfun(@(x) x(:), C,'uniformoutput',0);
    C = [C{:}];

    %This line computes a weighted average of all the pixel coordinates.
    %The weight is proportional to the pixel value.

    A = frame(:);
    CenterOfMass = A.'*C/sum(A,'double');

    %Encircled energy
    centerX = CenterOfMass(1);
    centerY = CenterOfMass(2);

else if strcmp(method, 'centralPx') == 1
        centerX = size(frame,1)/2;
        centerY = size(frame,2)/2;
end

% Reference to calculate EE is the central px of the frame.

dimy = size(frame,1);
dimx = size(frame,2);
circle = zeros(dimy, dimx);

for i = 1:dimy
    for j = 1:dimx
        dist = sqrt((j-centerX)^2+(i-centerY)^2);
        if dist <= radius
            circle(i,j) = 1;
        end
    end
end

EE = (sum(frame(find(circle==1)))/total_energy)*100;



end