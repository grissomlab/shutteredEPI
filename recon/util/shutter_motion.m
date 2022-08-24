function [img, params, rotDict] = shutter_motion(imMoved, imRef, ...
    angles, shifts, rotDict)

Nx = size(imMoved, 1);
Ny = size(imMoved, 2);

nShifts = length(shifts);
nAngles = length(angles);

imMovedSave = imMoved;

% define k-grid for rotation and translation
[kx, ky] = ndgrid(-1 / 2 : 1 / Nx : 1 / 2 - 1 / Nx, ...
    -1 / 2 : 1 / Ny : 1 / 2 - 1 / Ny);

% filter both the moved and reference images
% use the cone filter from Eddy et al, MRM 36:923-931 (1996)
% b = 1; % pixels, radius of filter
% a = sqrt(kx.^2 + ky.^2) * sqrt(Nx * Ny);
% w = 1 - b * a;
% w(w < 0) = 0;
% imMoved = conv2(imMoved, w, 'same');
% imRef = conv2(imRef, w, 'same');

% calculate a rotated reference image for each angle
if isempty(rotDict)
    
    rotDict = zeros(Nx, Ny, nAngles);
    
    parfor ii = 1 : nAngles
        
        %%%% Rotation Grot fatrix
        kxr = kx * cosd(angles(ii)) + ky * sind(angles(ii));
        kyr = -kx * sind(angles(ii)) + ky * cosd(angles(ii));
        Gtheta = Gmri([kxr(:) kyr(:)], true([Nx Ny]), 'fov', [Nx Ny]);
        
        rotDict(:, :, ii) = reshape(Gtheta * abs(imRef(:)), [Nx, Ny]);
        %imRotDict(:, :, ii) = ift2(rotDict(:, :, ii));

    end

end

% take the moved image to k-space
kMoved = ft2(abs(imMoved));

% loop over angles and calculate correlations with input image k-space
corrAngle = zeros(nAngles, 1);
for ii = 1 : nAngles
    
    corrAngle(ii) = sum(col(abs(rotDict(:, :, ii)) .* abs(kMoved)));
    %corrAngle(ii) = sum(col(abs(imRotDict(cols, :, ii)) .* abs(imMoved(cols, :))));
    
end
[~, bestAngleInd] = max(corrAngle);
params.rotation = angles(bestAngleInd);
%figure;plot(corrAngle)

% fix the rotation in the filtered image
kxr = kx * cosd(-params.rotation) + ky * sind(-params.rotation);
kyr = -kx * sind(-params.rotation) + ky * cosd(-params.rotation);
Gtheta = Gmri([kxr(:) kyr(:)], true([Nx Ny]), 'fov', [Nx Ny]);
kReg = reshape(Gtheta * abs(imMoved(:)), [Nx Ny]);

% loop over translations and find highest correlation
innProd = col(conj(ft2(abs(imRef))) .* kReg);
corrShift = zeros(nShifts, nShifts);
parfor ii = 1 : nShifts % x-shift
    
    for jj =  1 : nShifts % y-shift

        phs = exp(1i * 2 * pi * ...
            (shifts(ii) * kx(:) + shifts(jj) * ky(:)));
        corrShift(ii, jj) = phs.' * innProd;
        
    end
    
end
%figure;im(abs(corrShift))
[~, bestShiftInd] = max(abs(corrShift(:)));
[xShiftInd, yShiftInd] = ind2sub([nShifts, nShifts], bestShiftInd);
params.dx = shifts(xShiftInd);
params.dy = shifts(yShiftInd);

% fix the rotation and translation in the original moved image
kReg = reshape(Gtheta * imMovedSave(:), [Nx Ny]);
img = ift2(kReg .* exp(1i * 2 * pi * (params.dx * kxr + params.dy * kyr)));

