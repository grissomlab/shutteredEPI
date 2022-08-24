function G = Gmri_motion_object(N, Rotation, Translation, Type)

if nargin < 4
    Type = 1; % old type
end

%%%% Simulate Translation By FFT

[kx, ky] = ndgrid(-1 / 2 : 1 / N : 1 / 2 - 1 / N, ...
    -1 / 2 : 1 / N : 1 / 2 - 1 / N);

%%%% Rotation Grot fatrix
kxr = kx * cosd(Rotation) + ky * sind(Rotation);
kyr = -kx * sind(Rotation) + ky * cosd(Rotation);
Gtheta = Gmri([kxr(:) kyr(:)], true([N N]), 'fov', [N N]);


%%%% FFT Gcart fatrix
%Gcart = Gmri([kx(:) ky(:)], true([N N]), 'fov', [N N]); % REDFLAG - replace with Gmri_epi
Gcart = Gmri_epi(true([N N]), true([N N]), 0, 0, '2D', true);

if Type == 1 % Type 1, consitent with vuRigidRegistration

    %%%% Translation Matrix
    Translation_Matrix = diag_sp(1 / (N * N)^(2) * exp(-1i * 2 * pi * ...
        (-Translation(2) * kx(:) + -Translation(1) * ky(:)))); %%Original
    
    %      GFinal = block_fatrix({Gcart',Translation_Matrix,Gtheta},'type','mult');
    G = block_fatrix({Gcart', Gtheta, Gcart', Translation_Matrix, Gcart}, 'type', 'mult');

else % Type 2, consistent with shutter_motion
    
    %%%% Translation Matrix
    Translation_Matrix = diag_sp(1 / (N * N) * exp(-1i * 2 * pi * ...
        (Translation(1) * kxr(:) + Translation(2) * kyr(:)))); %%Original
    
    G = block_fatrix({Gcart', Translation_Matrix, Gtheta}, 'type', 'mult');
    
end
    
