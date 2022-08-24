function g = kernelfit_uniq(Nx, Ny, R, targdata, srcdata, lambda)

[Nkx, Nky, Ncoils] = size(targdata);

kypos = 1:R:Nky-R*(Ny-1);

targall = zeros((Nkx-Nx+1)*length(kypos), R-1, Ncoils);
ACSall = zeros((Nkx-Nx+1)*length(kypos), Nx * Ny * Ncoils);

for j = 1:R:Nky-R*(Ny-1) % loop over krenel ky starting locations

    % get the target data points for this line and offset
    targ = targdata((1:Nkx-Nx+1) + (Nx-1)/2,j + R*Ny/2 - (R - 1) - 1 + (1:R-1), 1:Ncoils);

    % reshape for fitting
    targflat = reshape(targ, Nkx-Nx+1, R-1, Ncoils);

    % concatenate into overall target vector
    targall = cat(1, targall, targflat);

    % get the A submatrix for this position
    ACS = zeros(Nx, Ny, Ncoils, Nkx-Nx+1);
    for i = 1:Nkx-Nx+1 % starting readout position in the kernel
        % grab the data around this target point into the ACS matrix
        ACS(:, :, :, i) = srcdata(i:i+Nx-1, j:R:j+R*(Ny-1), :);
    end

    % flatten the submatrix
    ACSflat = permute(reshape(ACS, [Nx * Ny * Ncoils, Nkx-Nx+1]), [2, 1]);

    % concatenate into overall ACS matrix
    ACSall = cat(1, ACSall, ACSflat);

end

AtA = ACSall'*ACSall;
s = svds(AtA);

ACSpinv = inv(ACSall'*ACSall + lambda*s(1)*eye(Nx*Ny*Ncoils)) * (ACSall');
g = reshape(ACSpinv * targall(:, :), [Nx, Ny, Ncoils, R-1, Ncoils]);
