function g = kernelfit_uniq(Nx, Ny, R, targdata, srcdata, vt, Nv, lambda)

[Nkx, Nky, Ncoils] = size(targdata);

kypos = 1:R:Nky-R*(Ny-1);

g = zeros(Nx, Ny, Ncoils, R-1, Ncoils, length(kypos));

% get the sub-v matrix
vsub = vt(1:Nv, :)';

for j = 1:R:Nky-R*(Ny-1) % loop over krenel ky starting locations

    % get the target data points for this line and offset
    targ = targdata((1:Nkx-Nx+1) + (Nx-1)/2,j + R*Ny/2 - (R - 1) - 1 + (1:R-1), 1:Ncoils);

    % reshape for fitting
    targflat = reshape(targ, Nkx-Nx+1, R-1, Ncoils);

    % get the A submatrix for this position
    ACS = zeros(Nx, Ny, Ncoils, Nkx-Nx+1);
    for i = 1:Nkx-Nx+1 % starting readout position in the kernel
        % grab the data around this target point into the ACS matrix
        ACS(:, :, :, i) = srcdata(i:i+Nx-1, j:R:j+R*(Ny-1), :);
    end

    % flatten the submatrix
    ACSflat = permute(reshape(ACS, [Nx * Ny * Ncoils, Nkx-Nx+1]), [2, 1]);

    % multiply them
    Avsub = ACSflat * vsub;
    if lambda > 0
        AvtAv = Avsub' * Avsub;
        evals = eig(AvtAv);
        lambdae = lambda*min(evals);
        AvsubInv = inv(AvtAv + lambdae .* eye(size(AvtAv, 1))) * Avsub';
    else
        AvsubInv = pinv(Avsub);
    end

    % get the coefficients and kernels
    gjind = find(kypos == j);
 
    coeffs = reshape(AvsubInv * targflat(:, :), [Nv, R-1, Ncoils]);
    g(:,:,:,:,:,gjind) = reshape(vsub * coeffs(:, :), [Nx, Ny, Ncoils, R-1, Ncoils]);

end
