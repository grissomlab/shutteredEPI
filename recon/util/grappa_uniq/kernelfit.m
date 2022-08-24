function g = kernelfit(Nx, Ny, R, data)

[Nkx, Nky, Ncoils] = size(data);

% get the target data points
targ = zeros(Nkx-Nx+1, Nky-R*(Ny-1), R-1, Ncoils);
for r = 1:R-1
    % grab these data point into the target matrix
    targ(:, :, r, :) = data((1:Nkx-Nx+1) + (Nx-1)/2, (1:Nky-R*(Ny-1)) + R*Ny/2 - (R - 1) - 1 + r, 1:Ncoils);
end

% reshape for fitting
targflat = reshape(targ, [(Nkx-Nx+1)*(Nky-R*(Ny-1)),R-1,Ncoils]);

% compute ACS pseudoinverse
ACSpinv = pinv(buildACS(Nx, Ny, R, data));

% Fit the kernels
g = reshape(ACSpinv * targflat(:, :), [Nx, Ny, Ncoils, R-1, Ncoils]);

