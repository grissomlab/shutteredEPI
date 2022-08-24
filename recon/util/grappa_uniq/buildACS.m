function ACSflat = buildACS(Nx, Ny, R, data)

[Nkx, Nky, Ncoils] = size(data);

% Build GRAPPA ACS data matrix
ACS = zeros(Nx, Ny, Ncoils, Nkx-Nx+1, Nky-R*(Ny-1));
for i = 1:Nkx-Nx+1 % starting readout position in the kernel
    for j = 1:Nky-R*(Ny-1) % starting phase encode position in the kernel
        % grab the data around this target point into the ACS matrix
        ACS(:, :, :, i, j) = data(i:i+Nx-1,j:R:j+R*(Ny-1),:);
    end
end

ACSflat = permute(reshape(ACS, [Nx * Ny * Ncoils, (Nkx-Nx+1)*(Nky-R*(Ny-1))]), [2, 1]);

