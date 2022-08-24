function recondata = predict(recondata, Nx, Ny, R, g)

[Nkx, Nky, Ncoils] = size(recondata);

gsz = size(g);

% flatten the kernel
g = reshape(permute(g, [4,5,1,2,3]), [gsz(4)*gsz(5), gsz(1)*gsz(2)*gsz(3)]);

% apply grappa prediction
for j = 1:R:Nky-R*(Ny-1) % the kernel's starting phase encode position
    for i = 1:Nkx-Nx+1 % the kernel's starting readout position
        srcdata = recondata(i:i+Nx-1, j:R:j+R*(Ny-1), :);
        recondata(i + (Nx-1)/2, j + R*Ny/2 - (R - 1) - 1 + (1:R-1), :) = ... 
            reshape(g * srcdata(:), [R-1, Ncoils]);
    end
end

