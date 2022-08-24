function recondata = predict_uniq(recondata, Nx, Ny, R, g)

[Nkx, Nky, Ncoils] = size(recondata);

kypos = 1:R:Nky-R*(Ny-1);

gsz = size(g);

% flatten the kernel
g = reshape(permute(g, [4,5,6,1,2,3]), [gsz(4), gsz(5), gsz(6), gsz(1)*gsz(2)*gsz(3)]);
g = permute(g, [4, 1, 2, 3]);
% also flatten R and coil dimensions of the kernel
g = reshape(permute(g, [1, 4, 2, 3]), [gsz(1)*gsz(2)*gsz(3), gsz(6), gsz(4)*gsz(5)]);
g = permute(g, [3, 1, 2]);

% apply grappa prediction
% for j = 1:R:Nky-R*(Ny-1) % the kernel's starting phase encode position
%     gjind = find(kypos == j);
%     gt = g(:,:,gjind);
%     predictions = zeros(Nkx-Nx+1, R-1, Ncoils);
%     parfor i = 1:Nkx-Nx+1 % the kernel's starting readout position
%         srcdata = recondata(i:i+Nx-1, j:R:j+R*(Ny-1), :);
%         predictions(i, :, :) = reshape(gt * srcdata(:), [R-1, Ncoils]);
%     end
%     recondata((1:Nkx-Nx+1) + (Nx-1)/2, j + R*Ny/2 - (R - 1) - 1 + (1:R-1), :) = predictions;
% end

for j = 1:R:Nky-R*(Ny-1) % the kernel's starting phase encode position
    gjind = find(kypos == j);
    gt = g(:,:,gjind);
    for i = 1:Nkx-Nx+1 % the kernel's starting readout position
        srcdata = recondata(i:i+Nx-1, j:R:j+R*(Ny-1), :);
        recondata(i + (Nx-1)/2, j + R*Ny/2 - (R - 1) - 1 + (1:R-1), :) = reshape(gt * srcdata(:), [R-1, Ncoils]);
    end
end