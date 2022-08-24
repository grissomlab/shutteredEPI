function output = ft2(in)
% Performs fftshift(ifft2(fftshift( input)))
% 2D inverse FT
dims = size(in);
in = in(:, :, :); % collapse non-space dims
output = zeros(dims);
for ii = 1 : size(in, 3)
    output(:, :, ii) = fftshift(fft2(fftshift(in(:, :, ii))));
end
