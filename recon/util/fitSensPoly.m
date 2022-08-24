function sensPoly = fitSensPoly(imgs, nComp, jLambda, genfig)

% probably need to actually multiply the sens into the ssq image,
% compare to the individual coil image!

[Nsq, nCoils, nShots] = size(imgs);
N = sqrt(Nsq);

[kxJ, kyJ] = meshgrid(-(nComp(1) - 1) / 2 : (nComp(1) - 1) / 2, ...
    -(nComp(2) - 1) / 2 : (nComp(2) - 1) / 2);
[x, y] = meshgrid(-N / 2 : N / 2 - 1);
AJ = exp(1i * 2 * pi / N * (x(:) * kxJ(:)' + y(:) * kyJ(:)'));
%wAJ = bsxfun(@times, wts, AJ);
%AJtwAJ = AJ' * wAJ;
%AJtAJ = AJ' * AJ;

imgssq = squeeze(ssq(ssq(imgs, 3), 2));
AJimg = bsxfun(@times, imgssq, AJ);
sc = pinv(AJimg' * AJimg + jLambda * eye(size(AJ, 2))) * ...
    (AJimg' * imgs(:, :));

sensPoly = reshape(AJ * sc, [Nsq, nCoils, nShots]);

if genfig
    figure(100);
    clf;
    im(reshape(squeeze(ssq(sensPoly, 2)), [N N nShots]));
    title('Poly Fit SENSE maps');
    drawnow
end


