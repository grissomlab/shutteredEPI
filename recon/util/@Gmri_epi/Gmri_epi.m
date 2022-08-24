function ob = Gmri_epi(kmask,immask,dk,dphi,oneDorTwoD,shotMask,sensMaps)

% kmask: mask of sampled k-space locations
if ~islogical(kmask)
    error 'kmask must be a logical nd matrix'
end
ob.kmask = kmask;
if ~islogical(immask)
    error 'immask must be a logical nd matrix'
end
ob.immask = immask;
if numel(immask) ~= numel(kmask)
    error 'size of immask must = size of kmask'
end
ob.dk = dk; % Frequency-dim delay for each shot
ob.dphi = dphi; % Phase offset for each shot
if length(dk) ~= length(dphi)
    error 'length of delays must equal length of phases'
end
ob.Nsh = (length(dk)+1)/2; % # shots is length of delays plus one, divided by two
if nargin < 5 % no 1D or 2D FFT's - we always have to do a DFT in the PE dim, 
    % but we can keep RO dim in the spatial domain. 
    oneDorTwoD = '2D'; % default is 2D FFT
end
ob.oneDorTwoD = oneDorTwoD;
if nargin < 6 % no shot mask provided
    shotMask = true(ob.Nsh,1); 
end
ob.shotMask = shotMask;
if nargin < 7 % no sens map provided
    sensMaps = 1;
end
ob.sensMaps = sensMaps;
[ob.M,ob.N] = size(kmask);
ob.C = conj(dftmtx(2*ob.Nsh)); % segment combination weights
% Phase encode-dependent phase terms
ob.pePhs = exp(-1i*2*pi*(0:2*ob.Nsh-1)'*(0:ob.N/(2*ob.Nsh)-1)./ob.N...
    +1i*[0;ob.dphi(:)]*ones(1,ob.N/(2*ob.Nsh)));
% Frequency encode-dependent phase terms (i.e. the delays in the FE dim)
if streq(ob.oneDorTwoD,'2D')
    ob.fePhs = exp(-1i*2*pi*fftshift(-ob.M/2:ob.M/2-1)'*[0;ob.dk(:)]'./ob.M);
else % 1D
    ob.fePhs = exp(-1i*2*pi*(-ob.M/2:ob.M/2-1)'*[0;ob.dk(:)]'./ob.M);
end
ob.allPhs = zeros(ob.M,ob.N/(2*ob.Nsh),2*ob.Nsh);
for ii = 1:2*ob.Nsh
    ob.allPhs(:,:,ii) = ob.fePhs(:,ii)*ob.pePhs(ii,:);
end
ob.adjoint = 0;
ob = class(ob,'Gmri_epi');
