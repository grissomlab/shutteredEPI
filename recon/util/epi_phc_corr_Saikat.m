function [outData,s] = epi_phc_corr_Saikat(imgData,phcData)


%%% Dimension 3 is 2 Averages
phcData = squeeze(mean(phcData,3));

if(length(size(phcData)) == 3)
    phcData = reshape(phcData,[size(phcData,1),1,size(phcData,2),size(phcData,3)]);        
end

phcData = squeeze(sum(phcData,3));

if(length(size(phcData)) == 2)
    phcData = reshape(phcData,[size(phcData,1),1,size(phcData,2)]);        
end

[Nro,Nc,nEchoes] = size(phcData);
nShots = size(imgData,2)/nEchoes;

phc_data2 = zeros(Nro,Nc,nShots*nEchoes);
for ii = 1:nShots
    phc_data2(:,:,ii:nShots:end) = phcData;
end


% for s = 1 : nShots
% % interleave even/odd echoes from each shot
% phc_data(:,:,s,s:nShots*2:end) = phcData(:,:,2:2:end);
% phc_data(:,:,s,s+nShots:nShots*2:end) = phcData(:,:,1:2:end);
% end


% estimate corrections from no-blip data

% phc data
% I *think* that for squeeze phc data, the dims are: 
%           ro, coil, shot, even/odd, echo
% so to get it in same format as no-blip (1DFT mode) image data,
% we have to interleave the even and odd echoes for each shot,
% then interleave the shots

% [Nro,Nc,nShots,~,nEchoes] = size(phcData);
% phc_data = zeros(Nro,Nc,nShots,nEchoes);
% 
% phc_data(:,:,:,1:2:end) = phcData(:,:,:,1,1:2:end);
% phc_data(:,:,:,2:2:end) = phcData(:,:,:,2,2:2:end);


% interleave shots (not tested with > 2 shots)
% phc_data2 = zeros(Nro,Nc,nShots*nEchoes);
% for ii = 1:nShots
%     phc_data2(:,:,1+(ii-1):nShots:end) = squeeze(phc_data(:,:,ii,:));
%     %phc_data2(:,:,2:Nshots:end) = squeeze(phc_data(:,:,2,:));
% end


bb = permute(phc_data2,[1 3 2]);

% 1DFT image data
%bb = data_noBlipShutter;


% loop over odd/even pairings
s = zeros(2,2*nShots-1); % the first and zeroth order phase coefficients
for ii = 1:nShots*2-1
    shiftedLines = bb(:,ii+1:2*nShots:end,:); % extract shifted lines
    refLines = bb(:,1:2*nShots:end,:); % extract reference lines
    refLines = refLines(:,1:size(shiftedLines,2),:); % remove last line if one too many
    [s(:,ii),~] = epiPhsDelayEst_nd3(refLines,shiftedLines,0.00001,40);
end
s

id = isnan(s); s(id) = 0;


% apply the corrections to with-blip data
bbCorr = imgData;
bbCorr = fftshift(ifft(fftshift(bbCorr,1),[],1),1); % ft in read dim to apply phase
Nro = size(bb,1); % # readout points
A = [2*pi*((0:Nro-1)'/Nro-1/2) ones(Nro,1)]; % first and zeroth order phases
for ii = 1:nShots*2-1
    bbCorr(:,ii+1:2*nShots:end,:) = bbCorr(:,ii+1:2*nShots:end,:).*repmat(exp(1i*A*s(:,ii)),...
        [1 size(bbCorr(:,ii+1:2*nShots:end,:),2) size(bbCorr,3)]);
end
outData = fftshift(fft(fftshift(bbCorr,1),[],1),1); % undo ft in read dim
%outData = fftshift(ifft(fftshift(bbCorr,2),[],2),2); % ft in PE dim to get 2D image
%bbssq = sqrt(sum(abs(bbCorr).^2,3)); % sum-of-squares image


% example 2x SENSE recon, using only one of the shots

