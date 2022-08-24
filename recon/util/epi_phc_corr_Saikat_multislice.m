function [outData, s] = epi_phc_corr_Saikat_multislice(imgData, phcData, echoSeparated)


%%% Dimension 3 is 2 Averages
phcData = squeeze(mean(phcData, 3));

if length(size(phcData)) == 3
    phcData = reshape(phcData, ...
        [size(phcData, 1), 1,size(phcData, 2), size(phcData, 3)]);        
end

phcData = squeeze(sum(phcData, 3));

if length(size(phcData)) == 2
    phcData = reshape(phcData, [size(phcData, 1), 1, size(phcData, 2)]);        
end

[Nro, Nc, nEchoes, Nsl] = size(phcData);
nShots = size(imgData, 2) / nEchoes;

phcData = permute(phcData, [1 3 2 4]); % move coils and slices to the end

if ~echoSeparated
    
    shiftedLines = phcData(:, 2 : 2 : end, :); % collapse slices and columns
    refLines = phcData(:, 1 : 2 : end, :);
    refLines = refLines(:, 1 : size(shiftedLines, 2), :); % remove last line if one too many
    [s, ~] = epiPhsDelayEst_nd3(refLines, shiftedLines, 0.00001, 40);
    s = [zeros(2, nShots - 1) repmat(s, [1 nShots])];
    
    id = isnan(s); s(id) = 0;

    % apply the corrections to with-blip data
    bbCorr = imgData;
    bbCorr = fftshift(ifft(fftshift(bbCorr, 1), [], 1), 1); % ft in read dim to apply phase
    Nro = size(phcData, 1); % # readout points
    A = [2 * pi * ((0 : Nro - 1)' / Nro - 1 / 2), ones(Nro, 1)]; % first and zeroth order phases
    for ii = 1 : nShots * 2 - 1
        bbCorr(:, ii + 1 : 2 * nShots : end, :) = bbCorr(:, ii + 1 : 2 * nShots : end, :) .* exp(1i * A * s(:, ii));
    end
    outData = fftshift(fft(fftshift(bbCorr , 1), [], 1), 1); % undo ft in read dim

else
    
    Nro = size(phcData, 1); % # readout points
    Nechoes = size(phcData, 2); % # of echoes
    A = [2 * pi * ((0 : Nro - 1)' / Nro - 1 / 2) ones(Nro, 1)]; % first and zeroth order phases
    bbCorr = imgData;
    bbCorr = fftshift(ifft(fftshift(bbCorr, 1), [], 1), 1); % ft in read dim to apply phase
    for ii = 1 : size(phcData, 2) - 1
        
        shiftedLines = phcData(:, ii + 1, :); % collapse slices and columns
        refLines = phcData(:, 1, :);
        [s, ~] = epiPhsDelayEst_nd3(refLines, shiftedLines, 0.00001, 40);
        id = isnan(s); s(id) = 0;

        % apply correction to every shot
        bbCorr(:, ii * nShots + 1 : (ii + 1) * nShots, :) = bbCorr(:, ii * nShots + 1 : (ii + 1) * nShots, :) .* ...
            repmat(exp(1i * A * s), [1, nShots, size(bbCorr, 3) * size(bbCorr, 4) * size(bbCorr, 5)]);
    
    end
    outData = fftshift(fft(fftshift(bbCorr , 1), [], 1), 1); % undo ft in read dim
    
end



% example 2x SENSE recon, using only one of the shots

