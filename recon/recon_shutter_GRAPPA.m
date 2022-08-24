%%%%%%%%%%% Shutter EPI Recon %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath util
addpath util/grappa_uniq

% recon parameters
refSaveData = false;
accSaveData = false;
refcorrect = true;
echoByEchoPHCCorr = false;
separateSlicePHCCorr = false;
moco = false 
motionAngles = -10 : 0.125 : 10;
motionShifts = -10 : 0.125 : 10;
phsCorr = 'after' %'after'; % estimate and perform phase correction 'after' GRAPPA, or 'none'
phsN = 8;

% Data selection
%dataset = '20190117_phantomMotion'; % 
%dataset = '20171102_Datalist_Human'; % original single-slice data
%dataset = 'datalist_20200902 copy pos 1'; % two-position 25-slice brain data
%dataset = 'datalist_20200902 FullEX pos 1'; % two-position 25-slice FULLEX brain data
%dataset = 'datalist_20190122_pos_1'; % two-position single-slice brain data
%dataset = 'datalist_20210511 FullEX pos 3';
%fNameSave = 'FullEX_20210511_pos_3_1111_nomoco_conv'
%dataset = 'datalist_20210511 pos 3'; % three-position brain data with 
%fNameSave = 'Shutter_20210511_pos_3_1111_nomoco_conv'
%swapRefDataShotsCoils = true;

%% 6/24 FullEX
%dataset = 'datalist_20210624 FullEX Series 13 Breathhold';
%fNameSave = 'FullEX_Series_13_Breathhold_moco_chi10ov100_phsN8';
%dataset = 'datalist_20210624 FullEX Series 15 Visual';
%fNameSave = 'FullEX_Series_15_Visual_moco_chi10ov100_phsN8';
%dataset = 'datalist_20210624 FullEX Series 21 Breathhold'; % Jonathan 
%fNameSave = 'FullEX_Series_21_Breathhold_moco_chi10ov100_phsN8';   
%dataset = 'datalist_20210624 FullEX Series 23 Visual'; % Jonathan 
%fNameSave = 'FullEX_Series_23_Visual_nophscorr_moco_chi10ov100_slice4'   

%% 6/24 Shuttered
%dataset = 'datalist_20210624 Shutter Series 9 Breathhold';
%fNameSave = 'Shutter_Series_9_Breathhold_moco_chi10ov100_phsN8'
%dataset = 'datalist_20210624 Shutter Series 17 Breathhold';
%fNameSave = 'Shutter_Series_17_Breathhold_moco_chi10ov100_phsN8'
%dataset = 'datalist_20210624 Shutter Series 19 Visual';
%fNameSave = 'Shutter_Series_19_Visual_nophscorr_moco_chi10ov100_slice4'
%swapRefDataShotsCoils = true;

%% 11/08 FullEX
%dataset = '20211108_Anderson_344229 FullEX Series 17 tSNR'; % Jonathan 
%fNameSave = '20211108_Anderson_344229_FullEX_Series_17_tSNR_norefcorr_test'   
%dataset = '20211108_Anderson_344229 FullEX Series 21 Visual'; % Jonathan 
%fNameSave = '20211108_Anderson_344229_FullEX_Series_21_Visual_moco_phscorr_conv'   

%dataset = '20211108_Anderson_344232 FullEX Series 9 tSNR'; % Abitha
%fNameSave = '20211108_Anderson_344232_FullEX_Series_9_tSNR_moco_phscorr_conv'
%dataset = '20211108_Anderson_344232 FullEX Series 11 Visual'; % Abitha
%fNameSave = '20211108_Anderson_344232_FullEX_Series_11_Visual_moco_phscorr_conv'

%% 11/08 Shuttered
dataset = '20211108_Anderson_344229 Shutter Series 24 tSNR'; % Jonathan
%fNameSave = '20211108_Anderson_344229_Shutter_Series_24_tSNR_moco_phscorr_conv'
%dataset = '20211108_Anderson_344229 Shutter Series 26 Visual'; % Jonathan
%fNameSave = '20211108_Anderson_344229_Shutter_Series_26_Visual_moco_phscorr_conv'

%dataset = '20211108_Anderson_344232 Shutter Series 14 tSNR' % Abitha
%fNameSave = '20211108_Anderson_344232_Shutter_Series_14_tSNR_moco_phscorr_conv_refcorr'
%dataset = '20211108_Anderson_344232 Shutter Series 16 Visual' % Abitha
%fNameSave = '20211108_Anderson_344232_Shutter_Series_16_Visual_moco_phscorr_conv'


%% 11/16 FullEX
%dataset = '20211116_Anderson_344321 FullEX Series 13 tSNR'
%fNameSave = '20211116_Anderson_344321_FullEX_Series_13_tSNR_moco_phscorr_conv'
%dataset = '20211116_Anderson_344321 FullEX Series 15 Visual'
%fNameSave = '20211116_Anderson_344321_FullEX_Series_15_Visual_moco_phscorr_conv'

%% 11/16 Shuttered
%dataset = '20211116_Anderson_344321 Shutter Series 18 tSNR'
%fNameSave = '20211116_Anderson_344321_Shutter_Series_18_tSNR_moco_phscorr_conv_refcorr'
%dataset = '20211116_Anderson_344321 Shutter Series 20 Visual'
%fNameSave = '20211116_Anderson_344321_Shutter_Series_20_Visual_moco_phscorr_conv_refcorr'


retrospective = false; 
% which slices to reconstruct?
%slices = 1 : 15; % motion cases
%slices = 1 : 7; % fMRI
slices = 4;
nSlices = length(slices);
% the following are only relevant if retrospective == false
accDataInd = [1]; % which accelerated data set(s) to use, within this data
shotSelect = [1 1 1 1]; % which index into accDataInd to take each shot from

[fnames, path_name, params] = shutterData(dataset, machine);
nCoils = params.nCoils;
nShots = params.nShots;
nSlicesData = params.nSlices;
sepRefScans = params.sepRefScans;
multislice = nSlicesData > 1;
R = params.R;
Nro = params.Nro;
Npe = params.Npe;
fullEXAcq = params.fullEXAcq;
if ~fullEXAcq 
    swapRefDataShotsCoils = params.swapRefDataShotsCoils;
end

% GRAPPA parameters
grappap.Nx = 3;     % kernel extent in frequency-encoded direction 
grappap.Ny = 4;     % kernel extent in primary phase-encoded direction   (a.k.a., Nblocks1)
grappap.code = 'convacclines'; % 'Jon', 'conv', 'uniq', or 'convacclines'
grappap.ndyn_fit = 10; % number of dynamic images to average for use in fitting the kernel
grappap.lambda_all_convacc = 10^-3;

%%%% Read Reference Shutter Data ----------------------------------------/
disp 'Reading Reference Shutter Data'
if ~fullEXAcq

    if sepRefScans % Each shot's reference scan was collected separately, in different raw data files
        data_Ref = zeros([N, N, nSlices, nCoils, nShots]);
        for n = 1 : nShots
            
            % read in and phase+delay-correct
            data_Ref(:, :, :, :, n) = readPhsCorrCropData( ...
                path_name, fnames.shutterRef{n}, refSaveData, multislice, ...
                echoByEchoPHCCorr, separateSlicePHCCorr, slices);

        end
    else % All reference data are in a single raw data file

        data_Ref = readPhsCorrCropData( ...
                path_name, fnames.shutterRef{1}, refSaveData, multislice, ...
                echoByEchoPHCCorr, separateSlicePHCCorr, slices);
        if swapRefDataShotsCoils
            data_Ref = reshape(data_Ref, [Nro, Npe * nShots * R, nSlices, nCoils, nShots]);
        else
            data_Ref = reshape(data_Ref, [Nro, Npe * nShots * R, nSlices, nShots, nCoils]);
            data_Ref = permute(data_Ref, [1 2 3 5 4]);
        end
        % we have to time-reverse it!
        data_Ref = data_Ref(:, :, :, :, nShots : -1 : 1); 

        if refcorrect == true 
            % correct the relative amplitudes and phases of each shot
            % this is based on the idea that if the shots are not properly equalized,
            % then the ACS matrix will have a high rank because there will not be 
            % a single shift-invariant kernel that works across all shots.
            % so if we force the ACS matrix to have lower rank then that should suppress 
            % the aliases somewhat and we can adjust the shot amplitudes and phases to be 
            % consistent with those lower aliases.
            %data_Ref_corr = zeros(size(data_Ref));
            for sl = 1 : nSlices 
                for nn = 1 : nShots 
                    disp(['Correcting reference data for shot ' num2str(nn) '/' num2str(nShots) ' of slice ' num2str(sl) '/' num2str(nSlices)]);
                    data_corr = sqz(data_Ref(:,:,sl,:,nn));
                    data_orig = sqz(data_Ref(:,:,sl,:,nn));
                    for kk = 1 : 5 % iterations
                        A = im2row(data_corr, [6 6]);
                        A = reshape(A, [size(A, 1), 6*6*nCoils]);
                        [u, s, v] = svd(A, 'econ');
                        s = max(0, s - s(3*nCoils, 3*nCoils));
                        At = u * s * v';
                        datat = row2im(reshape(At, [size(A, 1), 6*6, nCoils]), ...
                            [size(data_Ref, 1) size(data_Ref, 2)], [6, 6]);
                        % normalize overall scaling
                        datat = datat * (datat(:) \ data_orig(:)); 
                        % loop over shots and correct amplitude and phase shifts
                        for ii = 1 : R*nShots
                            data_corr(:, ii:R*nShots:end, :) = data_corr(:, ii:R*nShots:end, :) * ...
                                (col(data_corr(:, ii:R*nShots:end, :)) \ col(datat(:, ii:R*nShots:end, :)));
                        end
                    end
                    data_Ref(:, :, sl, :, nn) = data_corr;
                end
            end 
            %data_Ref = data_Ref_corr;
        %else
        %    load foo
        end 

    end

else % FullEX
    
    data_Ref = readPhsCorrCropData( ...
        path_name, fnames.shutterRef{1}, refSaveData, multislice, ...
        echoByEchoPHCCorr, separateSlicePHCCorr, slices);
    
    if refcorrect == true 
        % correct the relative amplitudes and phases of each shot
        % this is based on the idea that if the shots are not properly equalized,
        % then the ACS matrix will have a high rank because there will not be 
        % a single shift-invariant kernel that works across all shots.
        % so if we force the ACS matrix to have lower rank then that should suppress 
        % the aliases somewhat and we can adjust the shot amplitudes and phases to be 
        % consistent with those lower aliases.
        for sl = 1 : nSlices 
            disp(['Correcting reference data for slice ' num2str(sl) '/' num2str(nSlices)]);
            data_corr = sqz(data_Ref(:,:,sl,:));
            data_orig = sqz(data_Ref(:,:,sl,:));
            for kk = 1 : 5 % iterations
                A = im2row(data_corr, [6 6]);
                A = reshape(A, [size(A, 1), 6*6*nCoils]);
                [u, s, v] = svd(A, 'econ');
                s = max(0, s - s(3*nCoils, 3*nCoils));
                At = u * s * v';
                datat = row2im(reshape(At, [size(A, 1), 6*6, nCoils]), ...
                    [size(data_Ref, 1) size(data_Ref, 2)], [6, 6]);
                % normalize overall scaling
                datat = datat * (datat(:) \ data_orig(:)); 
                % loop over shots and correct amplitude and phase shifts
                for ii = 1 : R*nShots
                    data_corr(:, ii:R*nShots:end, :) = data_corr(:, ii:R*nShots:end, :) * ...
                        (col(data_corr(:, ii:R*nShots:end, :)) \ col(datat(:, ii:R*nShots:end, :)));
                end
            end
            data_Ref(:, :, sl, :) = data_corr;
        end 
    end  

end
img_Ref = ift2(data_Ref);
disp 'Combining Reference Shutter Images'
img_Ref_comb = ssq(img_Ref(:, :, :, :), 4); % ssq across shots & coils  

if fullEXAcq % duplicate reference to all shots
    data_Ref = repmat(data_Ref, [1 1 1 1 nShots]);
end   

%%%% Read Accelerated Shutter Ex Data ---------/
disp 'Reading Accelerated Shutter Data'

% decimate the Ref data to the same R as the prospective data, 
% to fill in data_Acc_ref, for phase correction or retrospective recon
nCoils = size(data_Ref, 4);
data_Acc_ref = zeros([Nro, nSlices, Npe * nShots, nCoils]);
for n = 1 : nShots 

    % interleave the undersampled Ref data
    start_ky = (n - 1) * R + 1;
    data_Acc_ref(:, :, n : nShots : end, :) = ...
        permute(data_Ref(:, start_ky : R * nShots : end, :, :, n), [1 3 2 4]);

end

if ~retrospective

    for ii = 1 : length(accDataInd)
        data_Acc_tmp = readPhsCorrCropData(path_name, fnames.shutterAcc{accDataInd(ii)}, ...
            accSaveData, multislice, echoByEchoPHCCorr, separateSlicePHCCorr, slices);
        if ~multislice
            % add a dimension where the slices would be
            data_Acc_tmp = permute(data_Acc_tmp, [1 4 2 3 5]);
        else
            % permute so slices is third dim
            data_Acc_tmp = permute(data_Acc_tmp, [1 3 2 4 5]);
        end
        if ii == 1
            % start by taking all data from first file
            data_Acc = data_Acc_tmp;
        else
            % fill in the data array with this file's data,
            % for each desired shot. There must be a simpler way to do this...
            for jj = 1 : nShots
                if shotSelect(jj) == ii % pull this shot data from this file
                    data_Acc(:, :, jj : nShots : end, :, :) = ...
                        data_Acc_tmp(:, :, jj : nShots : end, :, :);
                end
            end
        end
    end

else

    data_Acc = data_Acc_ref;

end

nDyns = size(data_Acc, 5); % number of time points 


%-------------------------------------------------------------------------%
%%%% RECONSTRUCTION  %%%%%%%%%%%%%%%%%%%%%

% make a window for low-res phase correction
window = padarray(col(tukeywin(phsN)) * col(tukeywin(phsN))', ...
    [Nro - phsN, Npe * nShots * R - phsN] ./ 2); 

imgGRAPPA = zeros(Nro, Npe * nShots * R, nSlices, nDyns);
imgGRAPPA_allData = zeros(Nro, Npe * nShots * R, nSlices, nDyns);
img_indShot_ref = zeros(Nro, Npe * nShots * R, nShots, nCoils, nSlices);
motionp = {};
for sl = 1 : nSlices

    %% Get Kernels and retrospective recons for this slice

    disp(['Getting individual-shot kernels for slice ' num2str(sl) '/' num2str(nSlices)]);
    % get kernels and retrospective individual shot and all-data recons for this slice, 
    % for motion and phase correction
    Gshot = {};
    for n = 1 : nShots 

        % zero-pad and circshift the fully-sampled data for this shot
        data_Ref_zpad = zeros([Nro + grappap.Nx - 1, Npe * nShots * R + grappap.Ny * R * nShots, nCoils]);
        data_Ref_zpad((grappap.Nx-1)/2+1:end-(grappap.Nx-1)/2, 2*R*nShots+1:end-2*R*nShots, :) = ...
            sqz(data_Ref(:, :, sl, :, n));
        data_Ref_zpad = circshift(data_Ref_zpad, [0, -(n - 1)*R, 0]);

        % get the zero-filled undersampled data for this shot
        data_zfilled_zpad = zeros([Nro + grappap.Nx - 1, Npe * nShots * R + grappap.Ny * R * nShots, nCoils]);
        if nDyns >= grappap.ndyn_fit % fit to dynamic data if more than one dynamic
            data_zfilled_zpad((grappap.Nx-1)/2+1:end-(grappap.Nx-1)/2, 2*R*nShots+1:R*nShots:end-2*R*nShots, :, :) = ... % already circshifted
                mean(sqz(data_Acc(:, sl, n : nShots : end, :, 1 : grappap.ndyn_fit)), 4);     
        else
            data_zfilled_zpad((grappap.Nx-1)/2+1:end-(grappap.Nx-1)/2, 2*R*nShots+1:R*nShots:end-2*R*nShots, :, :) = ... % already circshifted
                sqz(data_Acc_ref(:, sl, n : nShots : end, :));
        end

        Gshot{n} = kernelfit_uniq_2(grappap.Nx, grappap.Ny, R*nShots, data_Ref_zpad, data_zfilled_zpad, ...
            grappap.lambda_all_convacc);

        % make a prediction (currently unused)
        data_filled_zpad = predict(data_zfilled_zpad, grappap.Nx, grappap.Ny, R*nShots, Gshot{n});
        % move the data back
        data_filled_zpad = circshift(data_filled_zpad, [0, (n - 1)*R, 0]);
        % remove zero-padding
        data_filled = data_filled_zpad((grappap.Nx-1)/2+1:end-(grappap.Nx-1)/2, 2*R*nShots+1:end-2*R*nShots, :);
        % get the image 
        img_indShot_ref(:, :, n, :, sl) = ift2(data_filled);

    end
 
    % get kernel and retrospective all-data recons for this slice
    disp(['Getting all-shot kernel for slice ' num2str(sl) '/' num2str(nSlices)]);

    % zero-pad and circshift the fully-sampled data for this shot
    data_Ref_zpad = zeros([Nro + grappap.Nx - 1, Npe * nShots * R + grappap.Ny * R * nShots, nCoils, nShots]);
    for n = 1 : nShots 
        data_Ref_zpad((grappap.Nx-1)/2+1:end-(grappap.Nx-1)/2, 2*R*nShots+1:end-2*R*nShots, :, n) = ...
            sqz(data_Ref(:, :, sl, :, n));
        data_Ref_zpad(:,:,:,n) = circshift(data_Ref_zpad(:,:,:,n), [0, -(n - 1)*R, 0]);
    end

    % get the zero-filled undersampled data
    data_zfilled_zpad = zeros([Nro + grappap.Nx - 1, Npe * nShots * R + grappap.Ny * R * nShots, nCoils, nShots]);
    if nDyns >= grappap.ndyn_fit
        for n = 1 : nShots
            data_zfilled_zpad((grappap.Nx-1)/2+1:end-(grappap.Nx-1)/2, 2*R*nShots+1:R*nShots:end-2*R*nShots, :, n) = ... % already circshifted
                mean(sqz(data_Acc(:, sl, n : nShots : end, :, 1 : grappap.ndyn_fit)), 4);            
        end
    else
        for n = 1 : nShots
            data_zfilled_zpad((grappap.Nx-1)/2+1:end-(grappap.Nx-1)/2, 2*R*nShots+1:R*nShots:end-2*R*nShots, :, n) = ... % already circshifted
                sqz(data_Acc_ref(:, sl, n : nShots : end, :));                
        end
    end

    Gall = kernelfit_uniq_2(grappap.Nx, grappap.Ny, R*nShots, data_Ref_zpad(:, :, :), data_zfilled_zpad(:, :, :), ...
        grappap.lambda_all_convacc);

    if strcmp(phsCorr, 'after')
        % loop over retrospective shots, getting GRAPPA predictions for each,
        % for phase correction below
        data_filled_shot_ref = zeros(Nro, Npe * nShots * R, nCoils, nShots, nShots);
        for n = 1 : nShots
            % get the GRAPPA prediction for just this shot from ref data,
            % for phase correction later
            data_zfilled_shot = zeros(size(data_zfilled_zpad));
            data_zfilled_shot(:, :, :, n) = data_zfilled_zpad(:, :, :, n);
            data_zfilled_shot = predict(data_zfilled_shot(:,:,:), grappap.Nx, grappap.Ny, R*nShots, Gall);
            data_zfilled_shot = reshape(data_zfilled_shot, size(data_zfilled_zpad));
            % move the data back
            for nn = 2:nShots
                data_zfilled_shot(:,:,:,nn) = circshift(data_zfilled_shot(:,:,:,nn), [0, (nn - 1)*R, 0]);
            end
            % remove zero padding
            data_shot = data_zfilled_shot((grappap.Nx-1)/2+1:end-(grappap.Nx-1)/2, 2*R*nShots+1:end-2*R*nShots, :, :);
            data_filled_shot_ref(:, :, :, :, n) = data_shot;
        end
        % calculate low-res images from the ref data
        img_lowres_shot_ref = ift2(data_filled_shot_ref .* window);
        img_fullres_shot_ref = sum(ift2(data_filled_shot_ref), 5);
    else
        img_lowres_shot_ref = [];
    end 
    
    %% Recon each dynamic
    if moco 
        % initialize rotation dictionary for this slice
        rotDict = [];
    end
    for dd = 10 %1 : nDyns 

        disp(['Starting on dyn ' num2str(dd) '/' num2str(nDyns) '; slice ' num2str(sl) '/' num2str(nSlices)]);

        img_indShot = zeros(Nro, Npe * nShots * R, nShots, nCoils);
        for n = 1 : nShots  %%% For Reconstructing Shuttered scans
            
            disp(['GRAPPA for shot ' num2str(n) '/' num2str(nShots) ' of dyn ' num2str(dd) '/' num2str(nDyns) ...
                ' of slice ' num2str(sl) '/' num2str(nSlices)]);

            data_zfilled_zpad = zeros([Nro + grappap.Nx - 1, Npe * nShots * R + grappap.Ny * R * nShots, nCoils]);
            data_zfilled_zpad((grappap.Nx-1)/2+1:end-(grappap.Nx-1)/2, 2*R*nShots+1:R*nShots:end-2*R*nShots, :) = ... % already circshifted
                sqz(data_Acc(:, sl, n : nShots : end, :, dd));  

            % make the prediction
            data_filled_zpad = predict(data_zfilled_zpad, grappap.Nx, grappap.Ny, R*nShots, Gshot{n});
            % move the data back
            data_filled_zpad = circshift(data_filled_zpad, [0, (n - 1)*R, 0]);
            % remove zero-padding
            data_filled = data_filled_zpad((grappap.Nx-1)/2+1:end-(grappap.Nx-1)/2, 2*R*nShots+1:end-2*R*nShots, :);
            % get the image 
            img_indShot(:, :, n, :) = ift2(data_filled);
 
            if moco
                % get the reference image for this shot
                if ~fullEXAcq
                    % for motion correction
                    %tmp = img_indShot_ref(:, :, :, :, sl);
                    %img_Ref_moco = ssq(tmp(:, :, :), 3);
                    img_Ref_moco = img_Ref_comb(:, :, sl);
                    %img_Ref_moco = ssq(sqz(img_indShot_ref(:, :, n, :, sl)), 3);
                    % for fMRI Runs
                    %img_Ref_moco = ssq(sqz(img_Ref(:, :, sl, :, n)), 3);
                else 
                    img_Ref_moco = img_Ref_comb(:, :, sl);
                end

                % get this shot's image
                img_moved_shot = ssq(sqz(img_indShot(:, :, n, :)), 3);

                % compare images
                if fullEXAcq
                    [~, params, rotDict] = shutter_motion(img_moved_shot, img_Ref_moco, ...
                        motionAngles, motionShifts, rotDict);
                else
                    if dd == 1 % generate and save out rotated image dictionary
                        [~, params, rotDict(:,:,:,n)] = shutter_motion(img_moved_shot, img_Ref_moco, ...
                            motionAngles, motionShifts, []);
                    else % we already have the dictionary from first dynamic
                        [~, params, ~] = shutter_motion(img_moved_shot, img_Ref_moco, ...
                            motionAngles, motionShifts, rotDict(:,:,:,n));
                    end
                end
                motionp{sl, dd, n}.Rotation = params.rotation;
                motionp{sl, dd, n}.Translation = [params.dx params.dy];

                % correct the images
                [kx, ky] = ndgrid(-1 / 2 : 1 / Nro : 1 / 2 - 1 / Nro, ...
                    -1 / 2 : 1 / (Npe * nShots * R) : 1 / 2 - 1 / (Npe * nShots * R));
                kxr = kx * cosd(-motionp{sl, dd, n}.Rotation) + ky * sind(-motionp{sl, dd, n}.Rotation);
                kyr = -kx * sind(-motionp{sl, dd, n}.Rotation) + ky * cosd(-motionp{sl, dd, n}.Rotation);
                Gtheta = Gmri([kxr(:), kyr(:)], true([Nro, Npe * nShots * R]), 'fov', [Nro, Npe * nShots * R]);
                for ii = 1 : nCoils
                    kReg = reshape(Gtheta * col(img_indShot(:, :, n, ii)), [Nro, Npe * nShots * R]);
                    kReg = kReg .* exp(1i * 2 * pi * (motionp{sl, dd, n}.Translation(1) * kxr + motionp{sl, dd, n}.Translation(2) * kyr));
                    img_indShot(:, :, n, ii) = ift2(kReg);
                end

            end

        end      

        % combine individual shot recons across coils
        imgGRAPPA(:, :, sl, dd) = ssq(img_indShot(:, :, :), 3); 

        disp(['All-Data GRAPPA Recon for dyn ' num2str(dd) '/' num2str(nDyns) '; slice ' num2str(sl) '/' num2str(nSlices)]);
        % get the zero-filled undersampled data for this shot
        data_zfilled_zpad = zeros([Nro + grappap.Nx - 1, Npe * nShots * R + grappap.Ny * R * nShots, nCoils, nShots]);
        for n = 1 : nShots
            data_zfilled_zpad((grappap.Nx-1)/2+1:end-(grappap.Nx-1)/2, 2*R*nShots+1:R*nShots:end-2*R*nShots, :, n) = ... % already circshifted
                sqz(data_Acc(:, sl, n : nShots : end, :, dd));            
        end

        if ~moco 

            % recon the data all-together
            if strcmp(phsCorr, 'none')
                data_filled_zpad = predict(data_zfilled_zpad(:,:,:), grappap.Nx, grappap.Ny, R*nShots, Gall);
                % remove zero padding
                data_filled = data_filled_zpad((grappap.Nx-1)/2+1:end-(grappap.Nx-1)/2, 2*R*nShots+1:end-2*R*nShots, :, :);
                imgGRAPPA_allData(:, :, sl, dd) = ssq(ift2(data_filled), 3);
            else
                % loop over shots, getting predictions for each shot and applying 
                % phase correction to each
                data_filled = zeros(Nro, Npe * nShots * R, nCoils, nShots);
                for n = 1 : nShots 

                    % get the GRAPPA prediction for just this shot, from dyn data
                    data_zfilled_shot = zeros(size(data_zfilled_zpad));
                    data_zfilled_shot(:, :, :, n) = data_zfilled_zpad(:, :, :, n);
                    data_filled_shot = predict(data_zfilled_shot(:,:,:), grappap.Nx, grappap.Ny, R*nShots, Gall);
                    data_filled_shot = reshape(data_filled_shot, size(data_zfilled_zpad));
                    % move the data back
                    for nn = 2 : nShots
                        data_filled_shot(:,:,:,nn) = circshift(data_filled_shot(:,:,:,nn), [0, (nn - 1)*R, 0]);
                    end
                    % remove zero padding
                    data_filled_shot = data_filled_shot((grappap.Nx-1)/2+1:end-(grappap.Nx-1)/2, 2*R*nShots+1:end-2*R*nShots, :, :);

                    % loop over target coils and shots, applying phase correction
                    for ii = 1 : nCoils 

                        for jj = 1 : nShots 

                            kReg = data_filled_shot(:, :, ii, jj);

                            if strcmp(phsCorr, 'after')

                                % calculate low-res images from the motion-corrected dyn data, for phase correction
                                img_lowres_shot = ift2(window .* kReg);
                                phsDiff = angle(img_lowres_shot_ref(:, :, ii, jj, n) .* conj(img_lowres_shot));

                                % apply phase correction
                                kReg = ft2(exp(1i * phsDiff) .* ift2(kReg));

                            end

                            data_filled(:, :, ii, jj) = data_filled(:, :, ii, jj) + kReg;

                        end

                    end

                end

                imgGRAPPA_allData(:, :, sl, dd) = ssq(ift2(data_filled(:, :, :)), 3);

            end
 
        else % moco = True

            % loop over shots, applying kernel to each shot's data individually and motion-
            % correcting the result, before summing 
            data_filled = zeros(Nro, Npe * nShots * R, nCoils, nShots);
            for n = 1 : nShots 

                % get the GRAPPA prediction for just this shot, from dyn data
                data_zfilled_shot = zeros(size(data_zfilled_zpad));
                data_zfilled_shot(:, :, :, n) = data_zfilled_zpad(:, :, :, n);
                data_filled_shot = predict(data_zfilled_shot(:,:,:), grappap.Nx, grappap.Ny, R*nShots, Gall);
                data_filled_shot = reshape(data_filled_shot, size(data_zfilled_zpad));
                % move the data back
                for nn = 2 : nShots
                    data_filled_shot(:,:,:,nn) = circshift(data_filled_shot(:,:,:,nn), [0, (nn - 1)*R, 0]);
                end
                % remove zero padding
                data_filled_shot = data_filled_shot((grappap.Nx-1)/2+1:end-(grappap.Nx-1)/2, 2*R*nShots+1:end-2*R*nShots, :, :);

                % define k-grid for rotation and translation
                [kx, ky] = ndgrid(-1 / 2 : 1 / Nro : 1 / 2 - 1 / Nro, ...
                    -1 / 2 : 1 / (Npe * R * nShots) : 1 / 2 - 1 / (Npe * R * nShots));
                kxr = kx * cosd(-motionp{sl, dd, n}.Rotation) + ky * sind(-motionp{sl, dd, n}.Rotation);
                kyr = -kx * sind(-motionp{sl, dd, n}.Rotation) + ky * cosd(-motionp{sl, dd, n}.Rotation);
                Gtheta = Gmri([kxr(:), kyr(:)], true([Nro, Npe * nShots * R]), 'fov', [Nro, Npe * nShots * R]);

                % loop over target coils and shots, applying motion correction
                kReg = zeros(Nro, Npe * nShots * R, nCoils, nShots);
                for ii = 1 : nCoils 

                    for jj = 1 : nShots 

                        kReg(:, :, ii, jj) = reshape(Gtheta * col(ift2(data_filled_shot(:, :, ii, jj))), [Nro, Npe * nShots * R]);
                        kReg(:, :, ii, jj) = kReg(:, :, ii, jj) .* exp(1i * 2 * pi * (motionp{sl, dd, n}.Translation(1) * kxr + motionp{sl, dd, n}.Translation(2) * kyr));

                        if strcmp(phsCorr, 'after')

                            % calculate low-res images from the motion-corrected dyn data, for phase correction
                            img_lowres_shot = ift2(window .* kReg(:, :, ii, jj));
                            phsDiff = angle(img_lowres_shot_ref(:, :, ii, jj, n) .* conj(img_lowres_shot));

                            % apply phase correction
                            kReg(:, :, ii, jj) = ft2(exp(1i * phsDiff) .* ift2(kReg(:, :, ii, jj)));

                        end

                    end

                end

                data_filled = data_filled + kReg;

                disp(['Motion Params for Shot ' num2str(n) '; dyn ' num2str(dd) '/' num2str(nDyns) '; slice ' num2str(sl) '/' num2str(nSlices)]);
                disp(['    Rotation: ' num2str(motionp{sl, dd, n}.Rotation) ' degrees']);
                disp(['    X-Translation: ' num2str(motionp{sl, dd, n}.Translation(1)) ' voxels']);
                disp(['    Y-Translation: ' num2str(motionp{sl, dd, n}.Translation(2)) ' voxels']);

            end
            
            imgGRAPPA_allData(:, :, sl, dd) = ssq(ift2(data_filled(:, :, :)), 3);

        end 

    end % dynamic loop
        
end % slice loop

if ~isempty(fNameSave)
    save(fNameSave, 'img_Ref', 'img_Ref_comb', 'img_indShot', 'img_indShot_ref', ...
        'imgGRAPPA', 'imgGRAPPA_allData', 'grappap', 'moco', 'motionp', 'phsCorr');
end
