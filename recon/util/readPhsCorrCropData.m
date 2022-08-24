function data = readPhsCorrCropData(path_name, file_name, savemat, multislice, echoSeparated, sliceSeparated, slices)

% function to read in, phase correct, and then 2x crop EPI data in frequency encoding dim
% if savemat == true then it will save a .mat file containing the corrected
% data in the same path location as the raw data

coilSeparatedCorr = false; 

if multislice
    
    if savemat
        
        data = squeeze(f_read_listdata5([path_name, file_name], ...
            path_name, {'std'}));
        phc = squeeze(f_read_listdata5([path_name, file_name], ...
            path_name, {'phc'}));

        if(length(size(data)) == 3 && length(size(phc)) == 5)
            data = reshape(data, ...
                [size(data, 1), size(data, 2), 1, size(data, 3)]);
            phc = reshape(phc, ...
                [size(phc, 1), 1, size(phc, 2),...
                size(phc, 3), size(phc, 4), size(phc, 5)]);
        end
        
        if length(size(data)) == 4
            [ M, N_Full, NSl, NC] = size(data);
        else
            [ M, N_Full, NSl, nDyn, NC] = size(data);
            data = permute(data, [1 2 3 5 4]); % put dyn at end
        end
        
        if ~sliceSeparated
            
            data = epi_phc_corr_Saikat_multislice(data, ...
                permute(phc, [1 3 4 5 6 2]), echoSeparated);
            
        else
            
            for sl = 1 : NSl
                disp(num2str(sl));
                data(:, :, sl, :) = epi_phc_corr_Saikat_multislice( ...
                    squeeze(data(:, :, sl, :)), ...
                    squeeze(phc(:, sl, :, :, :, :)), echoSeparated);
            end
            
        end
        
        % crop data down
        dataAll = data(1 : 2 : end, :, :, :, :) + data(2 : 2 : end, :, :, :, :);
        
        for ii = 1 : NSl
            
            data = dataAll(:, :, ii, :, :);
            save([path_name, file_name, '_slice', num2str(ii), '.mat'], 'data');
            
        end
        
        % just return the slices the user wants
        data = dataAll(:, :, slices, :, :);
        
    else
        
        for ii = 1 : length(slices)
            
            % just load the data
            dataIn = load( ...
                [path_name, file_name, '_slice', num2str(slices(ii)), '.mat'], ...
                'data');
            if ii == 1
                data = zeros([size(dataIn.data, 1), size(dataIn.data, 2), ...
                    length(slices), size(dataIn.data, 4), size(dataIn.data, 5)]);
            end
            data(:, :, ii, :, :) = dataIn.data;
            
        end
        
    end
    
    
else
    
    %% REDFLAG May not work with dynamic data! See above for dynamic data handling

    if savemat
        
        % read raw data
        data = squeeze(f_read_listdata5([path_name, file_name], path_name, {'std'}));
        
        % read in phs corr data
        phc = squeeze(f_read_listdata5([path_name, file_name], path_name, {'phc'}));
        if length(size(phc)) == 4
            phc = reshape(phc, ...
                [size(phc,1), 1, size(phc, 2), size(phc, 3), size(phc, 4)]);
        end
        
        % apply phase correction
        if ~coilSeparatedCorr
            data = epi_phc_corr_Saikat_multislice(data, phc, echoSeparated);
        else
            for ii = 1 : size(data, 3)
                data(:,:,ii) = epi_phc_corr_Saikat_multislice( ...
                    data(:, :, ii), phc(:,ii,:,:,:), echoSeparated);
            end
        end
        
        % crop the data by 2x in freq enc dim
        data = data(1 : 2 : end, :, :) + data(2 : 2 : end, :, :);
        
        save([path_name, file_name, '.mat'], 'data');
        
    else
        
        % just load the data
        load([path_name, file_name, '.mat'], 'data');
        
    end
    
end