function [fnames, path_name, params] = shutterData(dataset, machine)

switch machine
    case 'waglaptop'
        stem = '/Users/guille/Dropbox/ShutterEPI/';
    case 'wageno'
        stem = '/nas/home/grissowa/shutteredEPIData/';
end

switch dataset
    
    case '20171102_Datalist_Human'
        fnames.shutterRef = {'raw_1102101', 'raw_110291', 'raw_110281', ...
            'raw_110261'}; % reference data
        fnames.shutterAcc = {'raw_1102121'};
        path_name = [stem dataset '/'];
        params.nCoils = 32;
        params.nSlices = 1;
        params.nShots = 4;
        params.N = 204;
        params.R = 3;
        params.slIndDisp = 1; % slice index for display 
        params.fullEXAcq = false;
    
    case 'datalist_20200902 copy pos 1'
        fnames.shutterRef = {'raw_90216', 'raw_90215', 'raw_90214', ...
            'raw_90213'}; % reference data
        fnames.shutterAcc = {'raw_90217', 'raw_90218', 'raw_90219', ... % Kinc ON, Kinc ON, Kinc OFF, pos 1
            'raw_90229', 'raw_90230', 'raw_90231'}; % Kinc ON, Kinc ON, Kinc OFF, pos 2
        path_name = [stem '20200920_7T_Human/'];
        params.nCoils = 32;
        params.nSlices = 25;
        params.nShots = 4;
        params.N = 228;
        params.R = 3;
        params.slIndDisp = 10; % slice index for display 
        params.fullEXAcq = false;
        params.sepRefScans = true;

    case 'datalist_20200902 copy pos 2'
        fnames.shutterRef = {'raw_90228', 'raw_90227', 'raw_90226', ...
            'raw_90225'}; % reference data
        fnames.shutterAcc = {'raw_90229', 'raw_90230', 'raw_90231', ... % Kinc ON, Kinc ON, Kinc OFF, pos 2
            'raw_90217', 'raw_90218', 'raw_90219'}; % Kinc ON, Kinc ON, Kinc OFF, pos 1
        path_name = [stem '20200920_7T_Human/'];
        params.nCoils = 32;
        params.nSlices = 25;
        params.nShots = 4;
        params.N = 228;
        params.R = 3;
        params.slIndDisp = 10; % slice index for display 
        params.fullEXAcq = false;
        params.sepRefScans = true;
        
    case '20190117_phantomMotion'
        fnames.shutterRef = {'raw_117121', 'raw_117111', 'raw_117101', ...
            'raw_117091'}; % reference data, position 1
        fnames.shutterAcc = {'raw_117131', 'raw_117211'}; % position 1, position 2
        path_name = [stem dataset '/'];
        params.nCoils = 32;
        params.nSlices = 1;
        params.nShots = 4;
        params.N = 228;
        params.R = 3;
        params.slIndDisp = 1; % slice index for display 
        params.fullEXAcq = false;

    case 'datalist_20190122_pos_1'
        fnames.shutterRef = {'raw_12222', 'raw_12221', 'raw_12220', ...
            'raw_12223'}; % reference data, position 1
        fnames.shutterAcc = {'raw_12224', 'raw_12232'}; % position 1, position 2
        path_name = [stem 'datalist_20190122/'];
        params.nCoils = 32;
        params.nSlices = 1;
        params.nShots = 4;
        params.N = 228;
        params.R = 3;
        params.slIndDisp = 1; % slice index for display 
        params.fullEXAcq = false;

    case 'datalist_20190122_pos_2'
        fnames.shutterRef = {'raw_12231', 'raw_12230', 'raw_12229', ...
            'raw_12228'}; % reference data, position 2
        fnames.shutterAcc = {'raw_12232', 'raw_12224'}; % position 2, position 1
        path_name = [stem 'datalist_20190122/'];
        params.nCoils = 32;
        params.nSlices = 1;
        params.nShots = 4;
        params.N = 228;
        params.R = 3;
        params.slIndDisp = 1; % slice index for display         
        params.fullEXAcq = false;

    % Full EX data
    case 'datalist_20200902 FullEX pos 1'
        fnames.shutterRef = {'raw_90211'}; % reference data
        fnames.shutterAcc = {'raw_90212', 'raw_90223'}; % pos 1, pos 2
        path_name = [stem '20200920_7T_Human/'];
        params.nCoils = 32;
        params.nSlices = 25;
        params.nShots = 4;
        params.N = 228;
        params.R = 3;
        params.slIndDisp = 10; % slice index for display 
        params.fullEXAcq = true;

    case 'datalist_20200902 FullEX pos 2'
        fnames.shutterRef = {'raw_90222'}; % reference data
        fnames.shutterAcc = {'raw_90223', 'raw_90212'}; % pos 1, pos 2
        path_name = [stem '20200920_7T_Human/'];
        params.nCoils = 32;
        params.nSlices = 25;
        params.nShots = 4;
        params.N = 228;
        params.R = 3;
        params.slIndDisp = 10; % slice index for display 
        params.fullEXAcq = true;

    case 'datalist_20210511 pos 1'
        fnames.shutterRef = {'raw_51119'}; % reference data
        fnames.shutterAcc = {'raw_51120', 'raw_51130', 'raw_51111'}; % pos 1, pos 2, pos 3
        %fnames.shutterAcc = {'raw_51120', 'raw_51122', 'raw_51130', 'raw_51111'}; % pos 1, pos 1 (3 dyn), pos 2, pos 3
        path_name = [stem 'datalist_20210511/'];
        params.nCoils = 32;
        params.nSlices = 15;
        params.nShots = 4;
        params.Nro = 228;
        params.Npe = 19;
        params.R = 3;
        params.slIndDisp = 10; % slice index for display 
        params.fullEXAcq = false;
        params.swapRefDataShotsCoils = false;
        params.sepRefScans = false;

    case 'datalist_20210511 FullEX pos 1'
        fnames.shutterRef = {'raw_51117'}; % reference data
        fnames.shutterAcc = {'raw_51118', 'raw_51128', 'raw_51109'}; % pos 1, pos 2, pos 3
        path_name = [stem 'datalist_20210511/'];
        params.nCoils = 32;
        params.nSlices = 15;
        params.nShots = 4;
        params.Nro = 228;
        params.Npe = 19;
        params.R = 3;
        params.slIndDisp = 6; % slice index for display 
        params.fullEXAcq = true;
        params.sepRefScans = false;

    case 'datalist_20210511 pos 2'
        fnames.shutterRef = {'raw_51129'}; % reference data
        fnames.shutterAcc = {'raw_51130', 'raw_51132', 'raw_51120', 'raw_51111'}; % pos 2, pos 2 (3 dyn), pos 1, pos 3
        path_name = [stem 'datalist_20210511/'];
        params.nCoils = 32;
        params.nSlices = 15;
        params.nShots = 4;
        params.Nro = 228;
        params.Npe = 19;
        params.R = 3;
        params.slIndDisp = 10; % slice index for display 
        params.fullEXAcq = false;
        params.swapRefDataShotsCoils = true;
        params.sepRefScans = false;

    case 'datalist_20210511 FullEX pos 2'
        fnames.shutterRef = {'raw_51127'}; % reference data
        fnames.shutterAcc = {'raw_51128', 'raw_51109', 'raw_51118'}; % pos 2, pos 3, pos 1
        path_name = [stem 'datalist_20210511/'];
        params.nCoils = 32;
        params.nSlices = 15;
        params.nShots = 4;
        params.Nro = 228;
        params.Npe = 19;
        params.R = 3;
        params.slIndDisp = 6; % slice index for display 
        params.fullEXAcq = true;
        params.sepRefScans = false;

    case 'datalist_20210511 pos 3'
        fnames.shutterRef = {'raw_51110'}; % reference data
        fnames.shutterAcc = {'raw_51111', 'raw_51113', 'raw_51120', 'raw_51120'}; % pos 3, pos 3 (3 dyn), pos 1, pos 2
        path_name = [stem 'datalist_20210511/'];
        params.nCoils = 32;
        params.nSlices = 15;
        params.nShots = 4;
        params.Nro = 228;
        params.Npe = 19;
        params.R = 3;
        params.slIndDisp = 10; % slice index for display 
        params.fullEXAcq = false;
        params.swapRefDataShotsCoils = true;
        params.sepRefScans = false;

    case 'datalist_20210511 FullEX pos 3'
        fnames.shutterRef = {'raw_51108'}; % reference data
        fnames.shutterAcc = {'raw_51109', 'raw_51118', 'raw_51128'}; % pos 3, pos 1, pos 2
        path_name = [stem 'datalist_20210511/'];
        params.nCoils = 32;
        params.nSlices = 15;
        params.nShots = 4;
        params.Nro = 228;
        params.Npe = 19;
        params.R = 3;
        params.slIndDisp = 6; % slice index for display 
        params.fullEXAcq = true;
        params.sepRefScans = false;

    case 'datalist_20210624 FullEX Series 13 Breathhold'
        fnames.shutterRef = {'raw_62412'};
        fnames.shutterAcc = {'raw_62413'};
        path_name = [stem '../../sengups2/ShutterEPI/20210624_Anderson_342770/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = true;
        params.sepRefScans = false; 

    case 'datalist_20210624 FullEX Series 21 Breathhold'
        fnames.shutterRef = {'raw_62420'};
        fnames.shutterAcc = {'raw_62421'};
        path_name = [stem '../../sengups2/ShutterEPI/20210624_Anderson_342770/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = true;
        params.sepRefScans = false; 

    case 'datalist_20210624 FullEX Series 15 Visual'
        fnames.shutterRef = {'raw_62414'};
        fnames.shutterAcc = {'raw_62415'};
        path_name = [stem '../../sengups2/ShutterEPI/20210624_Anderson_342770/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = true;
        params.sepRefScans = false; 

    case 'datalist_20210624 FullEX Series 23 Visual'
        fnames.shutterRef = {'raw_62422'};
        fnames.shutterAcc = {'raw_62423'};
        path_name = [stem '../../sengups2/ShutterEPI/20210624_Anderson_342770/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = true;
        params.sepRefScans = false; 

    case 'datalist_20210624 Shutter Series 9 Breathhold'
        fnames.shutterRef = {'raw_62408'};
        fnames.shutterAcc = {'raw_62409'};
        path_name = [stem '../../sengups2/ShutterEPI/20210624_Anderson_342770/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = false;
        params.sepRefScans = false; 

    case 'datalist_20210624 Shutter Series 11 Visual'
        % REDFLAG: Can't read it in!!
        fnames.shutterRef = {'raw_62410'};
        fnames.shutterAcc = {'raw_62411'};
        path_name = [stem '../../sengups2/ShutterEPI/20210624_Anderson_342770/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = false;
        params.sepRefScans = false; 

    case 'datalist_20210624 Shutter Series 17 Breathhold'
        fnames.shutterRef = {'raw_62416'};
        fnames.shutterAcc = {'raw_62417'};
        path_name = [stem '../../sengups2/ShutterEPI/20210624_Anderson_342770/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = false;
        params.sepRefScans = false;     
        
    case 'datalist_20210624 Shutter Series 19 Visual'
        fnames.shutterRef = {'raw_62418'};
        fnames.shutterAcc = {'raw_62419'};
        path_name = [stem '../../sengups2/ShutterEPI/20210624_Anderson_342770/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = false;
        params.sepRefScans = false; 
        params.swapRefDataShotsCoils = true;

    case '20211108_Anderson_344229 FullEX Series 17 tSNR'
        fnames.shutterRef = {'raw_110816'};
        fnames.shutterAcc = {'raw_110817'};
        path_name = [stem '../../sengups2/ShutterEPI/20211108_Anderson_344229/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = true;
        params.sepRefScans = false; 

    case '20211108_Anderson_344229 FullEX Series 21 Visual'
        fnames.shutterRef = {'raw_110820'};
        fnames.shutterAcc = {'raw_110821'};
        path_name = [stem '../../sengups2/ShutterEPI/20211108_Anderson_344229/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = true;
        params.sepRefScans = false; 

    case '20211108_Anderson_344232 FullEX Series 9 tSNR' % Abitha
        fnames.shutterRef = {'raw_110881'};
        fnames.shutterAcc = {'raw_110891'};
        path_name = [stem '../../sengups2/ShutterEPI/20211108_Anderson_344232/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = true;
        params.sepRefScans = false; 

    case '20211108_Anderson_344232 FullEX Series 11 Visual' % Abitha
        fnames.shutterRef = {'raw_1108101'};
        fnames.shutterAcc = {'raw_1108111'};
        path_name = [stem '../../sengups2/ShutterEPI/20211108_Anderson_344232/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = true;
        params.sepRefScans = false; 

    case '20211108_Anderson_344229 Shutter Series 24 tSNR'
        fnames.shutterRef = {'raw_110823'};
        fnames.shutterAcc = {'raw_110824'};
        path_name = [stem '../../sengups2/ShutterEPI/20211108_Anderson_344229/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = false;
        params.sepRefScans = false; 
        params.swapRefDataShotsCoils = true;

    case '20211108_Anderson_344229 Shutter Series 26 Visual'
        fnames.shutterRef = {'raw_110825'};
        fnames.shutterAcc = {'raw_110826'};
        path_name = [stem '../../sengups2/ShutterEPI/20211108_Anderson_344229/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = false;
        params.sepRefScans = false; 
        params.swapRefDataShotsCoils = true;

    case '20211108_Anderson_344232 Shutter Series 14 tSNR' % Abitha
        fnames.shutterRef = {'raw_1108131'};
        fnames.shutterAcc = {'raw_1108141'};
        path_name = [stem '../../sengups2/ShutterEPI/20211108_Anderson_344232/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = false;
        params.sepRefScans = false; 
        params.swapRefDataShotsCoils = true;

    case '20211108_Anderson_344232 Shutter Series 16 Visual' % Abitha
        fnames.shutterRef = {'raw_1108151'};
        fnames.shutterAcc = {'raw_1108161'};
        path_name = [stem '../../sengups2/ShutterEPI/20211108_Anderson_344232/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = false;
        params.sepRefScans = false; 
        params.swapRefDataShotsCoils = true;


    % 11/16/21 Charlotte
    case '20211116_Anderson_344321 FullEX Series 13 tSNR' % Charlotte
        fnames.shutterRef = {'raw_111612'};
        fnames.shutterAcc = {'raw_111613'};
        path_name = [stem '../../sengups2/ShutterEPI/20211116_Anderson_344321/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = true;
        params.sepRefScans = false; 

    case '20211116_Anderson_344321 FullEX Series 15 Visual' % Charlotte
        fnames.shutterRef = {'raw_111614'};
        fnames.shutterAcc = {'raw_111615'};
        path_name = [stem '../../sengups2/ShutterEPI/20211116_Anderson_344321/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = true;
        params.sepRefScans = false; 

    case '20211116_Anderson_344321 Shutter Series 18 tSNR' % Charlotte
        fnames.shutterRef = {'raw_111617'};
        fnames.shutterAcc = {'raw_111618'};
        path_name = [stem '../../sengups2/ShutterEPI/20211116_Anderson_344321/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = false;
        params.sepRefScans = false; 
        params.swapRefDataShotsCoils = true;

    case '20211116_Anderson_344321 Shutter Series 20 Visual' % Charlotte
        fnames.shutterRef = {'raw_111619'};
        fnames.shutterAcc = {'raw_111620'};
        path_name = [stem '../../sengups2/ShutterEPI/20211116_Anderson_344321/'];
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = false;
        params.sepRefScans = false; 
        params.swapRefDataShotsCoils = true;

    otherwise
        error 'Unrecognized dataset!'

end

disp(['Dataset ' dataset ' selected:'])
disp(['    nShots = ' num2str(params.nShots)])
disp(['         R = ' num2str(params.R)])
disp(['       Nro = ' num2str(params.Nro)])
disp(['       Npe = ' num2str(params.Npe)])
disp(['   nSlices = ' num2str(params.nSlices)])

%% Notes:
% 20210511 Data Table
% All raw data are named : raw_511XX
% Head Position 3
% 05 NoTRF R1
% 06 NoTRF R1
% 07 NoTRF R3
% 08 TRF No Shut R1
% 09 TRF No Shut R3
% 10 Shut Interleaved Ref
% 11 Shut R3 Kinc ON
% 12 Shut R3 Kinc OFF
% 13 Shut R3 Kinc ON, 3 Dynamics
% Head Position 1
% 15 NoTRF R1
% 16 NoTRF R3
% 17 TRF No Shut R1
% 18 TRF No Shut R3
% 19 Shut Interleaved Ref
% 20 Shut R3 Kinc ON
% 21 Shut R3 Kinc OFF
% 22 Shut R3 Kinc ON, 3 Dynamics
% 23 TRF No Shut R3 (Repeat of 18)
% Head Position 2
% 25 NoTRF R1
% 26 NoTRF R3
% 27 TRF No Shut R1
% 28 TRF No Shut R3
% 29 Shut Interleaved Ref
% 30 Shut R3 Kinc ON
% 31 Shut R3 Kinc OFF
% 32 Shut R3 Kinc ON, 3 Dynamics


% 20210624 data table:
% /nas/home/sengups2/ShutterEPI/20210624_Anderson_342770. 
% 4. T1 Structural
% 5. fmap-No Shim
% 6. TRF No Shutter R1
% 7. TRF No Shutter R4
% 8. ShutterRef
% 9. Shutter R4 Breathhold
% 10. Shutter Ref
% 11. Shutter R4 Visual (export failed somehow)
% 12. TRF No Shutter R1
% 13. TRF No Shutter R4 Breathhold
% 14. TRF No Shutter R1
% 15. TRF No Shutter R4 Visual
% 16. ShutterRef
% 17. Shutter R4 Breathhold
% 18. Shutter Ref
% 19. Shutter R4 Visual
% 20. TRF No Shutter R1
% 21. TRF No Shutter R4 Breathhold
% 22. TRF No Shutter R1
% 23. TRF No Shutter R4 Visual
% scans 8-15 and 16-23 are two runs of the entire set.

% 20211108 Shuttered EPI Scan 1  Anderson_344229
% 
% 14. TRF No Shut R1
% 15. fMap
% 16. TRF No Shut R1
% 17. TRF No Shut R4 tSNR
% 18. TRF No Shut R1
% 20. TRF No Shut R1
% 21. TRF No Shut R4 fMRI

% 22. Shutter R1
% 23. Shutter Ref Interleaved
% 24. Shut R4 tSNR
% 25. Shut Ref Interleaved
% 26. Shut R4 fMRI

% 20211108 Shuttered EPI Scan 2  Anderson_344232
% 
% 4. 3D T1
% 5. TRF No Shut R1
% 6. TRF No Shut R1
% 7. fMap
% 8. TRF No Shut R1
% 9. TRF No Shut R1 tSNR
% 10. TRF No Shut R1
% 11. TRF No Shut R1 fMRI


% 12. Shutter R1
% 13. Shutter Ref Interleaved
% 14. Shut R4 tSNR
% 15. Shut Ref Interleaved
% 16. Shut R4 fMRI


%20211116 Shuttered EPI Anderson_344321
%11. TRF No Shut R1
%12. TRF No Shut R1
%13. TRF No Shut R4 tSNR
%14. TRF No Shut R1
%15. TRF No Shut R4 fMRI
%16. Shutter R1
%17. Shutter Ref Interleaved
%18. Shut R4 tSNR
%19. Shut Ref Interleaved
%20. Shut R4 fMRI