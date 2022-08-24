function [fnames, params] = shutterData(dataset)

switch dataset

    case 'datalist_20210511 pos 1'
        fnames.shutterRef = {'raw_51119'}; % reference data
        fnames.shutterAcc = {'raw_51120', 'raw_51130', 'raw_51111'}; % pos 1, pos 2, pos 3
        %fnames.shutterAcc = {'raw_51120', 'raw_51122', 'raw_51130', 'raw_51111'}; % pos 1, pos 1 (3 dyn), pos 2, pos 3
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
        params.nCoils = 32;
        params.nSlices = 15;
        params.nShots = 4;
        params.Nro = 228;
        params.Npe = 19;
        params.R = 3;
        params.slIndDisp = 6; % slice index for display 
        params.fullEXAcq = true;
        params.sepRefScans = false;

    case '20211108_Anderson_344229 FullEX Series 17 tSNR'
        fnames.shutterRef = {'raw_110816'};
        fnames.shutterAcc = {'raw_110817'};
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
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = true;
        params.sepRefScans = false; 

    case '20211108_Anderson_344232 FullEX Series 9 tSNR' % Subject a
        fnames.shutterRef = {'raw_110881'};
        fnames.shutterAcc = {'raw_110891'};
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = true;
        params.sepRefScans = false; 

    case '20211108_Anderson_344232 FullEX Series 11 Visual' % Subject a
        fnames.shutterRef = {'raw_1108101'};
        fnames.shutterAcc = {'raw_1108111'};
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

    case '20211108_Anderson_344232 Shutter Series 14 tSNR' % Subject a
        fnames.shutterRef = {'raw_1108131'};
        fnames.shutterAcc = {'raw_1108141'};
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

    case '20211108_Anderson_344232 Shutter Series 16 Visual' % Subject a
        fnames.shutterRef = {'raw_1108151'};
        fnames.shutterAcc = {'raw_1108161'};
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

    case '20211116_Anderson_344321 FullEX Series 13 tSNR' % Subject b
        fnames.shutterRef = {'raw_111612'};
        fnames.shutterAcc = {'raw_111613'};
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = true;
        params.sepRefScans = false; 

    case '20211116_Anderson_344321 FullEX Series 15 Visual' % Subject b
        fnames.shutterRef = {'raw_111614'};
        fnames.shutterAcc = {'raw_111615'};
        params.nCoils = 32;
        params.nSlices = 7;
        params.nShots = 4;
        params.Nro = 276;
        params.Npe = 17; 
        params.R = 4;
        params.slIndDisp = 4;
        params.fullEXAcq = true;
        params.sepRefScans = false; 

    case '20211116_Anderson_344321 Shutter Series 18 tSNR' % Subject b
        fnames.shutterRef = {'raw_111617'};
        fnames.shutterAcc = {'raw_111618'};
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

    case '20211116_Anderson_344321 Shutter Series 20 Visual' % Subject b
        fnames.shutterRef = {'raw_111619'};
        fnames.shutterAcc = {'raw_111620'};
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

