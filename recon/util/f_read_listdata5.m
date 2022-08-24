
function [varargout] = f_read_listdata5(dataname_s,saveDir_s,varargin)
%[f_read_listdata5] parses LIST file and read DATA data.
%
% This file initially comes from [listparser.m].
%
% USAGE
%   [varargout] = f_read_listdata5(dataname_s,saveDir_s)
%   [std_data] = f_read_listdata5(dataname_s,saveDir_s,{'std'})
%   [phc_data,std_data] = f_read_listdata5(dataname_s,saveDir_s,{'phc','std'},loadopts)
%
% INPUT
%   dataname_s: 
%       Raw LIST/DATA filename, must be full filename e.g.,
%       'Z:\Research\PPE\DTI_NAV\Matlab\test\test_s_read_DWI2\raw_010\raw_0
%       10_cut'
%   saveDir_s:
%       Directory for saving *_attr_m.mat files.
%   varargin{1}: 
%       Cell array for data type (TYP) among ('std','phc','noi','nav','frc') to load
%   varargin{2}: 
%       loadopts, for selective loading of data type shown in varargin{1}
%   
%
% OUTPUT
%   varargout{1}: NOI data
%   varargout{2}: PHC data
%   varargout{3}: STD data
%   varargout{4}: NAV data
%   varargout{5}: FRC data
%
%
% Last modified
% 2011.04.08
%   This function is generated from [listparser.m].
% 2011.04.15.
%   Wrong 'offset' problem solved.
% 2011.04.16.
%   Add loadopts for selective loading.
% 2011.04.17.
%   This is generated from [f_read_listdata2_1.m].
% 2011.04.18.
%   This is generated from [f_read_listdata3.m].
%   This reads block of all channel data to reduce time spent on reading
%   .DATA data.
% 2011.10.31.
%   This is version 2 of [f_read_listdata4.m]. This takes care of
%   interleaved TYP such as PHC... STD... NAV... STD... NAV... etc.
% 2011.12.06.
%   Add selected_TYP of 'all' for reading all TYP data.
% 2011.12.09.
%   [f_read_listdata5.m] is generated. This is the same version as
%   [f_read_listdata4_v2.m], but generated to match the version number
%   shown in [s_read_DWI5.m]. Now [f_read_listdata5.m] is used in
%   [s_read_DWI5.m].
% 2012.02.03.
%   Save *_attr_m.mat under each raw data directory, e.g., raw_003_R3E13S6.
%   This will solve the problem reading wrong *_attr_m matrix when multiple
%   raw data were read.
%
% Ha-Kyu



%% Check input
fprintf('      function: [%s.m].\n',mfilename)
if nargin==2
    flag_loadopts = 0;
    selected_TYP = {'all'};
    fprintf('\nWhole loading for all TYP\n')
elseif nargin==3
    flag_loadopts = 0;
    selected_TYP = varargin{1}; % can be cell for multiple TYPs
    if ~iscell(selected_TYP)
        selected_TYP = {selected_TYP};
    end
    fprintf('\nWhole loading for TYP is\n')
    disp(upper(selected_TYP))
elseif nargin==4
    flag_loadopts = 1;
    selected_TYP = varargin{1};
    if ~iscell(selected_TYP)
        selected_TYP = {selected_TYP};
    end
    loadopts = varargin{2};
    fprintf('\nSelective loading for TYP is\n')
    disp(upper(selected_TYP))
else
    error('Number of input must be 2 to 4.')
end



%% Check output
if nargin==3 || nargin==4
    if nargout~=length(selected_TYP)
        error('Number of output must be the same as length(selected_TYP)')
    end
end



%% Check ' '_attr_m data if they are alread saved
 curDir_s = pwd;
% cd(saveDir_s)
% if exist('noi_attr_m.mat','file') && exist('phc_attr_m.mat','file') && ...
%         exist('std_attr_m.mat','file') && exist('nav_attr_m.mat','file') && ...
%         exist('frc_attr_m.mat','file')
%     flag_read_list = 0;
%     fprintf('Do not read LIST file\n')
% else
%     flag_read_list = 1;
%     fprintf('Read LIST file\n')
% end
% cd(curDir_s)

flag_read_list = 1; %%% SAIKAT



%% Attributes and data type
% 'size' is changed to 'siz'
% 'sign' is changed to 'sig'
% 'echo' is changed to 'ech'
attr = {'mix','dyn','card','ech','loca','chan','extr1','extr2','ky',...
    'kz','na','aver','sig','rf','grad','enc','rtop','rr','siz','offset'}; % size is changed to siz
% TYP_s = {'noi','phc','std','nav','frc'};
TYP_s = {'noi','phc','std','nav','frc','all'};



%% Read LIST file or load ' '_attr_m

if flag_read_list == 1 % read LIST file
    
    
    %% Count how many k-space lines are for each TYP
    
    % Read LIST file.
    listfilename_s = sprintf('%s.list',dataname_s);
    fid=fopen(listfilename_s,'r');
    
    % Read line.
    noi_count=0;
    phc_count=0;
    std_count=0;
    frc_count=0;
    nav_count=0;
    frx_count=0;  %%%% SAIKAT
    phx_count=0;  %%%% SAIKAT
    count=0;
    
    tic
    while ~feof(fid)
        tline = fgetl(fid);
        if strcmpi(tline(1),'#') || strcmpi(tline(1),'.')
            continue
        end
        typ = strtrim(tline(1:6));
        eval(sprintf('%s_count=%s_count+1;',lower(typ),lower(typ)))
        count=count+1;
    end
    t=toc;
    fprintf('Number of k-space lines for each TYP is counted in %.3f sec\n',t)
    
    % Number of each TYP.
    nNOI = noi_count;
    nPHC = phc_count;
    nSTD = std_count;
    nNAV = nav_count;
    nFRC = frc_count;
    
    % Check total and individual count.
    typ_count_v=[nNOI,nPHC,nSTD,nNAV,nFRC]; % [noi,phc,std,nav,frc]
    if sum(typ_count_v)~=count
        fclose(fid)
        error('Total k-space lines and sum of each TYP lines don''t match. Closing FID.')
    else
        % Rewind fid.
        frewind(fid);
    end
    fprintf('[nNOI,nPHC,nSTD,nNAV,nFRC]=[%d,%d,%d,%d,%d]\n', ...
        [nNOI,nPHC,nSTD,nNAV,nFRC])
    
    % Save typ_count_v.
    cd(saveDir_s)
    save  typ_count_v  typ_count_v
    
        
    
    %% Reserve attributes based on the number of k-space lines in each TYP
    
    % Get attributes for each TYP.
    noi_count=0;
    phc_count=0;
    std_count=0;
    nav_count=0;
    frc_count=0;
    count=0;
    noi_tmp_count=0;
    phc_tmp_count=0;
    std_tmp_count=0;
    nav_tmp_count=0;
    frc_tmp_count=0;    
    idx_first_wrong_offset=[];
    
    tic
    while ~feof(fid)
        tline=fgetl(fid);
        if strcmpi(tline(1),'#') || strcmpi(tline(1),'.')
            continue
        end
        typ = lower(strtrim(tline(1:6)));
        
        % Read each line of LIST file.
        %a=sscanf(tline(7:end),'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d');
        a=str2num(tline(7:end)); % faster than sscanf
        a = double(a);
        
        % Count.
        count=count+1; % total
        %if mod(count,1000)==0, fprintf('. '), end
        switch lower(typ) % each TYP
            case 'noi'
                noi_tmp_count = noi_tmp_count+1;
            case 'phc'
                phc_tmp_count = phc_tmp_count+1;
            case 'std'
                std_tmp_count = std_tmp_count+1;
            case 'nav'
                nav_tmp_count = nav_tmp_count+1;
            case 'frc'
                frc_tmp_count = frc_tmp_count+1;
            otherwise
                error('Unknown typ')
        end
        %tmp_count=tmp_count+1; % each TYP
        
        
        % Reserve runtime attr matrix for each typ.
        
        %**********
        % IMPORTANT
        % Never use 'single' precision zeros matrix. This will generate an
        % error for reading offset values.
        %**********
        
        %if tmp_count==1
            %fprintf('Generating attr matrix\n')
            if strcmpi(typ,'noi') && (noi_tmp_count==1)
                noi_tmp_attr_m = zeros(nNOI,length(attr)+1); % +1 for index in LIST file
                fprintf('  attr matrix is generated for %s\n',upper(typ))
            elseif strcmpi(typ,'phc') && (phc_tmp_count==1)
                phc_tmp_attr_m = zeros(nPHC,length(attr)+1);
                fprintf('  attr matrix is generated for %s\n',upper(typ))
            elseif strcmpi(typ,'std') && (std_tmp_count==1)
                std_tmp_attr_m = zeros(nSTD,length(attr)+1);
                fprintf('  attr matrix is generated for %s\n',upper(typ))
            elseif strcmpi(typ,'nav') && (nav_tmp_count==1)
                nav_tmp_attr_m = zeros(nNAV,length(attr)+1);
                fprintf('  attr matrix is generated for %s\n',upper(typ))
            elseif strcmpi(typ,'frc') && (frc_tmp_count==1)
                frc_tmp_attr_m = zeros(nFRC,length(attr)+1);
                fprintf('  attr matrix is generated for %s\n',upper(typ))
            %else
            %    error('Unknown TYP for reserving attr matrix.')
            end
        %end
        
        
        % Take each line of LIST file into attr matrix.
        switch lower(typ)
            case 'noi'
                noi_tmp_attr_m(noi_tmp_count,1) = count;
                noi_tmp_attr_m(noi_tmp_count,2:end) = a;
%                 if noi_tmp_count > 1
%                     noi_tmp_attr_m = f_validate_offset(noi_tmp_attr_m, ...
%                         noi_tmp_count,a,idx_first_wrong_offset,count);
%                 end
            case 'phc'
                phc_tmp_attr_m(phc_tmp_count,1) = count;
                phc_tmp_attr_m(phc_tmp_count,2:end) = a;
%                 if phc_tmp_count > 1
%                     phc_tmp_attr_m = f_validate_offset(phc_tmp_attr_m, ...
%                         noi_tmp_count,a,idx_first_wrong_offset,count);
%                 end
            case 'std'
                std_tmp_attr_m(std_tmp_count,1) = count;
                std_tmp_attr_m(std_tmp_count,2:end) = a;
%                 if std_tmp_count > 1
%                     std_tmp_attr_m = f_validate_offset(std_tmp_attr_m, ...
%                         noi_tmp_count,a,idx_first_wrong_offset,count);
%                 end
            case 'nav'
                nav_tmp_attr_m(nav_tmp_count,1) = count;
                nav_tmp_attr_m(nav_tmp_count,2:end) = a;
%                 if nav_tmp_count > 1
%                     nav_tmp_attr_m = f_validate_offset(nav_tmp_attr_m, ...
%                         noi_tmp_count,a,idx_first_wrong_offset,count);
%                 end
            case 'frc'
                frc_tmp_attr_m(frc_tmp_count,1) = count;
                frc_tmp_attr_m(frc_tmp_count,2:end) = a;
%                 if frc_tmp_count > 1
%                     frc_tmp_attr_m = f_validate_offset(frc_tmp_attr_m, ...
%                         noi_tmp_count,a,idx_first_wrong_offset,count);
%                 end
            otherwise
                error('Unknown typ')
        end
%         tmp_attr_m(tmp_count,1) = count;
%         tmp_attr_m(tmp_count,2:end) = a;
        

       %********** 
       % Don't use below or f_validate_offset. It will generate errors.
       %**********
       
%         % Check if offset value is correct.
%         % For reading a large LIST file, offset values are read wrong. This
%         % problem has been solved not using 'single' precision attr matrix.
%         % With 'single' precision, offset values are sometimes read wrong.
%         siz = a(19);
%         offset = a(20);
%         if tmp_count>1
%             siz_prev = tmp_attr_m(tmp_count-1,20);
%             offset_prev = tmp_attr_m(tmp_count-1,21);
%             if offset ~= (offset_prev+siz_prev)
%                 tmp_attr_m(tmp_count,21) = offset_prev+siz_prev;
%                 if isempty(idx_first_wrong_offset)
%                     idx_first_wrong_offset=count;
%                     fprintf('  WRONG first OFFSET FOUND\n')
%                     fprintf('  siz[%d],offset[%d],siz_prev[%d],offset_prev[%d],offset corrected[%d]\n',...
%                         siz,offset,siz_prev,offset_prev,(offset_prev+siz_prev))
%                 end
%             end
%         end
        
        
        % Final attr matrix.
        if strcmpi(typ,'noi') && noi_tmp_count==nNOI
            noi_attr_m = noi_tmp_attr_m;
            noi_count = noi_tmp_count; %tmp_count=0;
            %idx_first_wrong_offset = [];
            clear  noi_tmp_attr_m
        elseif strcmpi(typ,'phc') && phc_tmp_count==nPHC
            phc_attr_m = phc_tmp_attr_m;
            phc_count = phc_tmp_count; %tmp_count=0;
            %idx_first_wrong_offset =[];
            clear  tmp_attr_m
        elseif strcmpi(typ,'std') && std_tmp_count==nSTD
            std_attr_m = std_tmp_attr_m;
            std_count = std_tmp_count; %tmp_count=0;
            %idx_first_wrong_offset =[];
            clear  tmp_attr_m
        elseif strcmpi(typ,'nav') && nav_tmp_count==nNAV
            nav_attr_m = nav_tmp_attr_m;
            nav_count = nav_tmp_count; %tmp_count=0;
            %idx_first_wrong_offset =[];
            clear  tmp_attr_m
        elseif strcmpi(typ,'frc') && frc_tmp_count==nFRC
            frc_attr_m = frc_tmp_attr_m;
            frc_count = frc_tmp_count; %tmp_count=0;
            %idx_first_wrong_offset =[];
            clear  tmp_attr_m
        end
        
    end % while ~feof(fid)
    t=toc;
    
    fprintf('LIST file is parsed for each TYP in %.3f sec\n\n',t)    
    fclose(fid);
    
    
    % Save (TYP)_attr_m data.
    
    %**********
    % This is to read (TYP)_attr_m data for not parsing LIST file repeatedly for
    % selective loading. Because selective loading will work on the same raw
    % data, generate empty TYP data and save it such that selective loading can
    % recognize that it is reading the raw data more than once.
    %**********
    
    cd(saveDir_s)
    if exist('noi_attr_m','var')
        save  noi_attr_m  noi_attr_m
    else
        noi_attr_m = [];
        save  noi_attr_m  noi_attr_m
    end
    
    if exist('phc_attr_m','var')
        save  phc_attr_m  phc_attr_m
    else
        phc_attr_m = [];
        save  phc_attr_m  phc_attr_m
    end
    if exist('std_attr_m','var')
        save  std_attr_m  std_attr_m
    else
        std_attr_m = [];
        save  std_attr_m  std_attr_m
    end
    if exist('nav_attr_m','var')
        save  nav_attr_m  nav_attr_m
    else
        nav_attr_m = [];
        save  nav_attr_m  nav_attr_m
    end
    if exist('frc_attr_m','var')
        save  frc_attr_m  frc_attr_m
    else
        frc_attr_m = [];
        save  frc_attr_m  frc_attr_m
    end
    
else % Don't read LIST file. Just load (TYP)_attr_m files.
    
    % Load whole attr files.
    cd(saveDir_s)
    load  noi_attr_m
    load  phc_attr_m
    load  std_attr_m
    load  nav_attr_m
    load  frc_attr_m
    load  typ_count_v
    
end % if flag_read_list == 1



%% Generate selective attributes based on loadopts

% Do this only when there is loadopts for selective loading the raw data.
cd(curDir_s)
if flag_loadopts==1
    fprintf('Generate selective attr files\n')
    
    % Get attr name for selective loading.
    fnames = fieldnames(loadopts); % this must match attr
    
    
    % Modify field names of loadopts to match attr.
    
    %**********
    % This will include,
    % echo->ech, n.a.->na, sign->sig, size->siz
    % Then re-genrate loadopts using the modified field names.
    %**********
    
    for ind1=1:length(fnames)
        a=fnames{ind1};
        if strcmpi(a,'echo')
            fnames{ind1}='ech';
            a = getfield(loadopts,'echo');
            loadopts = rmfield(loadopts,'echo');
            loadopts.ech=a;
        elseif strcmpi(a,'n.a.')
            fnames{ind1}='na';
            a = getfield(loadopts,'n.a.');
            loadopts = rmfield(loadopts,'na');
            loadopts.ech=a;
        elseif strcmpi(a,'sign')
            fnames{ind1}='sig';
            a = getfield(loadopts,'sign');
            loadopts = rmfield(loadopts,'sig');
            loadopts.ech=a;
        elseif strcmpi(a,'size')
            fnames{ind1}='siz';
            a = getfield(loadopts,'size');
            loadopts = rmfield(loadopts,'siz');
            loadopts.ech=a;
        end
    end
    
    % Find index inside attr{}. This will be the column index in LIST file.
    % attr{fname_in_attr_index} == fnames.
    index_fname_in_attr = [];
    for ind1=1:length(fnames)
        for ind2=1:length(attr)
            a=fnames{ind1};
            b=attr{ind2};
            if strcmpi(b,a)
                index_fname_in_attr=[index_fname_in_attr;ind2];
            end
        end
    end
    
    % Check error.
    if length(fnames)~=length(index_fname_in_attr)
        error('Length of ''fnames'' and ''index_fname_in_attr'' doesn''t match')
    end
    
    
    % Re-generate (TYP)_attr_m.
    % This is to find common index for all loadopts fields.
    
    %**********
    % It applies to all TYP data.
    % (TYP)_attr_m: original
    % (TYP)_attr1_m: selected
    %**********
    
    if exist('noi_attr_m','var') && ~isempty(noi_attr_m) && ~isempty(strmatch('noi',selected_TYP))
        index_v = f_get_common_index(loadopts,noi_attr_m,fnames,index_fname_in_attr);
        noi_attr1_m = noi_attr_m(index_v,:);
    end
    if exist('phc_attr_m','var') && ~isempty(phc_attr_m) && ~isempty(strmatch('phc',selected_TYP))
        index_v = f_get_common_index(loadopts,phc_attr_m,fnames,index_fname_in_attr);
        phc_attr1_m = phc_attr_m(index_v,:);
    end
    if exist('std_attr_m','var') && ~isempty(std_attr_m) && ~isempty(strmatch('std',selected_TYP))
        index_v = f_get_common_index(loadopts,std_attr_m,fnames,index_fname_in_attr);
        std_attr1_m = std_attr_m(index_v,:);
    end
    if exist('nav_attr_m','var') && ~isempty(nav_attr_m) && ~isempty(strmatch('nav',selected_TYP))
        index_v = f_get_common_index(loadopts,nav_attr_m,fnames,index_fname_in_attr);
        nav_attr1_m = nav_attr_m(index_v,:);
    end
    if exist('frc_attr_m','var') && ~isempty(frc_attr_m) && ~isempty(strmatch('frc',selected_TYP))
        index_v = f_get_common_index(loadopts,frc_attr_m,fnames,index_fname_in_attr);
        frc_attr1_m = frc_attr_m(index_v,:);
    end
        
end



%% Set flag for block reading
% This block reading is over coil channels.
flag_blockread = 1; % default
if flag_loadopts==1
    if ~isempty(strmatch('chan',fnames)) % there is 'chan' for selective loading
        flag_blockread = 0;
    end
end
if flag_blockread==1
    fprintf('BLOCK READING\n')
end



%% Read DATA data for each TYP

% Read DATA file.
datafilename_s = sprintf('%s.data',dataname_s);
fid = fopen(datafilename_s,'r','ieee-le');
if fid==-1,
    error('cannot open %s for reading', datafilename_s);
end


% Get non-zero sized TYP.
typ_read_index = find(typ_count_v); % [noi,phc,std,nav,frc]


% Determine which TYP to read.
ind_match=[];
if strcmpi(selected_TYP,'all')
    ind_match = typ_read_index;
else
    for ind=1:length(selected_TYP)
        a = selected_TYP{ind};
        ind_match=[ind_match;strmatch(a,TYP_s)];
    end
end
typ_read_and_selected_index = intersect(typ_read_index,ind_match);


% Read each TYP of data.
for typ = typ_read_and_selected_index
    loadTYP = TYP_s{typ};
    
    % Get LIST file for each TYP.
    switch typ
        case 1
            %order_noi = {'kx','loca','chan'};            
            if flag_loadopts==1 && strcmpi(loadTYP,'noi')
                attr_m = noi_attr1_m;
            else
                attr_m = noi_attr_m;
            end
            clear  noi_attr*_m
        case 2
            %order_phc = {'kx','ky','kz','loca','dyn','card','echo','mix', ...
            %   'chan','aver','sign','grad'};            
            if flag_loadopts==1 && strcmpi(loadTYP,'phc')
                attr_m = phc_attr1_m; 
            else
                attr_m = phc_attr_m;
            end
            clear  phc_attr*_m
        case 3
            %order_std = {'kx','ky','kz','loca','dyn','card','echo','mix', ...
            %   'aver','chan','extr1','extr2'};                
            if flag_loadopts==1 && strcmpi(loadTYP,'std')
                attr_m = std_attr1_m;
            else
                attr_m = std_attr_m;
            end
            clear  std_attr*_m
        case 4
            %order_nav = {'kx','ky','kz','loca','dyn','card','echo','mix', ...
            %   'aver','chan','extr1','extr2'};
            if flag_loadopts==1 && strcmpi(loadTYP,'nav')
                attr_m = nav_attr1_m;
            else
                attr_m = nav_attr_m;
            end
            clear  nav_attr*_m
        case 5
            %order_frc = {'kx','loca','echo','mix','chan','sign'};
            if flag_loadopts==1 && strcmpi(loadTYP,'frc')
                attr_m = frc_attr1_m;
            else
                attr_m = frc_attr_m;
            end
            clear  frc_attr*_m
        otherwise
            error('Unknown typ')
    end
    
    % Report.
    fprintf('\nReading %s data, %d lines\n',upper(loadTYP),size(attr_m,1))
    
    
    % Get number of each attr and each attr's unique elements vector.
    % number of attr: nky, nkz, ... 
    % unique vector of attr: ky_unique_v, kz_unique_v, ... 
    for ind_col = [2:17,20] % [mix ~ enc,siz]
        
        % Get column index for attr{}.
        col_idx=ind_col-1;
        
        % Get the unique vector for each attr.
        v = attr_m(:,ind_col); 
        v_unique = unique(v);
        eval(sprintf('%s_unique_v = v_unique;',attr{col_idx}))        
        
        % Get the number for each attr.
        eval(sprintf('n%s = length(v_unique);',attr{col_idx}))        
        
        % Adjust the number and unique vector of each attr based on loadopts.
        if flag_loadopts==1 %&& strcmpi(TYP_s{typ},loadTYP)
            a = attr{col_idx};
            for ind1 = 1:length(fnames)
                if strcmpi(a,fnames{ind1})
                    eval(sprintf('%s_unique_v = loadopts.%s;',fnames{ind1},fnames{ind1}))
                    eval(sprintf('n%s = length(%s_unique_v);',fnames{ind1},fnames{ind1}))
                end
            end
        end        
        
        % Check error.
        if ind_col==20 && length(nsiz)~=1
            error('Length of unique size attr must be one')
        end
        if ind_col==20 % siz
            nkx = v_unique/4/2; % /4byte/2(real,imag)
        end
        clear  v_unique
    end
    
    % Report.
    %disp([nkx,nky,nkz,nloca,ndyn,ncard,nech,nmix,naver,nchan,nsig,nextr1,nextr2])
    
    
    % Reserve output data.
    switch typ
        case 1 % noi
            data = zeros(nkx,nloca,nchan,'single');
            fprintf('NOI data size is\n'), fprintf('  \n')
            disp([nkx,nloca,nchan])
        case 2 % phc
            data = zeros(nkx,nky,nkz,nloca,ndyn,ncard,nech,nmix, ...
                nchan,naver,nsig,ngrad,'single');
            fprintf('PHC data size is\n'), fprintf('  \n')
            disp([nkx,nky,nkz,nloca,ndyn,ncard,nech,nmix, ...
                nchan,naver,nsig,ngrad])
        case 3 % std
            data = zeros(nkx,nky,nkz,nloca,ndyn,ncard,nech,nmix, ...
                naver,nchan,nextr1,nextr2,'single');
            fprintf('STD data size is\n'), fprintf('  \n')
            disp([nkx,nky,nkz,nloca,ndyn,ncard,nech,nmix, ...
                naver,nchan,nextr1,nextr2])
        case 4 % nav
            data = zeros(nkx,nky,nkz,nloca,ndyn,ncard,nech,nmix, ...
                naver,nchan,nextr1,nextr2,'single');
            fprintf('NAV data size is\n'), fprintf('  \n')
            disp([nkx,nky,nkz,nloca,ndyn,ncard,nech,nmix, ...
                naver,nchan,nextr1,nextr2])
        case 5 % frc
            data = zeros(nkx,nloca,nech,nmix,nchan,nsig,'single');
            fprintf('FRC data size is\n'), fprintf('  \n')
            disp([nkx,nloca,nech,nmix,nchan,nsig])
        otherwise
            error('Unknown typ')
    end
    
    
    % Preliminary.
    clear size % remove overloaded variable
    N = size(attr_m,1);
    milepost = 0;
    addpercent = 0.05; % report at every addpercent % of reading
    
    % Set data index vector.
    if flag_blockread==0
        n_v = 1:N;
    else
        n_v = 1:nchan:N; % read channel block
    end
    
    % Read DATA data.
    tic
    for n=n_v
        
        % Get each attr value for each TYP.
        % attr value: ky, kz, loca, ...
        for ind_attr=2:21
            eval(sprintf('%s = attr_m(n,ind_attr);',attr{ind_attr-1}));
        end
        
        % Report.
        %disp([mix,dyn,card,ech,loca,chan,extr1,extr2,ky,kz,na,aver,...
        %   sig,rf,grad,enc,rtop,rr,siz,offset])
                
        %********** 
        % IMPORTANT
        % number of attr values: nky,nkz,nloca,...
        % vector of unique attr values: ky_unique_v,kz_unique_v,...
        % each of attr value: ky,kz,loca,card,...,siz,...
        %
        % For PHC,
        % sign -1 = index 1
        % sign +1 = index 2
        %**********
                
        % Read data.
        if fseek(fid,offset,'bof')==0
            
            % Measure the time for reading percentages of data.
            if flag_blockread==0
                if n==1
                    t0=clock;
                elseif n==round(N*addpercent)
                    t1=etime(clock,t0);
                    milepost = milepost+addpercent;
                    fprintf('%d%% of data are read in %.3f sec\n',milepost*100,t1)
                elseif mod(n,round(N*(milepost+addpercent)))==0
                    milepost = milepost+addpercent;
                    fprintf('%d%% ',uint8(milepost*100))
                end
            end
            
            % Read raw kx data.
            if flag_blockread==0
                tmpdata = fread(fid, siz/4,'float32');
            else
                tmpdata = fread(fid, siz/4*nchan,'float32');
            end
            tmpdata = tmpdata(1:2:end) + 1i*tmpdata(2:2:end);
            tmpdata = tmpdata * sig;
            
            % Make single complex data: size of double.
            tmpdata = single(tmpdata);
            if flag_blockread==1
                tmpdata = reshape(tmpdata,[nkx,nchan]);
            end
            
            % Generate index and read into the reserved space.
            switch typ
                case 1 % noi
                    % generate index.
                    loca_ind = loca-min(loca_unique_v)+1;
                    chan_ind = chan-min(chan_unique_v)+1;
                    if flag_blockread==1
                        chan_ind = chan_unique_v-min(chan_unique_v)+1;
                    end
                    
                    % get data.
                    data(1:length(tmpdata),loca_ind,chan_ind)=tmpdata;
                    
                case 2 % phc
                    % generate index.
                    ky_ind = ky-min(ky_unique_v)+1;
                    kz_ind = kz-min(kz_unique_v)+1;
                    loca_ind = loca-min(loca_unique_v)+1;
                    dyn_ind = dyn-min(dyn_unique_v)+1;
                    card_ind = card-min(card_unique_v)+1;
                    echo_ind = ech-min(ech_unique_v)+1;
                    mix_ind = mix-min(mix_unique_v)+1;
                    chan_ind = chan-min(chan_unique_v)+1;
                    aver_ind = aver-min(aver_unique_v)+1;
                    sign_ind = 0.5*sig+3/2; % y=1/2*x+3/2
                    grad_ind = grad-min(grad_unique_v)+1;
                    if flag_blockread==1
                        chan_ind = chan_unique_v-min(chan_unique_v)+1;
                    end
                    
                    % report.
                    %disp([ky,kz,loca,dyn,card,ech,mix,chan,aver,sig,grad])
                    %disp([ky_ind,kz_ind,loca_ind,dyn_ind,card_ind,...
                    %  echo_ind,mix_ind,chan_ind,aver_ind,sign_ind,grad_ind])
                    
                    % get data.
                    data(1:length(tmpdata),ky_ind,kz_ind,loca_ind,dyn_ind,card_ind, ...
                        echo_ind,mix_ind,chan_ind,aver_ind,sign_ind,grad_ind)=tmpdata;
                    
                case 3 % std
                    % generate index.
                    ky_ind = ky-min(ky_unique_v)+1;
                    kz_ind = kz-min(kz_unique_v)+1;
                    loca_ind = loca-min(loca_unique_v)+1;
                    dyn_ind = dyn-min(dyn_unique_v)+1;
                    card_ind = card-min(card_unique_v)+1;
                    echo_ind = ech-min(ech_unique_v)+1;
                    mix_ind = mix-min(mix_unique_v)+1;
                    aver_ind = aver-min(aver_unique_v)+1;
                    chan_ind = chan-min(chan_unique_v)+1;
                    extr1_ind = extr1-min(extr1_unique_v)+1;
                    extr2_ind = extr2-min(extr2_unique_v)+1;
                    if flag_blockread==1
                        chan_ind = chan_unique_v-min(chan_unique_v)+1;
                    end
                    
                    % report.
                    %disp([ky,kz,loca,dyn,card,ech,mix,aver,chan,extr1,extr2])
                    %disp([ky_ind,kz_ind,loca_ind,dyn_ind,card_ind, ...
                    %  echo_ind,mix_ind,aver_ind,chan_ind,extr1_ind,extr2_ind])
                    
                    % get data.
                    data(1:length(tmpdata),ky_ind,kz_ind,loca_ind,dyn_ind,card_ind, ...
                        echo_ind,mix_ind,aver_ind,chan_ind,extr1_ind,extr2_ind)=tmpdata;
                    
                case 4 % nav
                    % generate index.
                    ky_ind = ky-min(ky_unique_v)+1;
                    kz_ind = kz-min(kz_unique_v)+1;
                    loca_ind = loca-min(loca_unique_v)+1;
                    dyn_ind = dyn-min(dyn_unique_v)+1;
                    card_ind = card-min(card_unique_v)+1;
                    echo_ind = ech-min(ech_unique_v)+1;
                    mix_ind = mix-min(mix_unique_v)+1;
                    aver_ind = aver-min(aver_unique_v)+1;
                    chan_ind = chan-min(chan_unique_v)+1;
                    extr1_ind = extr1-min(extr1_unique_v)+1;
                    extr2_ind = extr2-min(extr2_unique_v)+1;
                    if flag_blockread==1
                        chan_ind = chan_unique_v-min(chan_unique_v)+1;
                    end
                    
                    % report.
                    %disp([ky_ind,kz_ind,loca_ind,dyn_ind,card_ind, ...
                    %    echo_ind,mix_ind,aver_ind,chan_ind,extr1_ind,extr2_ind])
                    
                    % get data.
                    data(1:length(tmpdata),ky_ind,kz_ind,loca_ind,dyn_ind,card_ind, ...
                        echo_ind,mix_ind,aver_ind,chan_ind,extr1_ind,extr2_ind)=tmpdata;
                    
                case 5 % frc data
                    % generate index.
                    loca_ind = loca-min(loca_unique_v)+1;
                    echo_ind = ech-min(ech_unique_v)+1;
                    mix_ind = mix-min(mix_unique_v)+1;
                    chan_ind = chan-min(chan_unique_v)+1;
                    sign_ind = 0.5*sig+3/2; % y=1/2*x+3/2
                    if flag_blockread==1
                        chan_ind = chan_unique_v-min(chan_unique_v)+1;
                    end
                    
                    % report.
                    %disp([loca_ind,echo_ind,mix_ind,chan_ind])
                    
                    % get data.
                    data(1:length(tmpdata),loca_ind,echo_ind,mix_ind,chan_ind,sign_ind)=tmpdata;
                    
                otherwise
                    error('Unknown TYP')
            end
        else
            error('Cannot FSEEK to offset=[%f] in data file',offset)
        end % if fseek(fid,offset,'bof')==0                
    end % for n=1:N
    fprintf('\n')
    
    % Clear unused data.
    clear attr_m  tmpdata
    
    % Rewind FID.
    frewind(fid);
    
    % Reserve data.
    if typ==1
        data_noi = data;
    elseif typ==2
        data_phc = data;
    elseif typ==3
        data_std = data;
    elseif typ==4
        data_nav = data;
    elseif typ==5
        data_frc = data;
    else
        error('Unknown TYP')
    end
    clear  data
    
end % for typ=typ_read_index
fprintf('\n')

% Close FID.
fclose(fid);

% Output.
if strcmpi(selected_TYP,'all')
    for ind=1:nargout%length(typ_read_and_selected_index)
        if ind <= length(typ_read_and_selected_index)
            a = TYP_s{ind};
            eval(sprintf('varargout{ind} = data_%s;',a))
        else
            varargout{ind} = [];
        end
        eval(sprintf('clear  data_%s',a))
    end
else
    for ind=1:length(selected_TYP)
        a = selected_TYP{ind};
        eval(sprintf('varargout{ind} = data_%s;',a))
        eval(sprintf('clear  data_%s',a))
    end
end

% Clear data.
clear  *attr*  *unique*  index*  fid  tmp*







%% subfunctions

% This is replaced by [f_get_common_index.m]

% function index0_v = get_common_index(loadopts,attr_m,fnames,index_fname_in_attr)
% %[get_common_index] finds common index in attr_m based on all loadopts fields.
% % USAGE:
% %   index0_v = get_common_index(loadopts,attr_m,fnames,index_fname_in_attr)
% %
% %
% % Last modified
% % 2011.04.16.
%
% for ind=1:length(index_fname_in_attr)
%     v = attr_m(:,index_fname_in_attr(ind)+1);
%     index_v = find(v==loadopts.(sprintf('%s',fnames{ind})));
%     if ind==1
%         index0_v=index_v;
%     end
%     index1_v = intersect(index0_v,index_v);
%     index0_v = index1_v;
%     clear  index1_v  index_v  v
% end
%
% return


% This is replaced by [f_validate_offset.m].
% function tmp_attr_m = validate_offset(tmp_attr_m,tmp_count,a,idx_first_wrong_offset,count)
% % Check if offset value is correct.
% % For reading a large LIST file, offset values are read wrong. This
% % problem has been solved not using 'single' precision attr matrix.
% % With 'single' precision, offset values are sometimes read wrong.
% siz = a(19);
% offset = a(20);
% if tmp_count>1
%     siz_prev = tmp_attr_m(tmp_count-1,20);
%     offset_prev = tmp_attr_m(tmp_count-1,21);
%     if offset ~= (offset_prev+siz_prev)
%         tmp_attr_m(tmp_count,21) = offset_prev+siz_prev;
%         if isempty(idx_first_wrong_offset)
%             idx_first_wrong_offset=count;
%             fprintf('  WRONG first OFFSET FOUND\n')
%             fprintf('  siz[%d],offset[%d],siz_prev[%d],offset_prev[%d],offset corrected[%d]\n',...
%                 siz,offset,siz_prev,offset_prev,(offset_prev+siz_prev))
%         end
%     end
% else
%     error('tmp_count must be > 1')
% end





