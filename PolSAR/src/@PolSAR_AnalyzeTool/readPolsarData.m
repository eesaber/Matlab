function readPolsarData(obj, varargin)
    % READPOLSARDATA read polarimetric SAR data. 
    %
    % Syntax: READPOLSARDATA()
    % Other m-files required: none
    % Subfunctions: UAVSAR(obj), PALSAR(obj), ERS(obj)
    % MAT-files required: none
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------
    disp('loading data...')
    switch obj.PLATFORM
        case 'ALOS PALSAR'
            palsar(obj);
        case 'UAVSAR'
            uavsar(obj, varargin{:});
        case 'ERS-2'
            ers2(obj);
        otherwise
            error('The platform is not support!')
    end
end

%% READ ALOS-PALSAR DATA
function palsar(obj)
    % Check the image size is correct
    type([obj.INPUT_PATH '/config.txt'])
    % 1/2 |S_{hh} + S_{vv}|^2
    fid = fopen([obj.INPUT_PATH '/' 'T11.bin'],'r','ieee-le'); 
    obj.T_11 = (single(fread(fid, obj.IMAGE_SIZE,'real*4')));
    % 1/2 |S_{hh} - S_{vv}|^2
    fid = fopen([obj.INPUT_PATH '/' 'T22.bin'],'r','ieee-le'); 
    obj.T_22 = (single(fread(fid, obj.IMAGE_SIZE,'real*4')));
    % 2|S_{hv}|^2
    fid = fopen([obj.INPUT_PATH '/' 'T33.bin'],'r','ieee-le'); 
    obj.T_33 = (single(fread(fid, obj.IMAGE_SIZE,'real*4')));
    % 1/2 (S_{hh}+S_{vv}) (S_{hh}-S_{vv})^*
    fid = fopen([obj.INPUT_PATH '/' 'T12_real.bin'],'r','ieee-le'); 
    obj.T_12 = (single(fread(fid, obj.IMAGE_SIZE,'real*4')));
    fid = fopen([obj.INPUT_PATH '/' 'T12_imag.bin'],'r','ieee-le'); 
    obj.T_12 = obj.T_12 + 1j*(single(fread(fid, obj.IMAGE_SIZE,'real*4')));
    % (S_{hh}+S_{vv}) S_{hv} 
    fid = fopen([obj.INPUT_PATH '/' 'T13_real.bin'],'r','ieee-le'); 
    obj.T_13 = (single(fread(fid, obj.IMAGE_SIZE,'real*4')));
    fid = fopen([obj.INPUT_PATH '/' 'T13_imag.bin'],'r','ieee-le'); 
    obj.T_13 = obj.T_13 + 1j*(single(fread(fid, obj.IMAGE_SIZE,'real*4')));
    % (S_{hh}-S_{vv}) S_{hv} 
    fid = fopen([obj.INPUT_PATH '/' 'T23_real.bin'],'r','ieee-le'); 
    obj.T_23 = (single(fread(fid, obj.IMAGE_SIZE,'real*4')));
    fid = fopen([obj.INPUT_PATH '/' 'T23_imag.bin'],'r','ieee-le'); 
    obj.T_23 = obj.T_23 + 1j*(single(fread(fid, obj.IMAGE_SIZE,'real*4')));
    
    % Convert coherency matrix to covariance matrix
    obj.coh2cov()
end

%% READ ERS-2 DATA
function ers2(obj)
    [obj.vv_vv,~] = geotiffread(obj.INPUT_PATH);
    %obj.vv_vv = flipud(rot90(obj.vv_vv));
    obj.IMAGE_SIZE = size(obj.vv_vv);
end

%% READ UAVSAR DATA
function uavsar(obj, varargin)
    % parse input argument
    cd(obj.INPUT_PATH)
    parse_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'numeric'},{'scalar'});
    addParameter(parse_, 'offset', 0, validationFcn_1_);
    addParameter(parse_, 'size', -1, validationFcn_1_);
    parse(parse_,varargin{:})
    p_ = parse_.Results;
    % read set
    keySet = {'slc','mlc','grd'};
    valueSet = {{'slc_mag.set_rows', 'slc_mag.set_cols'}, ...
                {'mlc_pwr.set_rows', 'mlc_pwr.set_cols'},...
                {'grd_pwr.set_rows', 'grd_pwr.set_cols'}};
    dict= containers.Map(keySet,valueSet);
    %f_ext = input('Input file extension [slc, mlc, grd]:  ', 's');
    f_ext = 'mlc';
    dict = dict(f_ext);
    assert(size(split(strtrim(ls(['*.' f_ext]))),1)>0, ['Need .', f_ext, ' file.'])
    % Read UAVSAR data.
    META_NAME = split(strtrim(ls('*.ann')));
    assert(size(META_NAME,1)<2, 'There are two or more .ann file.')
    assert(size(META_NAME,1)>0, 'Need .ann file.')
    fid = fopen(cell2mat(META_NAME),'r');
    obj.IMAGE_SIZE = zeros(1,2);
    while ~feof(fid)
        tline = fgetl(fid);        
        if(isempty(tline) || tline(1)==';'), continue; end
        token = split(regexprep(tline,'\s+',' '));
        if(strcmp(token(1), dict(1)))
            obj.IMAGE_SIZE(2) = str2num(cell2mat(token(4)));
        elseif (strcmp(token(1), dict(2)))
            obj.IMAGE_SIZE(1) = str2num(cell2mat(token(4)));
        else 
            continue;
        end
        if all(obj.IMAGE_SIZE), break; end
    end
    fclose(fid);
    % check the offset and image size is correct
    assert(mod(p_.offset, obj.IMAGE_SIZE(1))<1, 'Error offset size')
    if p_.size == -1
         p_.size =  [obj.IMAGE_SIZE(1), (prod(obj.IMAGE_SIZE)-p_.offset)/obj.IMAGE_SIZE(1)];
    else
        assert(mod(p_.size, obj.IMAGE_SIZE(1))<1, 'Error image size')
        p_.size = [obj.IMAGE_SIZE(1), p_.size/obj.IMAGE_SIZE(1)];
    end
    % Read image based on image size
    fname = split(strtrim(ls(['*.' f_ext])));
    if strcmp('slc', f_ext)
        error('read slc function is under construction.')
        %obj.hh = fread(fid, [9900*2 95000],'real*4=>single');
        %hh = hh(1:2:end, :).*exp(1j*hh(2:2:end, :));
    else
        
        % |S_hh|^2
        idx = ~cellfun(@isempty, strfind(fname, 'HHHH'));
        fid = fopen(cell2mat(fname(idx)));
        fseek(fid, p_.offset, 'bof');
        obj.hh_hh = fread(fid, p_.size,'real*4=>single');
        fclose(fid);
        %}
        % |S_hv|^2
        idx = ~cellfun(@isempty, strfind(fname, 'HVHV'));
        fid = fopen(cell2mat(fname(idx)));        
        fseek(fid, p_.offset, 'bof');
        obj.hv_hv = fread(fid, p_.size,'real*4=>single');
        fclose(fid);
        % |S_vv|^2
        idx = ~cellfun(@isempty, strfind(fname, 'VVVV'));
        fid = fopen(cell2mat(fname(idx)));
        fseek(fid, p_.offset, 'bof');
        obj.vv_vv = fread(fid, p_.size,'real*4=>single');
        fclose(fid);
        % S_hh S_hv^*
        idx = ~cellfun(@isempty, strfind(fname, 'HHHV'));
        fid = fopen(cell2mat(fname(idx)));
        fseek(fid, 2*p_.offset, 'bof');
        obj.hh_hv = fread(fid, [2,1].*p_.size,'real*4=>single');
        obj.hh_hv = obj.hh_hv(1:2:end, :).*exp(1j*obj.hh_hv(2:2:end, :));
        fclose(fid);
        % S_hh S_vv^*
        idx = ~cellfun(@isempty, strfind(fname, 'HHVV'));
        fid = fopen(cell2mat(fname(idx)));
        fseek(fid, 2*p_.offset, 'bof');
        obj.hh_vv = fread(fid, [2,1].*p_.size,'real*4=>single');
        obj.hh_vv = obj.hh_vv(1:2:end, :).*exp(1j*obj.hh_vv(2:2:end, :));
        fclose(fid);
        % S_hv S_hv^*
        idx = ~cellfun(@isempty, strfind(fname, 'HVVV'));
        fid = fopen(cell2mat(fname(idx)));
        fseek(fid, 2*p_.offset, 'bof');
        obj.hv_vv = fread(fid, [2,1].*p_.size,'real*4=>single');
        obj.hv_vv = obj.hv_vv(1:2:end, :).*exp(1j*obj.hv_vv(2:2:end, :));
        fclose(fid);
        
    end
    
end