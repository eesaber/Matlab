function [hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO(varargin)
	% Data spec:
	% [row, col]: 8735*23499, without header 
	% mlc_mag.row_mult (m/pixel) = 4.99654098 (7.2) ; MLC S (azimuth) Slant Post Spacing
	% mlc_mag.col_mult (m/pixel) = 7.2 (4.99654098) ; MLC C (range) Slant Post Spacing
	% NOTICE! The data is rotated for 90 degree ! 
	parse_ = inputParser;
	validationFcn_1_ = @(x) validateattributes(x,{'double'},{'nonnegative'}); 
	validationFcn_2_ = @(x) validateattributes(x,{'char'},{});
    validationFcn_3_ = @(x) validateattributes(x,{'logical'},{'scalar'});
    validationFcn_4_ = @(x) validateattributes(x,{'double'},{'scalar'});
    validationFcn_5_ = @(x) validateattributes(x,{'char'},{});
	addParameter(parse_,'CutBatch',[],validationFcn_1_);
	addParameter(parse_,'Test','',validationFcn_2_);
    addParameter(parse_,'ReadNewFile',0,validationFcn_3_);
    addParameter(parse_,'MissionNum',0,validationFcn_4_);
    addParameter(parse_,'Type','grd',validationFcn_5_);
	parse(parse_,varargin{:})

    % file type   
    typ = parse_.Results.Type;
    disp(['Using ' typ ' as input!'])    
	% mission
    global dir task im_size
    mission_num = parse_.Results.MissionNum;
	switch mission_num
		case 1
			disp('UAVSAR Mission: Aso volcano, 熊本、日本')
			im_size = [23499,8735];
			im_size_c = im_size.*[2 ,1];
			dir = '/media/akb/2026EF9426EF696C/raw_data/Aso_Kumamoto/';
			task = 'PiSAR2_07507_13170_009_131109_L090';
			verid = '01';
		case 2
			disp('UAVSAR Mission: HaywardFault, CA USA')
			%im_size = [21798, 13827];
			%im_size_c = [im_size(1)*2, im_size(2)];
            im_size = [3300, 15900];
			im_size_c = im_size.*[2 ,1];
			dir = '/media/akb/2026EF9426EF696C/raw_data/HaywardFault_CA/';
			task = 'Haywrd_23501_17114_003_171019_L090';
			verid = '01';
		case 3
			disp('UAVSAR Mission: SMAPVEX12, CAN')
			im_size = [17020, 11274];
			im_size_c = im_size.*[2 ,1];
			dir = '/media/akb/2026EF9426EF696C/raw_data/SMAPVEX12/';
			task = 'winnip_31605_12056_002_120705_L090';
        case 4
            disp('UAVSAR Mission: Aleutian Volcanoes')
			im_size = [53672, 17186];
			im_size_c = im_size.*[2 ,1];
			dir = '/media/akb/2026EF9426EF696C/raw_data/Aleutian/';
			task = 'aleutn_09103_09077_000_090930_L090';
			verid = '01';
        case 5
            disp('UAVSAR Mission: Manu National Park')
			im_size = [21947, 24787];
			im_size_c = im_size.*[2 ,1];
			dir = '/media/akb/2026EF9426EF696C/raw_data/ManuNationalPark/';    
			task = 'ManuNP_22024_14057_004_140501_L090';
			verid = '01';
		case 6
			disp('UAVSAR Mission: PPA')
			im_size = [4920, 34191];
			im_size_c = im_size.*[2 ,1];
			dir = '/media/akb/2026EF9426EF696C/raw_data/PPA/';
			task = 'PPAtst_18013_11046_003_110707_L090';
        case 7
            disp('UAVSAR Mission: Beaufort')
            %im_size = [17186 53672];
            im_size = [17186 53672/4];
            im_size_c = im_size.*[2 ,1];
			dir = '/media/akb/2026EF9426EF696C/raw_data/Beaufort/';
			task = 'beaufo_01105_15148_004_151006_L090';
			verid = '01';
        case 8 
            disp('UAVSAR Mission: Gulf of Mexico Scen.1')
            im_size = [3300, 31696];
            im_size_c = im_size.*[2,1];
            dir = '/media/akb/2026EF9426EF696C/raw_data/Mexico_Gulf_1/';
			task = 'GOMoil_07601_10052_102_100622_L090';
			verid = '01';
		case 9
			disp('UAVSAR Mission: Gulf of Mexico Scen.2')
            im_size = [3300, 36842];
            im_size_c = im_size.*[2,1];
            dir = '/media/akb/2026EF9426EF696C/raw_data/Mexico_Gulf_2/';
			task = 'gulfco_14010_10054_100_100623_L090';
			verid = '05';
		case 10
			disp('UAVSAR Mission: Gulf of Mexico Scen.4')
            im_size = [3300, 38688];
            im_size_c = im_size.*[2,1];
            dir = '/media/akb/2026EF9426EF696C/raw_data/Mexico_Gulf_4/';
			task = 'GOMoil_14201_10053_000_100622_L090';
			verid = '02';
		otherwise 
			error('You need to select a mission')
	end
	
	if(parse_.Results.ReadNewFile)
		disp('Parsing input file...')
		fid = fopen([dir task 'HHHH_CX_' verid '.' typ ],'r','ieee-le'); 
		hh_hh = (single(fread(fid, im_size,'real*4')));
		%hh_hh = single((fread(fid, im_size,'real*4')));
        %hh_hh = hh_hh(:,im_size(2)/2);
        
		fid = fopen([dir task 'HVHV_CX_' verid '.' typ ],'r','ieee-le'); 
		hv_hv = (single(fread(fid, im_size,'real*4')));
        %hv_hv = single((fread(fid, im_size,'real*4')));
        %hv_hv = hv_hv(:,im_size(2)/2);
        
		fid = fopen([dir task 'VVVV_CX_' verid '.' typ ],'r','ieee-le'); 
		vv_vv = (single(fread(fid, im_size,'real*4')));
        %vv_vv = single((fread(fid, im_size,'real*4')));
		%vv_vv = vv_vv(:,im_size(2)/2);
        
		fid = fopen([dir task 'HVVV_CX_' verid '.' typ ],'r','ieee-le'); 
		hv_vv = fread(fid,im_size_c,'real*4');
		hv_vv = (single(hv_vv(1:2:end, :) + 1j*hv_vv(2:2:end, :)));
        %hv_vv = single((hv_vv(1:2:end, :) + 1j*hv_vv(2:2:end, :)));
        %hv_vv = hv_vv(:,im_size(2)/2);
		
		fid = fopen([dir task 'HHVV_CX_' verid '.' typ ],'r','ieee-le'); 
		hh_vv = single(fread(fid,im_size_c,'real*4'));
		hh_vv = (hh_vv(1:2:end, :) + 1j*hh_vv(2:2:end, :));
		%hh_vv = (hh_vv(1:2:end, :) + 1j*hh_vv(2:2:end, :));
        %hh_vv = hh_vv(:,im_size(2)/2);
		
		fid = fopen([dir task 'HHHV_CX_' verid '.' typ ],'r','ieee-le'); 
		hh_hv = single(fread(fid,im_size_c,'real*4'));
		hh_hv = (hh_hv(1:2:end, :) + 1j*hh_hv(2:2:end, :));
		%hh_hv = (hh_hv(1:2:end, :) + 1j*hh_hv(2:2:end, :));
        %hh_hv = hh_hv(:,im_size(2)/2);
        
		fclose(fid) ;
		clear fid
		disp('Saving parsed file...')
		save([dir 'Covariance_s.mat'],'-v7.3', 'hh_hh', 'hv_hv', 'vv_vv', 'hh_hv', 'hh_vv', 'hv_vv');
	else
		if numel(parse_.Results.Test) ~= 0
			if exist([dir parse_.Results.Test '.mat'], 'file')
				disp(['Loading ' parse_.Results.Test '.mat  ...'])
				load([dir parse_.Results.Test '.mat']);
                im_size = size(hh_hh);
			else
				error('You have to do CutBatch First!')	
			end
		else
			disp('Loading whole image...')
			load([dir 'Covariance_s.mat']);
			%load([temp 'Covariance_d.mat']);
			%load([temp 'Covariance_ds.mat']);
		end
	end	
	fprintf('\n')
	if numel(parse_.Results.CutBatch)
		disp('This will replace test.mat\n')
    	a = parse_.Results.CutBatch;
        r_1 = a(1);
        r_2 = a(2);
        c_1 = a(3);
        c_2 = a(4);
		hh_hh = hh_hh(r_1:r_2, c_1:c_2);
		hv_hv = hv_hv(r_1:r_2, c_1:c_2);
		vv_vv = vv_vv(r_1:r_2, c_1:c_2);
		hh_hv = hh_hv(r_1:r_2, c_1:c_2);
		hh_vv = hh_vv(r_1:r_2, c_1:c_2);
		hv_vv = hv_vv(r_1:r_2, c_1:c_2);
		save([dir 'test.mat'],'-v6', 'hh_hh', 'hv_hv', 'vv_vv', 'hh_hv', 'hh_vv', 'hv_vv');
	end
	if(0)
		figure(1)
			imagesc(10*log10(hh_hh))
			set(gca,'Ydir','normal')
			title('$10 log |S_{hh}|^2$', 'Interpreter', 'latex'), colorbar, colormap jet;
			caxis([-30 20])
		figure(2)
			imagesc(10*log10(hv_hv))
			set(gca,'Ydir','normal')
			title('$10 \log |S_{hv}|^2$', 'Interpreter', 'latex'), colorbar, colormap jet;
			caxis([-30 10])
		figure(3)
			imagesc(10*log10(vv_vv))
			set(gca,'Ydir','normal')
			title('$10 \log |S_{vv}|^2$', 'Interpreter', 'latex'), colorbar, colormap gray;
			caxis([-30 20])
		figure(4)
			imagesc(10*log10(vv_vv./hh_hh))
			set(gca,'Ydir','normal')
			title('$10 \log \frac{|S_{vv}|^2}{|S_{hh}|^2}$', 'Interpreter', 'latex'), colorbar, colormap jet;
			caxis([-8 4])
			plot_para('Filename','criteria','Maximize',true)
		Pauli_g = 10*log10(hh_hh+2*hv_hv+vv_vv);
		figure(5)
			imagesc(Pauli_g)
			xlabel('Azimuth')
			set(gca,'Ydir','normal'), colorbar, colormap gray;
			caxis([-30 20])
			plot_para('Filename','My_pauli_g','Maximize',true)
		clear Pauli_g;
	end
end