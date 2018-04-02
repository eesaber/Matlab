function [hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO(varargin)
	% Data spec:
	% [row, col]: 8735*23499, without header 
	% mlc_mag.row_mult (m/pixel) = 4.99654098 (7.2) ; MLC S (azimuth) Slant Post Spacing
	% mlc_mag.col_mult (m/pixel) = 7.2 (4.99654098) ; MLC C (range) Slant Post Spacing
	% NOTICE! The data is rotated for 90 degree ! 
	parse_ = inputParser;
	validationFcn_1_ = @(x) validateattributes(x,{'double'},{'nonnegative'}); 
	validationFcn_2_ = @(x) validateattributes(x,{'logical'},{'scalar'});
    validationFcn_3_ = @(x) validateattributes(x,{'logical'},{'scalar'});
	addParameter(parse_,'CutBatch',[],validationFcn_1_);
	addParameter(parse_,'Test',0,validationFcn_2_);
    addParameter(parse_,'ReadNewFile',0,validationFcn_3_);
	parse(parse_,varargin{:})

	chk_pw()
    % file type   
    typ = '.mlc';
    disp(['Using ' typ ' as input!'])
	% mission
	mission_num = 2;
	switch mission_num
		case 1
			disp('Mission: Aso volcano, 熊本、日本')
			im_size = [23499,8735];
			im_size_c = im_size.*[2 ,1];
			dir = '/media/akb/2026EF9426EF696C/raw_data/Aso_Kumamoto/';
			task = 'PiSAR2_07507_13170_009_131109_L090';
		case 2
			disp('Mission: HaywardFault, CA USA')
			%im_size = [21798, 13827];
			%im_size_c = [im_size(1)*2, im_size(2)];
            im_size = [3300, 15900];
			im_size_c = im_size.*[2 ,1];
			dir = '/media/akb/2026EF9426EF696C/raw_data/HaywardFault_CA/';
			task = 'Haywrd_23501_17114_003_171019_L090';
		case 3
			disp('Mission: SMAPVEX12, CAN')
			im_size = [17020, 11274];
			im_size_c = im_size.*[2 ,1];
			dir = '/media/akb/2026EF9426EF696C/raw_data/SMAPVEX12/';
			task = 'winnip_31605_12056_002_120705_L090';
        case 4
            disp('Mission: Aleutian Volcanoes')
			im_size = [53672, 17186];
			im_size_c = im_size.*[2 ,1];
			dir = '/media/akb/2026EF9426EF696C/raw_data/Aleutian/';
			task = 'aleutn_09103_09077_000_090930_L090';
        case 5
            disp('Mission: Beaufort')
            
            dir = '/media/akb/2026EF9426EF696C/raw_data/ Beaufort/';
            task = 'beaufo_01105_15148_004_151006_L090';
		otherwise 
			error('You need to select a mission')
	end
	
	if(parse_.Results.ReadNewFile)
		disp('Parsing input file...')
		fid = fopen([dir task 'HHHH_CX_01' typ],'r','ieee-le'); 
		hh_hh = single(rot90(fread(fid, im_size,'real*4')));
		
		%hh_hh = sparse(rot90(fread(fid, im_size,'real*4')));
		fid = fopen([dir task 'HVHV_CX_01' typ],'r','ieee-le'); 
		hv_hv = single(rot90(fread(fid, im_size,'real*4')));
		%hv_hv = sparse(rot90(fread(fid, im_size,'real*4')));
		fid = fopen([dir task 'VVVV_CX_01' typ],'r','ieee-le'); 
		vv_vv = single(rot90(fread(fid, im_size,'real*4')));
		%vv_vv = sparse(rot90(fread(fid, im_size,'real*4')));
		
		fid = fopen([dir task 'HVVV_CX_01' typ],'r','ieee-le'); 
		hv_vv = fread(fid,im_size_c,'real*4');
		hv_vv = single(rot90(hv_vv(1:2:end, :) + 1j*hv_vv(2:2:end, :)));
		%hv_vv = sparse((rot90(hv_vv(1:2:end, :) + 1j*hv_vv(2:2:end, :))));
		
		fid = fopen([dir task 'HHVV_CX_01' typ],'r','ieee-le'); 
		hh_vv = fread(fid,im_size_c,'real*4');
		hh_vv = single(rot90(hh_vv(1:2:end, :) + 1j*hh_vv(2:2:end, :)));
		%hh_vv = sparse(rot90(hh_vv(1:2:end, :) + 1j*hh_vv(2:2:end, :)));
		
		fid = fopen([dir task 'HHHV_CX_01' typ],'r','ieee-le'); 
		hh_hv = fread(fid,im_size_c,'real*4');
		hh_hv = single(rot90(hh_hv(1:2:end, :) + 1j*hh_hv(2:2:end, :)));
		%hh_hv = sparse(rot90(hh_hv(1:2:end, :) + 1j*hh_hv(2:2:end, :)));
		fclose(fid) ;
		clear fid
		disp('Saving parsed file...')
		save([dir 'Covariance_s.mat'],'-v7.3', 'hh_hh', 'hv_hv', 'vv_vv', 'hh_hv', 'hh_vv', 'hv_vv');
	else
		if parse_.Results.Test
			if exist([dir 'test.mat'], 'file')
				disp('Loading test.mat  ...')
				load([dir 'test.mat']);
			else
				error('You have to do CutBatch First!')	
			end
		else
			disp('Loading image...')
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