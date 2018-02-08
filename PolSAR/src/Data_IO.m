function [hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO(varargin)
	% Data spec:
	% [row, col]: 8735*23499, without header 
	% mlc_mag.row_mult (m/pixel) = 4.99654098 (7.2) ; MLC S (azimuth) Slant Post Spacing
	% mlc_mag.col_mult (m/pixel) = 7.2 (4.99654098) ; MLC C (range) Slant Post Spacing
	% NOTICE! The data is rotated for 90 degree ! 
	parse_ = inputParser;
	validationFcn_1_ = @(x) validateattributes(x,{'double'},{'nonnegative'}); 
	validationFcn_2_ = @(x) validateattributes(x,{'logical'},{'scalar'});
	addParameter(parse_,'CutBatch',[],validationFcn_1_);
	addParameter(parse_,'Test',0,validationFcn_2_);
	parse(parse_,varargin{:})

	% Data IO
	if isunix
		cd /home/akb/Code/Matlab
		temp = '/media/akb/2026EF9426EF696C/raw_data/PiSAR2_07507_13170_009_131109_L090_CX_01_grd/';
	end
	if(0)
		fprintf('Parsing input file...')
		fid = fopen([temp 'PiSAR2_07507_13170_009_131109_L090HHHH_CX_01.grd'],'r','ieee-le'); 
		hh_hh = single(rot90(fread(fid,[23499,8735],'real*4')));
		
		%hh_hh = sparse(rot90(fread(fid,[23499,8735],'real*4')));
		fid = fopen([temp 'PiSAR2_07507_13170_009_131109_L090HVHV_CX_01.grd'],'r','ieee-le'); 
		hv_hv = single(rot90(fread(fid,[23499,8735],'real*4')));
		%hv_hv = sparse(rot90(fread(fid,[23499,8735],'real*4')));
		fid = fopen([temp 'PiSAR2_07507_13170_009_131109_L090VVVV_CX_01.grd'],'r','ieee-le'); 
		vv_vv = single(rot90(fread(fid,[23499,8735],'real*4')));
		%vv_vv = sparse(rot90(fread(fid,[23499,8735],'real*4')));
		
		fid = fopen([temp 'PiSAR2_07507_13170_009_131109_L090HVVV_CX_01.grd'],'r','ieee-le'); 
		hv_vv = fread(fid,[23499*2,8735],'real*4');
		hv_vv = single(rot90(hv_vv(1:2:end, :) + 1j*hv_vv(2:2:end, :)));
		%hv_vv = sparse((rot90(hv_vv(1:2:end, :) + 1j*hv_vv(2:2:end, :))));
		
		fid = fopen([temp 'PiSAR2_07507_13170_009_131109_L090HHVV_CX_01.grd'],'r','ieee-le'); 
		hh_vv = fread(fid,[23499*2,8735],'real*4');
		hh_vv = single(rot90(hh_vv(1:2:end, :) + 1j*hh_vv(2:2:end, :)));
		%hh_vv = sparse(rot90(hh_vv(1:2:end, :) + 1j*hh_vv(2:2:end, :)));
		
		fid = fopen([temp 'PiSAR2_07507_13170_009_131109_L090HHHV_CX_01.grd'],'r','ieee-le'); 
		hh_hv = fread(fid,[23499*2,8735],'real*4');
		hh_hv = single(rot90(hh_hv(1:2:end, :) + 1j*hh_hv(2:2:end, :)));
		%hh_hv = sparse(rot90(hh_hv(1:2:end, :) + 1j*hh_hv(2:2:end, :)));
		fclose(fid) ;
		clear fid
		save([temp 'Covariance_ds.mat'],'-v7.3', 'hh_hh', 'hv_hv', 'vv_vv', 'hh_hv', 'hh_vv', 'hv_vv');
	else
		if parse_.Results.Test
			if exist([temp 'test.mat'], 'file')
				fprintf('Loading test.mat  ...')
				load([temp 'test.mat']);
			else
				error('You have to do CutBatch First!')	
			end
		else
			fprintf('Loading image...')
			load([temp 'Covariance.mat']);
			%load([temp 'Covariance_d.mat']);
			%load([temp 'Covariance_ds.mat']);
		end
	end	
	fprintf('\n')
	if numel(parse_.Results.CutBatch)
		disp('This will replace test.mat\n')
    	[r_1, r_2, c_1, c_2] = parse_.Results.CutBatch;
		hh_hh = hh_hh(r_1:r_2, c_1:c_2);
		hv_hv = hv_hv(r_1:r_2, c_1:c_2);
		vv_vv = vv_vv(r_1:r_2, c_1:c_2);
		hh_hv = hh_hv(r_1:r_2, c_1:c_2);
		hh_vv = hh_vv(r_1:r_2, c_1:c_2);
		hv_vv = hv_vv(r_1:r_2, c_1:c_2);
		save([temp 'test.mat'],'-v7.3', 'hh_hh', 'hv_hv', 'vv_vv', 'hh_hv', 'hh_vv', 'hv_vv');
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