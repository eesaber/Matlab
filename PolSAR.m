% Data size: 8735*23499, no header
% Data IO
clear,clc
read = 0;
if(read)
	temp = '/media/akb/2026EF9426EF696C/raw_data/PiSAR2_07507_13170_009_131109_L090_CX_01_grd/';
	fid = fopen([temp 'PiSAR2_07507_13170_009_131109_L090HHHH_CX_01.grd'],'r','ieee-le'); 
	hh_hh = single(rot90(fread(fid,[23499,8735],'real*4')));
	fid = fopen([temp 'PiSAR2_07507_13170_009_131109_L090HVHV_CX_01.grd'],'r','ieee-le'); 
	hv_hv = single(rot90(fread(fid,[23499,8735],'real*4')));
	fid = fopen([temp 'PiSAR2_07507_13170_009_131109_L090VVVV_CX_01.grd'],'r','ieee-le'); 
	vv_vv = single(rot90(fread(fid,[23499,8735],'real*4')));

	fid = fopen([temp 'PiSAR2_07507_13170_009_131109_L090HVVV_CX_01.grd'],'r','ieee-le'); 
	hv_vv = fread(fid,[23499*2,8735],'real*4');
	hv_vv = single(rot90(hv_vv(1:2:end, :) + 1j*hv_vv(2:2:end, :)));
	fid = fopen([temp 'PiSAR2_07507_13170_009_131109_L090HHVV_CX_01.grd'],'r','ieee-le'); 
	hh_vv = fread(fid,[23499*2,8735],'real*4');
	hh_vv = single(rot90(hh_vv(1:2:end, :) + 1j*hh_vv(2:2:end, :)));
	fid = fopen([temp 'PiSAR2_07507_13170_009_131109_L090HHHV_CX_01.grd'],'r','ieee-le'); 
	hh_hv = fread(fid,[23499*2,8735],'real*4');
	hh_hv = single(rot90(hh_hv(1:2:end, :) + 1j*hh_hv(2:2:end, :)));
	fclose(fid) ;
	clear fid
	save('Covariance.mat', 'hh_hh', 'hv_hv', 'vv_vv', 'hh_hv', 'hh_vv', 'hv_vv');
else
	load('Covariance.mat');
end	
%%
if(0)
	figure(1)
		imagesc(10*log10(hh_hh))
		set(gca,'Ydir','normal')
		title('$10 log |S_{hh}|^2$', 'Interpreter', 'latex'), colorbar, colormap gray;
		caxis([-30 20])
	figure(2)
		imagesc(10*log10(hv_hv))
		set(gca,'Ydir','normal')
		title('$10 \log |S_{hv}|^2$', 'Interpreter', 'latex'), colorbar, colormap gray;
		caxis([-30 20])
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
	Pauli_g = 10*log10(hh_hh+hv_hv+vv_vv);
	figure(5)
		imagesc(Pauli_g)
		xlabel('azimuth')
		set(gca,'Ydir','normal'), colorbar, colormap gray;
		caxis([-30 20])
		plot_para('Filename','My_pauli_g','Maximize',true)
	clear Pauli_g;
end

%% Pauli decomposition
[row, col] = size(hh_hh);
if(1)
	Pauli = ones(row, col, 3);
	Pauli(:,:,1) = sqrt(hh_hh + vv_vv - hh_vv - conj(hh_vv))/max(max(sqrt(hh_hh + vv_vv - hh_vv - conj(hh_vv))));	% |S_vv - S_hh|^2 -> double bounce scattering 
	Pauli(:,:,2) = sqrt(hv_hv)/max(max(hv_hv));									% |S_hv|^2 -> volume scattering
	Pauli(:,:,3) = sqrt(hh_hh + vv_vv + hh_vv + conj(hh_vv))/max(max(sqrt(hh_hh + vv_vv + hh_vv + conj(hh_vv))));	% |S_vv + S_hh|^2 -> single scattering
	Color = ones(row, col, 3);
	Color(:,:,1) = sqrt(hh_hh);	% |S_vv - S_hh|^2 -> double bounce scattering 
	Color(:,:,2) = sqrt(0.5*sqrt(hv_hv));									% |S_hv|^2 -> volume scattering
	Color(:,:,3) = sqrt(vv_vv);	% |S_vv + S_hh|^2 -> single scattering
	figure(6)
		image(Pauli)
		set(gca,'Ydir','normal')
		xlabel('azimuth')
		plot_para('Filename','My_pauli_c','Maximize',true)
	clear Pauli
	figure(7)
		image(Color)
		set(gca,'Ydir','normal')
		xlabel('azimuth')
		plot_para('Maximize',true)
	clear Color;
end

%% Four-component decomposition
%{
% Covariance matrix 
C_avg = zeros(3,3);
C_avg(1,1) = mean(mean(hh_hh));	C_avg(1,2) = sqrt(2)*mean(mean(hh_hv));	C_avg(1,3) = mean(mean(hh_vv));
C_avg(2,1)= conj(C_avg(1,2));	C_avg(2,2) = 2*mean(mean(hv_hv));	C_avg(2,3) = sqrt(2)*mean(mean(hv_vv));
C_avg(3,1) = conj(C_avg(1,3));	C_avg(3,2) = conj(C_avg(2,3));	C_avg(3,3) = mean(mean(vv_vv));
%}

% f_ is the scattering matrix coefficient. The subscript
% f_s: surface, f_d: double-bounce, f_v: volume, f_c: helix 
f_c = 2 * imag(hh_hv + hv_vv);

% Decide which volume scattering model is used.
% With three interval (-infty,-2dB) [-2dB, 2dB) [2dB, infty)
C_vol_ratio = 10*log10(mean(mean(vv_vv./((hh_hh==0)+hh_hh)))); 
fprintf('The ratio is %fdB \n', C_vol_ratio);
f_v = 15/2*(hv_hv - f_c/4);

% Decide double scattering or single scattering domainate 
 B = vv_vv - 0.2*f_v - 0.25*f_c; C = hh_hh - 8/15*f_v - 0.25*f_c;
 D = hh_vv - 2/15*f_v - 0.25*f_c;
fill_z = logical(zeros(row, col));
fill_z( hh_hh == 0) = 1;

% hh_hh + vv_vv - hh_vv - conj(hh_vv) > hh_hh + vv_vv + hh_vv + conj(hh_vv)
	beta = (hh_hh + vv_vv - hh_vv - conj(hh_vv)) > (hh_hh + vv_vv + hh_vv + conj(hh_vv));	
	alpha = beta + beta.*conj((C-B)./(D-B+fill_z));
	%alpha = beta + beta.*conj(( hh_hh - 8/15*f_v - 0.25*f_c - (vv_vv - 0.2*f_v - 0.25*f_c) )./(hh_vv - 2/15*f_v - 0.25*f_c - (vv_vv - 0.2*f_v - 0.25*f_c) ));
	f_d = beta.*(C-B)./(abs(alpha).^2-1);
	%f_d = beta.*(hh_vv - 2/15*f_v - 0.25*f_c - (vv_vv - 0.2*f_v - 0.25*f_c) )/(alpha-1);
	f_s = beta.*(B - f_d);
	%f_s = beta.*(vv_vv - 0.2*f_v - 0.25*f_c - f_d);
% hh_hh + vv_vv - hh_vv - conj(hh_vv) > hh_hh + vv_vv + hh_vv + conj(hh_vv)
	temp_a = (hh_hh + vv_vv - hh_vv - conj(hh_vv)) < (hh_hh + vv_vv + hh_vv + conj(hh_vv));
	temp_b = -temp_a + temp_a.*conj((C-B)./(B+D +fill_z));
	%temp_b = -temp_a + temp_a.*conj(( hh_hh - 8/15*f_v - 0.25*f_c - (vv_vv - 0.2*f_v - 0.25*f_c) )./(hh_vv - 2/15*f_v - 0.25*f_c + (vv_vv - 0.2*f_v - 0.25*f_c) ));
	temp_fs = temp_a.*(C-B)./(abs(temp_b).^2 -1);
	%temp_fs = temp_a.*((hh_vv - 2/15*f_v - 0.25*f_c + (vv_vv - 0.2*f_v - 0.25*f_c) )./(1+temp_b));
	
	f_d = temp_a.*(B-temp_fs) + f_d;
	f_s = temp_fs + f_s; 
	clear temp_fs;
	alpha = alpha + temp_a;
	clear temp_a
	beta = beta + temp_b;
	clear temp_b;
% The contribution from each scattering mechanism
P_s = f_s.*(1+abs(beta).^2); P_d = f_d.*(1+abs(alpha).^2); 
clear alpha beta B C D;
P_s(isnan(P_s)) = 0; P_d(isnan(P_d)) = 0;
f_v(isnan(f_v)) = 0; f_c(isnan(f_c)) = 0;
% P_v = f_v; P_c = f_c;
%%
%clear f_s f_d f_c hh_hh hv_hv vv_vv hh_hv hh_vv hv_vv;
FourCompo = ones(row, col, 3);
FourCompo(:,:,1) = P_d;
%clear P_d
FourCompo(:,:,2) = f_v;
%clear f_v
FourCompo(:,:,3) = P_s;
%clear P_s;
if(1)
	figure(10)
	image(FourCompo)
		set(gca,'Ydir','normal')
		xlabel('azimuth')
		plot_para('Filename','Four_compo','Maximize',true)
	clear FourCompo
end
%%
close all
qq(:,:,1) = -1;
qq(:,:,2) = 1.5;
qq(:,:,3) = 1.5;
image(qq)
img = imread([temp 'PiSAR2_07507_13170_009_131109_L090_CX_01_pauli.tif']);