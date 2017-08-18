% Data size: 8735*23499, no header
% Data IO
clear,clc
read = 0;
temp = '/media/akb/2026EF9426EF696C/raw_data/PiSAR2_07507_13170_009_131109_L090_CX_01_grd/';
if(read)
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
	Pauli_g = 10*log10(hh_hh+2*hv_hv+vv_vv);
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
non_z_ind = (hh_hh ~= 0);
are_z = ~non_z_ind;
n_std = 1.1;
if(0)
	Pauli = ones(row, col, 3);
	temp = sqrt(hh_hh + vv_vv - hh_vv - conj(hh_vv));	% |S_vv - S_hh|^2 -> double bounce scattering 
	m = mean(mean(temp(non_z_ind))) + n_std*std2(temp(non_z_ind));
	Pauli(:,:,1) = temp/m;
	
	temp = sqrt(hv_hv);									% |S_hv|^2 -> volume scattering
	qq = reshape(temp,[1,row*col]);
	m = mean(mean(temp(non_z_ind))) + n_std*std2(temp(non_z_ind));
	Pauli(:,:,2) = temp/m;
	
	temp = sqrt(hh_hh + vv_vv + hh_vv + conj(hh_vv));	% |S_vv + S_hh|^2 -> single scattering
	qq = reshape(temp,[1,row*col]);
	m = mean(mean(temp(non_z_ind))) + n_std*std2(temp(non_z_ind));
	Pauli(:,:,3) = temp/m;
	figure(6)
		image(Pauli)
		set(gca,'Ydir','normal')
		xlabel('azimuth')
		plot_para('Filename','My_pauli_c','Maximize',true)
	clear Pauli
	figure(7)
		image(imread([temp 'PiSAR2_07507_13170_009_131109_L090_CX_01_pauli.tif']))
		plot_para('Maximize',true)
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
f_c = 2 *abs(imag(hh_hv + hv_vv));

% Decide which volume scattering model is used.
% With three interval (-infty,-2dB) [-2dB, 2dB) [2dB, infty)

crit_fv = 10*log10(vv_vv./(hh_hh + are_z*0.001)); 
f_v = (8*hv_hv - 2*f_c) .*(crit_fv > -2) .*(crit_fv < 2);
f_v = f_v + 15/2*(hv_hv - f_c/4).*(f_v == 0);
clear crit_fv;

% Decide double scattering or single scattering domainate 
B = vv_vv - 0.2*f_v - 0.25*f_c; C = hh_hh - 8/15*f_v - 0.25*f_c;
D = hh_vv - 2/15*f_v - 0.25*f_c;

figure
imagesc((hh_hh<2*hv_hv).*(vv_vv<2*hv_hv))
set(gca,'Ydir','normal'), colormap gray;

% hh_hh + vv_vv - hh_vv - conj(hh_vv) > hh_hh + vv_vv + hh_vv + conj(hh_vv)
	beta = (hh_hh + vv_vv - hh_vv - conj(hh_vv)) > (hh_hh + vv_vv + hh_vv + conj(hh_vv));	
	alpha = beta + beta.*conj((C-B)./(D-B+are_z));
	%alpha = beta + beta.*conj(( hh_hh - 8/15*f_v - 0.25*f_c - (vv_vv - 0.2*f_v - 0.25*f_c) )./(hh_vv - 2/15*f_v - 0.25*f_c - (vv_vv - 0.2*f_v - 0.25*f_c) ));
	f_d = beta.*(C-B)./(abs(alpha).^2-1);
	%f_d = beta.*(hh_vv - 2/15*f_v - 0.25*f_c - (vv_vv - 0.2*f_v - 0.25*f_c) )/(alpha-1);
	f_s = beta.*(B - f_d);
	%f_s = beta.*(vv_vv - 0.2*f_v - 0.25*f_c - f_d);
% hh_hh + vv_vv - hh_vv - conj(hh_vv) > hh_hh + vv_vv + hh_vv + conj(hh_vv)
	temp_a = (hh_hh + vv_vv - hh_vv - conj(hh_vv)) < (hh_hh + vv_vv + hh_vv + conj(hh_vv));
	temp_b = -temp_a + temp_a.*conj((C-B)./(B+D +are_z));
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
P_t = P_s+P_d+f_s+f_v;
P_s(P_s < 0) = 0; P_d( P_d < 0) = 0; f_v(f_v < 0) = 0;
clear alpha beta B C D;
if(0)
	figure(8)
		imagesc(10*log10(P_t))
		xlabel('azimuth')
		set(gca,'Ydir','normal'), colorbar, colormap gray;
		caxis([-30 20])
		plot_para('Maximize',true)
	figure(9)
		imagesc(10*log10(f_c))
		xlabel('azimuth')
		set(gca,'Ydir','normal'), colorbar, colormap gray;
		caxis([-30 20])
		plot_para('Maximize',true)
end

%%
n_std = 1;
FourCompo = single(zeros(row, col, 3));
up_ = 15; low_ = -30;
t_p = 10*log10(P_d);
t_p(t_p < low_) = low_;
t_p(t_p > up_ ) = up_;
FourCompo(:,:,1) = (t_p-low_)/(up_-low_);
%clear P_d
t_p = 10*log10(f_v);
t_p(t_p < low_) = low_;
t_p(t_p > up_ ) = up_;
FourCompo(:,:,2) = (t_p-low_)/(up_-low_);
%clear f_v
t_p = 10*log10(P_s);
t_p(t_p < low_) = low_;
t_p(t_p > up_ ) = up_;
FourCompo(:,:,3) = (t_p-low_)/(up_-low_);
close all
%clear P_s;
if(1)
	figure(10)
	image(FourCompo)
		set(gca,'Ydir','normal')
		xlabel('azimuth')
		plot_para('Filename','Four_compo','Maximize',true)
	clear FourCompo
end


	
	
