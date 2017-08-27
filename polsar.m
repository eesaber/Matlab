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

if(1)	% Plot the Pauli-decomposition.
	up_ = 10; low_ = -20;
	Pauli = ones(row, col, 3);	
	t_p = 10*log10(sqrt(hh_hh + vv_vv - hh_vv - conj(hh_vv)));	% |S_vv - S_hh|^2 -> double bounce scattering 
	t_p(t_p < low_) = low_;
	t_p(t_p > up_ ) = up_;
	Pauli(:,:,1) = (t_p-low_)/(up_-low_);	
	t_p= 10*log10(4*sqrt(hv_hv));									% |S_hv|^2 -> volume scattering
	t_p(t_p < low_) = low_;
	t_p(t_p > up_ ) = up_;
	Pauli(:,:,2) = (t_p-low_)/(up_-low_);
	t_p = 10*log10(sqrt(hh_hh + vv_vv + hh_vv + conj(hh_vv)));	% |S_vv + S_hh|^2 -> single scattering
	t_p(t_p < low_) = low_;
	t_p(t_p > up_ ) = up_;
	Pauli(:,:,3) = (t_p-low_)/(up_-low_);
	figure(6)
		image(Pauli)
		set(gca,'Ydir','normal')
		%xlim([7500 8200])
		%ylim([2000 2600])
		xlabel('azimuth')
		plot_para('Filename','My_pauli_c','Maximize',true)
	clear Pauli t_p
	%{
	figure(7)
		image(imread([temp 'PiSAR2_07507_13170_009_131109_L090_CX_01_pauli.tif']))
		plot_para('Maximize',true)
	%}
end



%% Four-component decomposition
% f_ is the scattering matrix coefficient. The subscript
% f_s: surface, f_d: double-bounce, f_v: volume, f_c: helix 
f_c = 2 *abs(imag(hh_hv + hv_vv));

% Decide which volume scattering model is used.
% With three interval (-infty,-2dB) [-2dB, 2dB) [2dB, infty)

crit_fv = 10*log10(vv_vv./(hh_hh + are_z*0.001)); 
f_v = (8*hv_hv - 2*f_c) .*(crit_fv > -2) .*(crit_fv < 2);
f_v = f_v + 15/2*(hv_hv - f_c/4).*(f_v == 0);

% Decide double scattering or single scattering domainate 
B = vv_vv - 0.2*f_v - 0.25*f_c; C = hh_hh - 8/15*f_v - 0.25*f_c;
D = hh_vv - 2/15*f_v - 0.25*f_c;

% hh_hh + vv_vv - hh_vv - conj(hh_vv) > hh_hh + vv_vv + hh_vv + conj(hh_vv)
	%beta = (hh_hh + vv_vv - hh_vv - conj(hh_vv)) > (hh_hh + vv_vv + hh_vv + conj(hh_vv));
	beta = real(hh_vv < 0);
	alpha = beta + beta.*abs( (C-B)./(D-B+are_z));
	f_d = beta.*(C-B)./(abs(alpha).^2-1);
	f_s = beta.*(B - f_d);
% hh_hh + vv_vv - hh_vv - conj(hh_vv) > hh_hh + vv_vv + hh_vv + conj(hh_vv)
	%temp_a = (hh_hh + vv_vv - hh_vv - conj(hh_vv)) < (hh_hh + vv_vv + hh_vv + conj(hh_vv));
	temp_a = real(hh_vv > 0);
	temp_b = temp_a + real(temp_a.*(C-B)./(B+D +are_z));
	temp_fs = temp_a.*(C-B)./(abs(temp_b).^2 -1);

	f_d = temp_a.*(B-temp_fs) + f_d;
	f_s = temp_fs + f_s; 
	clear temp_fs;
	alpha = alpha - temp_a;
	clear temp_a
	beta = beta + temp_b;
	clear temp_b;
% The contribution from each scattering mechanism
P_s = f_s.*(1+abs(beta).^2); P_d = f_d.*(1+abs(alpha).^2);
P_t = P_s+P_d+f_s+f_v;
% The power contribution should be positive. Let the negative power be zero.
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

if(1)	% Plot the 4-component decomposition.
	FourCompo = single(zeros(row, col, 3));
	up_ = 20; low_ = -30;
	FourCompo(:,:,1) = 10*log10(P_d);
	FourCompo(:,:,2) = 10*log10(f_v);
	FourCompo(:,:,3) = 10*log10(P_s);
	clear P_d  f_v P_s f_c f_d f_s;
	FourCompo(FourCompo < low_) = low_;
	FourCompo(FourCompo > up_) = up_;
	FourCompo = (FourCompo-low_)/(up_-low_);
	figure(10)
		image(FourCompo)
		set(gca,'Ydir','normal')
		xlabel('azimuth')
		%plot_para('Filename','Four_compo','Maximize',true)
		plot_para('Maximize',true)
	clear FourCompo
end

%% Non-negative eigenvalue model-based decomposition
% Decide which volume scattering model is used.
% With three interval (-infty,-2dB) [-2dB, 2dB) [2dB, infty)hh_hh + vv_vv - hh
xi = 8; zet = 3; eta = 4; rho = 2;
Z = (hh_hh*zet + vv_vv*xi - hh_vv*conj(rho) - conj(hh_vv)*rho);
a = min( 2*hv_hv/eta, 1/2/(xi*zet - rho^2)*(Z - sqrt(Z.^2 - 4*(xi*zet - rho^2)*(hh_hh.*vv_vv - abs(hh_vv).^2))) ).*(crit_fv<-2);

xi = 3; zet = 3; eta = 2; rho = 1;
a = a + min( 2*hv_hv./ eta, 1/2/(xi*zet - rho^2)*(Z - sqrt(Z.^2 - 4*(xi*zet - rho^2)*(hh_hh.*vv_vv - abs(hh_vv).^2))) ).*(crit_fv>-2).*(crit_fv<2);
xi = 3; zet = 8; eta = 4; rho = 2;
a = a + min( 2*hv_hv./ eta, 1/2/(xi*zet - rho^2)*(Z - sqrt(Z.^2 - 4*(xi*zet - rho^2)*(hh_hh.*vv_vv - abs(hh_vv).^2))) ).*(crit_fv>2);
clear Z xi zet eta rho;
