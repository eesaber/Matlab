% Data size: 8735*23499, no header
% Data IO
clear,clc
cd /home/akb/Code/Matlab
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

if(0)	% Plot the Pauli-decomposition.
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
																																																													
%% Four-component decomposition (option: compensate the orientation)
fprintf('Compensate the oriented angel or not?\n')
if(1)
	fprintf('Yes. CAUTION! THIS WILL CHANGE hh_hh, vv_vv, ETC PERMANENTLY. \n')
	filename = 'Four_compoR';
	%R = [1, 0, 0; 0, cos(2*x), sin(2*x); 0, -sin(2*x), cos(2*x)];
	%R = 0.5*[1+cos(2*x), sqrt(2)*sin(2*x), 1-cos(2*x); -sqrt(2)*sin(2*x), 2*cos(2*x), sqrt(2)*sin(2*x); 1-cos(2*x), -sqrt(2)*sin(2*x), 1+cos(2*x)];
	x = -0.5*atan((2*sqrt(2)*real(hv_vv-hh_hv))./(are_z+hh_hh+vv_vv-2*real(hh_vv)-2*hv_hv)).*non_z_ind;
	x(x>pi/2) = x(x>pi/2)-pi;
	x(x<-pi/2) = x(x<-pi/2)+pi;
	%temp = R*[hh_hh, hh_hv, hh_vv; conj(hh_vv), hv_hv, hv_vv; conj(hh_vv), conj(hv_vv), vv_vv]*R.';
	C = cell(3); R = cell(3); C_n = cell(3);
	C{1,1} = hh_hh; C{1,2} = hh_hv; C{1,3} = hh_vv;
	C{2,2} = hv_hv; C{2,3} = hv_vv; C{3,3} = vv_vv;
	clear hh_hh hh_hv hh_vv hv_hv hv_vv vv_vv 
	R{1,1} = 0.5*(1+cos(2*x)); R{1,2} = -0.5*sqrt(2)*sin(2*x); R{1,3} = 0.5*(1-cos(2*x));
	R{2,2} = cos(2*x); R{2,3} = -0.5*sqrt(2)*sin(2*x); R{3,3} = 0.5*(1+cos(2*x));
	clear x
	for n = 1 : 3
		for m = n : 3
			C_n{n,m} = single(zeros(row,col));
		end
	end
	
	for n = 1 : 3
		for m = n : 3
			for q = 1 : 3
				if n>q && q>m
					C_n{n,m} = C_n{n,m} + conj(C{q,n})*(1-mod(m+q,2)).*R{m,q};
				elseif n>q
					C_n{n,m} = C_n{n,m} + conj(C{q,n}).*R{q,m};
				elseif q>m
					C_n{n,m} = C_n{n,m} + C{n,q}*(1-mod(m+q,2)).*R{m,q};
				else
					C_n{n,m} = C_n{n,m} + C{n,q}.*R{q,m};
				end
			end
		end
	end
	for n = 1 : 3
		for m = n : 3
			C{n,m} = single(zeros(row,col));
		end
	end
	
	R{1,2} = -R{1,2}; R{2,3} = -R{2,3};
	for n = 1 : 3
		for m = n : 3
			for q = 1 : 3
				if n>q && q>m
					C{n,m} = C{n,m} + (1-2*mod((q+n),2))*R{q,n}.*conj(C_n{m,q});
				elseif q>m
					C{n,m} = C{n,m} + R{n,q}.*conj(C_n{m,q});
				elseif n>q
					C{n,m} = C{n,m} + (1-2*mod((q+n),2))*R{q,n}.*C_n{q,m};
				else
					C{n,m} = C{n,m} + R{n,q}.*C_n{q,m};
				end
			end
		end
	end
	hh_hh = real(C{1,1}); hh_hv = C{1,2}; hh_vv = C{1,3};
	hv_hv = real(C{2,2}); hv_vv = C{2,3}; vv_vv = real(C{3,3});
	clear C C_n R x
else
	filename = 'Four_compo';
	fprintf('No. \n')
end
%%
% f_ is the scattering matrix coefficient. The subscript
% f_s: surface, f_d: double-bounce, f_v: volume, f_c: helix 
f_c = 2 *abs(imag(hh_hv + hv_vv));

% Decide which volume scattering model is used.
% With three interval (-infty,-2dB) [-2dB, 2dB) [2dB, infty)

crit_fv = 10*log10(vv_vv./(hh_hh + are_z*0.001)); 
temp_domain = (crit_fv >= -2) .*(crit_fv <= 2);
f_v = (8*hv_hv - 2*f_c) .*temp_domain;
temp = f_v<0;
f_v(temp) = f_v(temp) + 2*f_c(temp);
f_temp = 15/2*(hv_hv - f_c/4).*(f_v == 0);
temp = f_temp < 0;
f_temp(temp) = f_temp(temp) + 15/8*f_c(temp);
f_v = f_temp + f_v;
clear f_temp temp

% Decide double scattering or single scattering domainate 
S = hh_hh + vv_vv + hh_vv + conj(hh_vv) - 0.5*f_v;
D = (hh_hh + vv_vv - hh_vv - conj(hh_vv)) - 2*hv_hv.*temp_domain ......
	-(7/30*f_v+0.5*f_c).*(~temp_domain); 
C = vv_vv - hh_hh + hh_vv - conj(hh_vv) - 1/6*f_v.*(crit_fv<-2) + 1/6*f_v.*(crit_fv>2);
clear temp_domain crit_fv

% surface dominates alpha = -1
temp_dom = 4*real(hh_vv) - 2* hv_hv - f_c > 0;
P_s = (S + abs(C).^2./S).*temp_dom;
P_d = (D - abs(C).^2./S).*temp_dom;
% Double-dominates beta = 1
temp_dom = 4*real(hh_vv) - 2* hv_hv - f_c < 0;
P_s = P_s + (D+abs(C).^2./D).*temp_dom;
P_d = P_d + (S-abs(C).^2./D).*temp_dom;
clear S D C
% The contribution from each scattering mechanism

% The power contribution should be positive. Let the negative power be zero.
P_t = hh_hh + vv_vv + 2*hv_hv;
temp_dom = P_s < 0;
P_s(temp_dom) = 0; 
P_d(temp_dom) = P_t(temp_dom) - f_c(temp_dom) - f_v(temp_dom);
temp_dom = P_d < 0;
P_d(temp_dom) = 0;
P_s(temp_dom) = P_t(temp_dom) - f_c(temp_dom) - f_v(temp_dom);

temp_dom = f_c+f_v > P_t;
P_s(temp_dom) = 0; P_d(temp_dom) = 0;
f_v(temp_dom) = P_t(temp_dom) - f_c(temp_dom);
f_v(f_v<0) = 0;
clear temp_dom

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
	clear P_d  f_v P_s f_c f_d f_s P_t
	FourCompo(FourCompo < low_) = low_;
	FourCompo(FourCompo > up_) = up_;
	FourCompo = (FourCompo-low_)/(up_-low_);
	figure(10)
		image(FourCompo)
		set(gca,'Ydir','normal')
		xlabel('azimuth')
		%xlim([8100 8900])
		%ylim([2800 3100])
		plot_para('Filename',filename,'Maximize',true)
		%plot_para('Maximize',true)
	%clear FourCompo
end

%% eigenvalue model-based 4-component decomposition
% Build table 

% Query

