%% Data IO
[hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = data_io();

%% Pauli decomposition
[row, col] = size(hh_hh);
non_z_ind = (hh_hh ~= 0);
are_z = ~non_z_ind;

% Plot the Pauli-decomposition.
if(1)	
	up_ = 10; low_ = -20;
	Pauli = ones(row, col, 3);	
	% |S_vv - S_hh|^2 -> double bounce scattering 
	t_p = 10*log10(sqrt(hh_hh + vv_vv - hh_vv - conj(hh_vv)));	
	t_p(t_p < low_) = low_;
	t_p(t_p > up_ ) = up_;
	Pauli(:,:,1) = (t_p-low_)/(up_-low_);	
	% |S_hv|^2 -> volume scattering
	t_p= 10*log10(sqrt(4*hv_hv));
	t_p(t_p < low_) = low_;
	t_p(t_p > up_ ) = up_;
	Pauli(:,:,2) = (t_p-low_)/(up_-low_);
	% |S_vv + S_hh|^2 -> single scattering
	t_p = 10*log10(sqrt(hh_hh + vv_vv + hh_vv + conj(hh_vv)));	
	t_p(t_p < low_) = low_;
	t_p(t_p > up_ ) = up_;
	Pauli(:,:,3) = (t_p-low_)/(up_-low_);
	figure(7)
		%image(Pauli(2800:3400, 8300:8855,:))
		image(Pauli)
		set(gca,'Ydir','normal')
		%xlim([8300 8855])
		%ylim([2800 3400])
		%set(gca,'xtick',0:111:555,'ytick',0:150:600,'XTickLabel', 0:800:4000, 'YTickLabel', [0:800:2400 3000] );
		xlabel('Azimuth (m)', 'Fontsize', 40)
		ylabel('Range (m)', 'Fontsize', 40)
		plot_para('Filename','My_pauli_c','Maximize',true, 'Ratio', [4 3 1]);
	%clear Pauli t_p
	%{
	figure(7)
		image(imread([temp 'PiSAR2_07507_13170_009_131109_L090_CX_01_pauli.tif']))
		plot_para('Maximize',true)
	%}
end
																																																									
%% Four-component decomposition (option: compensate the orientation)
fprintf('Compensate the oriented angle or not? CAUTION! THIS WILL CHANGE hh_hh, vv_vv, etc PERMANENTLY. \n')
if(1)
	fprintf('Yes.  \n')
	filename = 'Four_compoR';
	%R = [1, 0, 0; 0, cos(2*x), sin(2*x); 0, -sin(2*x), cos(2*x)];
	%R = 0.5*[1+cos(2*x), sqrt(2)*sin(2*x), 1-cos(2*x); -sqrt(2)*sin(2*x), 2*cos(2*x), sqrt(2)*sin(2*x); 1-cos(2*x), -sqrt(2)*sin(2*x), 1+cos(2*x)];
	x = 0.5*atan((4*real(hv_vv - hh_hv))./(are_z + 4*hv_hv - hh_hh - vv_vv + 2*real(hh_vv)));
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
		xlabel('Azimuth')
		set(gca,'Ydir','normal'), colorbar, colormap jet;
		caxis([-30 20])
		plot_para('Maximize',true)
	figure(9)
		imagesc(10*log10(f_c))
		title('$P_h$ (dB)', 'Interpreter', 'latex')
		set(gca,'Ydir','normal'), colorbar, colormap jet;
		caxis([-30 10])
		plot_para('Maximize',true)
end

if(1)	% Plot the 4-component decomposition.
	FourCompo = single(zeros(row, col, 3));
	up_ = 20; low_ = -30;
	FourCompo(:,:,1) = 10*log10(P_d);
	FourCompo(:,:,2) = 10*log10(f_v);
	FourCompo(:,:,3) = 10*log10(P_s);
	%%clear P_d  f_v P_s f_c f_d f_s 
	f_d = P_d; f_s = P_s; 
	clear P_d  P_s   
	FourCompo(FourCompo < low_) = low_;
	FourCompo(FourCompo > up_) = up_;
	FourCompo = (FourCompo-low_)/(up_-low_);
	figure(10)
		image(FourCompo(2800:3400, 8300:8855,:))
		set(gca,'Ydir','normal')
		set(gca,'xtick',0:111:555,'ytick',0:150:600,'XTickLabel', 0:800:4000, 'YTickLabel', [0:800:2400 3000] );
		xlabel('Azimuth (m)', 'Fontsize', 40)
		ylabel('Range (m)', 'Fontsize', 40)
		plot_para('Filename',filename,'Maximize',true,'Ratio', [4 3 1]);
		%plot_para('Maximize',true)
	clear FourCompo filename
end

%% eigenvalue model-based 4-component decomposition
% Build table
if(0) 
	fprintf('Building table...')
	f_s = 0:0.1:1; f_d = 0:0.1:1; theta = 0:pi/12:pi/4;
	beta = 0.2:0.2:1; 
	chi = [];
	for n = -1 : 0.2: -0.2
		for m = 0.8: -0.2: 0
			if abs(n + 1j*m)<=1 
				if m~=0
					chi = [chi, n + 1j*m, n - 1j*m];
				else
					chi = [chi, n + 1j*m];
				end
			end
			
		end
	end

	C_v(:,:,1) = 1/8*[3, 0, 1; 0, 2, 0; 1, 0, 3];
	C_v(:,:,2) = 1/15*[8, 0, 2; 0, 4, 0; 2, 0, 3];
	C_v(:,:,3) = 1/15*[3, 0, 2; 0, 4, 0; 2, 0, 8];
	table = zeros(138600,9); 
	qq = 1;
	for k = 1 : numel(theta)
		for n = 1 : numel(f_s)
			for m = 1 : numel(f_d)
				if f_s(n) + f_d(m) > 1
					continue
				else
					f_v = 1 - f_s(n) - f_s(m);
				end
				R = 0.5*[1+cos(2*theta(k)), sqrt(2)*sin(2*theta(k)), 1-cos(2*theta(k)); -sqrt(2)*sin(2*theta(k)), 2*cos(2*theta(k)), sqrt(2)*sin(2*theta(k));......
					1-cos(2*theta(k)), -sqrt(2)*sin(2*theta(k)), 1+cos(2*theta(k))];
				for beta_it = 1: numel(beta)
					for chi_it = 1: numel(chi)
						for cv_it = 1 : 3
							C_s = 1/(abs(chi(chi_it))^2 + 1)*[abs(chi(chi_it))^2, 0, chi(chi_it); 0,0,0; conj(chi(chi_it)),0,1];
							C_d = 1/(abs(beta(beta_it)^2) + 1)*[beta(beta_it)^2, 0, beta(beta_it); 0,0,0; conj(beta(beta_it)),0,1];
							C = f_s(n)*C_s + f_d(m)*R*C_d*R.' + f_v*C_v(:,:,cv_it);
							[V, p] = eig(C);
							%alpha = acos(abs(V(1,:) + V(3,:))/sqrt(2));
							alpha = acos(abs(V(3,:)));
							p = diag(p)/sum(diag(p));
							table(qq,:) = [f_s(n), f_d(m), f_v, p.', alpha];
							qq=qq+1;
						end
					end
				end
			end
		end
		fprintf('.')
	end
	[r_,~]= size(table);
	r_ = r_/4;
	temp = zeros(r_,9,4);
	for n = 1 : 4
		temp(:,:,n) = table(1+(n-1)*r_: n*r_,:);
	end
	table = single(temp);
	save('4comp_table.mat', 'table')
	clear 
else
	fprintf('Loading table...')
	load('4comp_table.mat')
end
fprintf('\n')
%% Query
P_t = hh_hh + vv_vv + 2*hv_hv;
P_h = single(zeros(row, col));
P_s = P_h; P_d = P_h; P_v = P_h;
[t_row,~,~] = size(table);

D_1 = table(:,:,1);
D_2 = table(:,:,2);
D_3 = table(:,:,3);
D_4 = table(:,:,4);
co_1 = repmat([zeros(1,3), 0.5, 0.5, zeros(1,4)], t_row, 1); 
left = [2800 3400];
down = [8300 8855]; down_1 = down(1); down_2 = down(2);
parfor m = left(1) : left(2)
	for n = down_1: down_2
		C = [hh_hh(m,n), sqrt(2)*hh_hv(m,n), hh_vv(m,n); sqrt(2)*conj(hh_hv(m,n)), 2*hv_hv(m,n), sqrt(2)*hv_vv(m,n);......
			conj(hh_vv(m,n)), sqrt(2)*conj(hv_vv(m,n)), vv_vv(m,n)];
		[V, p] = eig(C);
		%alpha = acos(abs(V(1,:) + V(3,:))/sqrt(2));
		alpha = acos(abs(V(3,:)));
		delta = acos(abs((V(3,:)-V(1,:)+1j*sqrt(2)*V(2,:))/2)./sin(alpha));
		hel = sin(alpha.^2).*(cos(delta.^2)- sin(delta.^2));
		% f_s(n), f_d(m), f_v, p(3:-1:2).', alpha.'
		P_h(m,n) = abs(sum(hel*p));
		p = diag(p)/sum(diag(p));
		P_r = P_t(m,n) - P_h(m,n);

		eta = 0.25*atan((-2*sqrt(2)*real((V(3,:)-V(1,:)).*conj(V(2,:)))) ./ (abs(V(3,:)-V(1,:)).^2 - abs(V(2,:)).^2))+pi/4;
		eta(eta>pi/4) = eta(eta>pi/4) - pi/2;
		d_theta = abs(eta(3)-eta(2));

		if d_theta > pi/4
			d_theta = -d_theta + pi/2;
		end
		if d_theta < pi/24
			Daba = D_1;
		elseif pi/24 <= d_theta && d_theta < pi/8
			Daba = D_2;
		elseif pi/8 <= d_theta && d_theta < 5*pi/24
			Daba = D_3;
		elseif 5*pi/24 <= d_theta && d_theta <= pi/4
			Daba = D_4;
		end
		S = sum(0.5*(ones(t_row,1)*p(2:3).' - Daba(:,5:6)).^2, 2);
		[~, ind_] = min(S + sum(Daba(:,4:6).*((ones(t_row,1)*alpha - Daba(:,7:9))*2/pi).^2, 2));
		if P_r > 0
			P_s(m,n) = P_r*Daba(ind_, 1); 
			P_d(m,n) = P_r*Daba(ind_, 2); 
			P_v(m,n) = P_r*Daba(ind_, 3);
		end
	end
end
clear table D_1 D_2 D_3 D_4
%% Plot the 4-component + eigenvalue decomposition.
if(1)	
	FourCompo_eig = single(zeros(row, col, 3));
	up_ = 20; low_ = -30;
	FourCompo_eig(:,:,1) = 10*log10(P_d);
	FourCompo_eig(:,:,2) = 10*log10(P_v);
	FourCompo_eig(:,:,3) = 10*log10(P_s);
	%clear P_d  P_v P_s P_v P_h
	FourCompo_eig(FourCompo_eig < low_) = low_;
	FourCompo_eig(FourCompo_eig > up_) = up_;
	FourCompo_eig = (FourCompo_eig-low_)/(up_-low_);
	figure(18)
		image(FourCompo_eig(left(1):left(2),down(1):down(2),:))
		set(gca,'Ydir','normal')
		set(gca,'xtick',0:111:555,'ytick',0:150:600,'XTickLabel', 0:800:4000, 'YTickLabel', [0:800:2400 3000] );
		xlabel('Azimuth (m)', 'Fontsize', 40)
		ylabel('Range (m)', 'Fontsize', 40)
		plot_para('Filename','4com_eigen','Maximize',true,'Ratio', [4 3 1]);
		%plot_para('Maximize',true)
	clear FourCompo_eig
end
%%
%{
a = left(1)+317; b = down(1)+503;
labels = {'surface','double','volume','helix'};
Ceig_4 = double([P_s(a,b) P_d(a,b) P_v(a,b) P_h(a,b)]);
Com_4 = double([f_s(a,b) f_d(a,b) f_v(a,b) f_c(a,b)]);
figure(98)
ax1 = subplot(1,2,1);
p = pie(Ceig_4/sum(Ceig_4));
legend(labels(Ceig_4>0),'Location','southoutside','Orientation','horizontal')
set(gca,'FontSize',20)
title(ax1,'4-comp. with eigen','Interpreter', 'latex')
ax2 = subplot(1,2,2);
pie(Com_4/sum(Com_4))
title(ax2,'4-comp.','Interpreter', 'latex')
legend(labels(Com_4>0),'Location','southoutside','Orientation','horizontal')
set(gca,'FontSize',20)
plot_para('Maximize',true)
%}
if(1)
	figure(14)
		imagesc(10*log10(P_s(left(1):left(2),down(1):down(2),:)))
		set(gca,'Ydir','normal')
		title('$P_{s}$ (dB)', 'Interpreter', 'latex'), colorbar, colormap jet, caxis([-30 10]), plot_para('Filename','a','Maximize',true,'Ratio', [7 5 1]);
	figure(15)
		imagesc(10*log10(P_d(left(1):left(2),down(1):down(2),:)))
		set(gca,'Ydir','normal')
		title('$P_{d}$ (dB)', 'Interpreter', 'latex'), colorbar, colormap jet, caxis([-30 10]), plot_para('Filename','b','Maximize',true,'Ratio', [7 5 1]);
		
	figure(16)
		imagesc(10*log10(abs(P_v(left(1):left(2),down(1):down(2),:))))
		set(gca,'Ydir','normal')
		title('$P_{v}$ (dB)', 'Interpreter', 'latex'), colorbar, colormap jet, caxis([-30 10]), plot_para('Filename','c','Maximize',true,'Ratio', [7 5 1]);
		
	figure(17)
		imagesc(10*log10(P_h(left(1):left(2),down(1):down(2),:)))
		set(gca,'Ydir','normal')
		title('$P_{h}$ (dB)', 'Interpreter', 'latex'), colorbar, colormap jet, caxis([-30 10]), plot_para('Filename','d','Maximize',true,'Ratio', [7 5 1]);
end