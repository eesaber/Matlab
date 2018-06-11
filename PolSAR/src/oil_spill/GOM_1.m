% This file implements the detection method of oil-spill by using PolSAR.
% Target area is in the Gulf of Mexico, mission 8
clear 
clc
%% read data
disp('loading data...')
[hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',8);
[N_az, N_ra] = size(hh_hh);
size_N = numel(hh_hh);
span = hh_hh+vv_vv+2*hv_hv;
%% 4-component decomposition
FourComp_decomp(hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv, '4decomp')
%% Span
figure
    imagesc(10*log10(span))
    Plotsetting_GOM1([-40 0],'Colorbar_unit',[40 -70])
    annotation('rectangle',[0.125 0.4 0.04 0.525],'Color','k','Linewidth',2)
        annotation('textbox',[0.17 0.71 0.1 0.1],'String','$A_1$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.7 0.525 0.06 0.4],'Color','k','Linewidth',2)
        annotation('textbox',[0.66 0.71 0.1 0.1],'String','$A_2$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.78 0.2 0.06 0.725],'Color','k','Linewidth',2)
        annotation('textbox',[0.835 0.71 0.1 0.1],'String','$A_3$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.68 0.2 0.1 0.325],'Color','k','Linewidth',2)
        annotation('textbox',[0.64 0.4 0.1 0.1],'String','$A_4$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    plot_para('Filename','output/span', 'Maximize',true)
%% Pauli decomposition
figure
Pauli_decomp((hh_hh+vv_vv-hh_vv-conj(hh_vv)), 2*hv_hv,(hh_hh+vv_vv+hh_vv+conj(hh_vv)),'Filename','Pauli_decomp')
close all

%% Use moving Average T_mn 
mask = ones(3,3)/9;
T_11 = conv2((hh_hh+vv_vv+hh_vv+conj(hh_vv))/2, mask, 'same');
T_22 = conv2((hh_hh+vv_vv-hh_vv-conj(hh_vv))/2, mask, 'same');
T_33 = conv2(2*hv_hv, mask, 'same');
T_12 = conv2((hh_hh-vv_vv-hh_vv+conj(hh_vv))/2, mask, 'same');
T_13 = conv2(hh_hv + conj(hv_vv), mask, 'same');
T_23 = conv2(hh_hv - conj(hv_vv), mask, 'same');
%clear  hh_hh hv_hv vv_vv hh_hv hh_vv hv_vv
close all
temp_T = cat(1,cat(2, reshape(T_11,[1,1,size_N]), reshape(T_12,[1,1,size_N]), reshape(T_13,[1,1,size_N])), ......
             cat(2,reshape(conj(T_12),[1,1,size_N]), reshape(T_22,[1,1,size_N]), reshape(T_23,[1,1,size_N])),.......
             cat(2,reshape(conj(T_13),[1,1,size_N]), reshape(conj(T_23),[1,1,size_N]), reshape(T_33,[1,1,size_N])));
clear T_11 T_22 T_33 T_12 T_13 T_23 
%% Eigen-decomposition
read = 0;
if read 
    [H, alpha_bar] = Eigen_decomp('FileName', 'eigen');
else
    [~] = Eigen_decomp('T',temp_T,'Calculate',true,'FileName','eigen_avg');
end
close all
%% Obtain terrain slope in azimuth and induced angle by terrain slope. 
get_polangle(temp_T);
clear temp_T
T_11 = (hh_hh+vv_vv+hh_vv+conj(hh_vv))/2;
T_22 = (hh_hh+vv_vv-hh_vv-conj(hh_vv))/2;
T_33 = 2*hv_hv;
T_12 = (hh_hh-vv_vv-hh_vv+conj(hh_vv))/2;
T_13 = hh_hv + conj(hv_vv);
T_23 = hh_hv - conj(hv_vv);
close all
%% Indicator
%{
figure
    imagesc(H+alpha_bar/180*pi+A_1+abs(hh_vv)./sqrt(hh_hh.*vv_vv))
    Plotsetting_GOM1([0 3])
    xlabel('Azimuth (km)')
    ylabel('Range (km)') 
    plot_para('Filename','output/para_F', 'Maximize',true)
%}
%%
figure
    imagesc(2*(hv_hv-real(hh_vv))./(hh_hh + 2*hv_hv + vv_vv))
    Plotsetting_GOM1([-1 1])
    annotation('rectangle',[0.125 0.4 0.04 0.525],'Color','k','Linewidth',2)
        annotation('textbox',[0.17 0.71 0.1 0.1],'String','$A_{\mu 1}$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.7 0.525 0.06 0.4],'Color','k','Linewidth',2)
        annotation('textbox',[0.65 0.71 0.1 0.1],'String','$A_{\mu 2}$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.78 0.2 0.06 0.725],'Color','k','Linewidth',2)
        annotation('textbox',[0.8 0.9 0.1 0.1],'String','$A_{\mu 3}$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.68 0.2 0.1 0.325],'Color','k','Linewidth',2)
        annotation('textbox',[0.63 0.4 0.1 0.1],'String','$A_{\mu 4}$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    plot_para('Filename','output/para_confomty', 'Maximize',true)
%%
figure
    imagesc(-(hv_hv>abs(real(hh_vv))))
    Plotsetting_GOM1([-1 0])
    annotation('rectangle',[0.125 0.4 0.04 0.525],'Color','k','Linewidth',2)
        annotation('textbox',[0.17 0.71 0.1 0.1],'String','$A_{\mu 1}$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.7 0.525 0.06 0.4],'Color','k','Linewidth',2)
        annotation('textbox',[0.65 0.71 0.1 0.1],'String','$A_{\mu 2}$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.78 0.2 0.06 0.725],'Color','k','Linewidth',2)
        annotation('textbox',[0.8 0.9 0.1 0.1],'String','$A_{\mu 3}$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.68 0.2 0.1 0.325],'Color','k','Linewidth',2)
        annotation('textbox',[0.63 0.4 0.1 0.1],'String','$A_{\mu 4}$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    colormap gray; colorbar off
    plot_para('Filename','output/para_muller', 'Maximize',true)
%%
figure
    imagesc(abs(T_12)./sqrt(T_11.*T_22))
    Plotsetting_GOM1([0 1])
    annotation('rectangle',[0.125 0.7 0.04 0.225],'Color','k','Linewidth',2)
        annotation('textbox',[0.17 0.71 0.1 0.1],'String','$A_{\gamma 1}$','Linestyle','none','Color','w','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.7 0.525 0.06 0.4],'Color','k','Linewidth',2)
        annotation('textbox',[0.65 0.71 0.1 0.1],'String','$A_{\gamma 2}$','Linestyle','none','Color','w','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.77 0.525 0.07 0.4],'Color','k','Linewidth',2)
        annotation('textbox',[0.82 0.42 0.1 0.1],'String','$A_{\gamma 3}$','Linestyle','none','Color','w','Fontsize',40,'Interpreter', 'latex')
    plot_para('Filename','output/para_corelation12', 'Maximize',true)
%%
figure
    imagesc((T_11-T_22)./(T_11+T_22))
    Plotsetting_GOM1([0 1])
    annotation('rectangle',[0.125 0.7 0.04 0.225],'Color','k','Linewidth',2)
        annotation('textbox',[0.17 0.71 0.1 0.1],'String','$A_{k 1}$','Linestyle','none','Color','w','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.7 0.525 0.06 0.4],'Color','k','Linewidth',2)
        annotation('textbox',[0.65 0.71 0.1 0.1],'String','$A_{k 2}$','Linestyle','none','Color','w','Fontsize',40,'Interpreter', 'latex')
    plot_para('Filename','output/para_moisture', 'Maximize',true)
%%
figure
    imagesc((T_22-T_33)./(T_22+T_33))
    Plotsetting_GOM1([0 1])
    xlabel('Azimuth (km)')
    ylabel('Range (km)') 
    plot_para('Filename','output/para_roughness', 'Maximize',true)
%%
figure
    %imagesc(atand((T_22+T_33)./T_11))
    imagesc((T_22+T_33)./T_11)
    Plotsetting_GOM1([0 1])
    plot_para('Filename','output/para_dielectric', 'Maximize',true)
%%
figure
    imagesc(sqrt(atand((T_22+T_33)./T_11)))
    Plotsetting_GOM1([0 10])
    plot_para('Filename','output/para_braggalpha', 'Maximize',true)
%%
figure
    imagesc(10*log10(conv2(vv_vv./hh_hh, ones(3,3), 'same')/9))
    Plotsetting_GOM1([0 15],1)


%% Incident angle and Bragg wavenumber
global im_size
r_phi = atand(linspace(4602.29004, 4602.29004+22490.8262, im_size(1))/12497);
figure
    plot(r_phi,'k', 'Linewidth',3)
    xlim([1 3300])
    ylim([20 65])
    set(gca,'XTick',0:600:3300, 'Xticklabel',cellstr(int2str((0:600*5:3300*5)'/1000))')
    hold on
    scatter(500,r_phi(500),100,'o','r','filled')
    scatter(1600,r_phi(1600),100,'o','r','filled')
    scatter(3100,r_phi(3100),100,'o','r','filled')
    hold off
    xlabel('Range (km)')
    ylabel('$\phi$ (deg)','Interpreter', 'latex')
    grid on
    pos = get(gca, 'Position');
    yoffset = +0.05;
    pos(2) = pos(2) + yoffset;
    set(gca, 'Position', pos)
    plot_para('Filename','incidence_angle', 'Maximize',true, 'Ratio',[4 3 1])
%% Bragg scattering coefficient 
epsilon_oil = 2.3-1j*0.02;
epsilon_sea = 70-1j*50;
k = 2*pi/(physconst('LightSpeed')/1.2575e9);
phi = 0:1:89;
theta = 15.2;
zeta = theta;
phi_i = acosd(cosd(phi+theta)*cosd(zeta));
B_hh_oil = (cosd(phi)-sqrt(epsilon_oil-sind(phi).^2))./(cosd(phi)+sqrt(epsilon_oil-sind(phi).^2));
B_vv_oil = ((epsilon_oil-1)*(sind(phi).^2 - epsilon_oil*(1 + sind(phi).^2)))./(epsilon_oil*cosd(phi)+sqrt(epsilon_oil-sind(phi).^2)).^2;
B_hh_sea = (cosd(phi)-sqrt(epsilon_sea-sind(phi).^2))./(cosd(phi)+sqrt(epsilon_sea-sind(phi).^2));
B_vv_sea = ((epsilon_sea-1)*(sind(phi).^2 - epsilon_sea*(1 + sind(phi).^2)))./(epsilon_sea*cosd(phi)+sqrt(epsilon_sea-sind(phi).^2)).^2;

S_hh_oil = cosd(phi_i).^4.* abs((sind(phi+theta)*cosd(zeta)./sind(phi_i)).^2.*B_hh_oil + (sind(zeta)./sind(phi_i)).^2.*B_vv_oil).^2;
S_vv_oil = cosd(phi_i).^4.* abs((sind(phi+theta)*cosd(zeta)./sind(phi_i)).^2.*B_vv_oil + (sind(zeta)./sind(phi_i)).^2.*B_hh_oil).^2;
S_hh_sea = cosd(phi_i).^4.* abs((sind(phi+theta)*cosd(zeta)./sind(phi_i)).^2.*B_hh_sea + (sind(zeta)./sind(phi_i)).^2.*B_vv_sea).^2;
S_vv_sea = cosd(phi_i).^4.* abs((sind(phi+theta)*cosd(zeta)./sind(phi_i)).^2.*B_vv_sea + (sind(zeta)./sind(phi_i)).^2.*B_hh_sea).^2;
figure
    a = plot(phi, abs(S_hh_sea),'b--', phi, abs(S_vv_sea),'b','Linewidth',3);
    hold on 
    b = plot(phi, abs(S_hh_oil),'k-.', phi, abs(S_vv_oil),'k:','Linewidth',3);
    grid on
    hold off
    
	set(gca,'Xlim',[0 80],'Ylim',[0,1])
    xlabel('$\phi$ (deg)','Interpreter', 'latex')
    ylabel('$S_{vv}$, $S_{hh}$','Interpreter', 'latex')
    plot_para('Filename','die_angle', 'Maximize',true, 'Ratio',[4 3 1])

figure
    plot(phi, 20*log10(abs(S_vv_sea./S_hh_sea)),'k--',......
    phi, 20*log10(abs(S_vv_oil./S_hh_oil)),'k','Linewidth',2.5)
    xlabel('$\phi$ (deg)','Interpreter', 'latex')
    ylabel('$|S_{vv}|^2 / |S_{hh}|^2$ (dB)','Interpreter', 'latex')
    set(gca,'Xlim',[0 80],'Ylim',[0 30],'xGrid','on','yGrid','on')
    plot_para('Filename','hv_ratio_thm', 'Maximize',true, 'Ratio',[4 3 1])
%%
figure
    plot(phi,atand(abs((B_hh-B_vv)./(B_hh+B_vv))),'k')
    grid on
    xlabel('angle (deg)')
%%  
figure
    plot(phi,2*pi./((physconst('LightSpeed')/1.2575e9)/2./sind(phi)), 'k', 'Linewidth',3)
    grid on
    xlabel('$\phi$ (deg)','Interpreter', 'latex')
    ylabel('$(2 \pi)/\Lambda$ (rad/m)','Interpreter', 'latex')
    xlim([20, 65])
    plot_para('Filename','braggwavenum', 'Maximize',true, 'Ratio',[4 3 1])