%% The code implements PolSAR decomposition by using non-negative matrix factorization
clear 
clc
%% read data
disp('loading data...')
area = 'area1';
[hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO();
[N_az, N_ra] = size(hh_hh);
size_N = numel(hh_hh);
span = hh_hh+vv_vv+2*hv_hv;
% 4-component decomposition
FourComp_decomp(hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv, '4decomp')
%% Span
figure
    imagesc(10*log10(span))
    Plotsetting_1([-40 0],'Colorbar_unit',[40 -70])
    annotation('rectangle',[0.125 0.4 0.04 0.525],'Color','k','Linewidth',2)
        annotation('textbox',[0.17 0.71 0.1 0.1],'String','$A_1$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.78 0.2 0.06 0.725],'Color','k','Linewidth',2)
        annotation('textbox',[0.835 0.71 0.1 0.1],'String','$A_2$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.68 0.2 0.1 0.325],'Color','k','Linewidth',2)
        annotation('textbox',[0.66 0.71 0.1 0.1],'String','$A_3$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.7 0.525 0.06 0.4],'Color','k','Linewidth',2)
        annotation('textbox',[0.64 0.4 0.1 0.1],'String','$A_4$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    xlabel('Azimuth (km)')
    ylabel('Range (km)')
    plot_para('Filename','output/span', 'Maximize',true)
%%
T_11 = (hh_hh+vv_vv+hh_vv+conj(hh_vv))/2;
T_22 = (hh_hh+vv_vv-hh_vv-conj(hh_vv))/2;
T_33 = 2*hv_hv;
T_12 = (hh_hh-vv_vv-hh_vv+conj(hh_vv))/2;
T_13 = hh_hv + conj(hv_vv);
T_23 = hh_hv - conj(hv_vv);
%clear  hh_hh hv_hv vv_vv hh_hv hh_vv hv_vv
close all
temp_T = cat(1,cat(2, reshape(T_11,[1,1,size_N]), reshape(T_12,[1,1,size_N]), reshape(T_13,[1,1,size_N])), ......
             cat(2,reshape(conj(T_12),[1,1,size_N]), reshape(T_22,[1,1,size_N]), reshape(T_23,[1,1,size_N])),.......
             cat(2,reshape(conj(T_13),[1,1,size_N]), reshape(conj(T_23),[1,1,size_N]), reshape(T_33,[1,1,size_N])));
%% Pauli decomposition
figure
Pauli_decomp(2*T_22, T_33, 2*T_11, 'Filename','Pauli_decomp','saibu',false)
%Pauli_decomp((hh_hh+vv_vv-hh_vv-conj(hh_vv))/2, 2*hv_hv,(hh_hh+vv_vv+hh_vv+conj(hh_vv))/2,'Filename','Pauli_decomp')
close all
%% Eigen-decomposition
[H, alpha_bar] = Eigen_decomp(temp_T, span);
%Eigen_decomp(1,1);
close all
%% Indicator
%{
figure
    imagesc(H+alpha_bar/180*pi+A_1+abs(hh_vv)./sqrt(hh_hh.*vv_vv))
    Plotsetting_1([0 3])
    xlabel('Azimuth (km)')
    ylabel('Range (km)') 
    plot_para('Filename','output/para_F', 'Maximize',true)
%}
%%
figure
    imagesc(2*(hv_hv-real(hh_vv))./(hh_hh + 2*hv_hv + vv_vv))
    Plotsetting_1([-1 1])
    xlabel('Azimuth (km)')
    ylabel('Range (km)') 
    plot_para('Filename','output/para_confomty', 'Maximize',true)
%%
figure
    imagesc(-(hv_hv>abs(real(hh_vv))))
    Plotsetting_1([-1 0])
    colormap gray; colorbar off
    xlabel('Azimuth (km)')
    ylabel('Range (km)') 
    plot_para('Filename','output/para_muller', 'Maximize',true)
%%
figure
    imagesc(abs(T_12)./sqrt(T_11.*T_22))
    Plotsetting_1([0 1])
    xlabel('Azimuth (km)')
    ylabel('Range (km)') 
    plot_para('Filename','output/para_corelation12', 'Maximize',true)
%%
figure
    imagesc((T_11-T_22)./(T_11+T_22))
    Plotsetting_1([0 1])
    xlabel('Azimuth (km)')
    ylabel('Range (km)') 
    plot_para('Filename','output/para_moisture', 'Maximize',true)
%%
figure
    imagesc((T_22-T_33)./(T_22+T_33))
    Plotsetting_1([0 1])
    xlabel('Azimuth (km)')
    ylabel('Range (km)') 
    plot_para('Filename','output/para_roughness', 'Maximize',true)
%%
figure
    %imagesc(atand((T_22+T_33)./T_11))
    imagesc((T_22+T_33)./T_11)
    Plotsetting_1([0 1])
    xlabel('Azimuth (km)')
    ylabel('Range (km)') 
    plot_para('Filename','output/para_dielectric', 'Maximize',true)
%%
figure
    imagesc(sqrt(atand((T_22+T_33)./T_11)))
    Plotsetting_1([0 10])
    xlabel('Azimuth (km)')
    ylabel('Range (km)') 
    plot_para('Filename','output/para_braggalpha', 'Maximize',true)
%% Obtain terrain slope in azimuth and induced angle by terrain slope. 
close all
Find_angle(temp_T);
%clear T_11 T_22 T_33 T_12 T_13 T_23 
clear temp_T 
%% NMF
 
%% Incident angle and Bragg wavenumber
global im_size
figure
    plot(atand(linspace(4602.29004, 4602.29004+22490.8262, im_size(1))/12497),'k', 'Linewidth',3)
    xlim([1 3300])
    ylim([20 65])
    xlabel('Range')
    ylabel('incidence angle (deg)')
    grid on 
    hold on 
    plot([1500 1500], [20 ,50],'-.r', 'Linewidth',2)
    hold off
    annotation('textbox',[.3 .2 .3 .2],'String', 'small','FitBoxToText','on','FontSize', 40);
    annotation('textbox',[.6 .2 .6 .2],'String', 'big','FitBoxToText','on','FontSize', 40);
    plot_para('Filename','output/incidence_angle', 'Maximize',true)
%% Bragg scattering coefficient 
epsilon_oil = 2;
epsilon_sea = 40;
theta = 0:1:90;
beta = 2*pi/(physconst('LightSpeed')/1.2575e9);
B_hh_oil = (cosd(theta)-sqrt(epsilon_oil-sind(theta).^2))./(cosd(theta)+sqrt(epsilon_oil-sind(theta).^2));
B_vv_oil = ((epsilon_oil-1)*(sind(theta).^2 - epsilon_oil*(1 + sind(theta).^2)))./(epsilon_oil*cosd(theta)+sqrt(epsilon_oil-sind(theta).^2)).^2;
B_hh_sea = (cosd(theta)-sqrt(epsilon_sea-sind(theta).^2))./(cosd(theta)+sqrt(epsilon_sea-sind(theta).^2));
B_vv_sea = ((epsilon_sea-1)*(sind(theta).^2 - epsilon_sea*(1 + sind(theta).^2)))./(epsilon_sea*cosd(theta)+sqrt(epsilon_sea-sind(theta).^2)).^2;
figure
    plot(theta, abs(B_hh_sea),'b-.', theta, abs(B_vv_sea),'b','Linewidth',3)
    hold on 
    plot(theta, abs(B_hh_oil),'k-.', theta, abs(B_vv_oil),'k','Linewidth',3)
    plot([65, 65], [0,10],'r', 'Linewidth',3)
    plot([25, 25], [0,10],'r', 'Linewidth',3)
    grid on
    hold off
    ylim([0,10])
    xlim([0,90])
    xlabel('angle (deg)')
    plot_para('Filename','output/die_angle', 'Maximize',true)
%%
    figure
    plot(theta,atand(abs((B_hh-B_vv)./(B_hh+B_vv))),'k')
    grid on
    xlabel('angle (deg)')

figure
    plot(theta,(physconst('LightSpeed')/1.2575e9)/2./sind(theta))
    grid on
    xlabel('angle (deg)')