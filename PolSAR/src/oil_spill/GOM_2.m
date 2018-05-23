% This file implements the detection method of oil-spill by using PolSAR.
% Target area is in the Gulf of Mexico, mission 9
clear 
clc
chk_pw('output_GulfMexico2/area1/')
%% read data
disp('loading data...')
sub_map = 1;
if sub_map == 0
    [hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',9);
else
    [hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',9,'Test','area1');
end
[N_az, N_ra] = size(hh_hh);
size_N = numel(hh_hh);
span = hh_hh+vv_vv+2*hv_hv;
% 4-component decomposition
FourComp_decomp(hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv, '4decomp')
%% Span
figure
    imagesc(10*log10(span))
    Plotsetting_GOM2([-40 0],1,'Colorbar_unit',[40 -70])
    annotation('rectangle',[0.3 0.15 0.04 0.775],'Color','k','Linewidth',2)
        annotation('textbox',[0.25 0.5 0.1 0.1],'String','$A_1$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.4 0.15 0.18 0.775],'Color','k','Linewidth',2)
        annotation('textbox',[0.6 0.5 0.1 0.1],'String','$A_2$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    plot_para('Filename','span', 'Maximize',true)
figure
    imagesc(10*log10(hh_hh))
    Plotsetting_GOM2([-40 0],1,'Colorbar_unit',[40 -70])
    plot_para('Filename','S_hh', 'Maximize',true, 'Ratio',[4 3 1])
figure
    imagesc(10*log10(vv_vv))
    Plotsetting_GOM2([-40 0],1,'Colorbar_unit',[40 -70])
    plot_para('Filename','S_vv', 'Maximize',true, 'Ratio',[4 3 1])
figure
    imagesc(10*log10(hv_hv))
    Plotsetting_GOM2([-40 0],1,'Colorbar_unit',[40 -70])
    plot_para('Filename','S_hv', 'Maximize',true, 'Ratio',[4 3 1])

%% Pauli decomposition
figure
    %Pauli_decomp(2*T_22, T_33, 2*T_11, 'Filename','Pauli_decomp','saibu',false)
    Pauli_decomp((hh_hh+vv_vv-hh_vv-conj(hh_vv)), 2*hv_hv,(hh_hh+vv_vv+hh_vv+conj(hh_vv)),'Filename','Pauli_decomp')
close all
%% Eigen-decomposition
if sub_map == 0
    [H, alpha_bar] = Eigen_decomp('FileName','eigen');
else
    [H, alpha_bar] = Eigen_decomp('FileName','eigen_area1');
end
%% Obtain terrain slope in azimuth and induced angle by terrain slope. 
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
clear T_11 T_22 T_33 T_12 T_13 T_23 
get_polangle(temp_T);
clear temp_T
T_11 = (hh_hh+vv_vv+hh_vv+conj(hh_vv))/2;
T_22 = (hh_hh+vv_vv-hh_vv-conj(hh_vv))/2;
T_33 = 2*hv_hv;
T_12 = (hh_hh-vv_vv-hh_vv+conj(hh_vv))/2;
T_13 = hh_hv + conj(hv_vv);
T_23 = hh_hv - conj(hv_vv);

%% Indicator
%{
figure
    imagesc(H+alpha_bar/180*pi+A_1+abs(hh_vv)./sqrt(hh_hh.*vv_vv))
    Plotsetting_GOM2([0 3])
    xlabel('Azimuth (km)')
    ylabel('Range (km)') 
    plot_para('Filename','para_F', 'Maximize',true)
%}
%%
figure
    imagesc(2*(hv_hv-real(hh_vv))./(hh_hh + 2*hv_hv + vv_vv))
    Plotsetting_GOM2([-1 1],1)
    annotation('rectangle',[0.3 0.15 0.04 0.775],'Color','k','Linewidth',2)
        annotation('textbox',[0.25 0.5 0.1 0.1],'String','$A_{\mu 1}$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.4 0.15 0.18 0.775],'Color','k','Linewidth',2)
        annotation('textbox',[0.6 0.5 0.1 0.1],'String','$A_{\mu 2}$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    plot_para('Filename','para_confomty', 'Maximize',true)
%%
figure
    imagesc(-(hv_hv>abs(real(hh_vv))))
    Plotsetting_GOM2([-1 0],1)
    annotation('rectangle',[0.13 0.8 0.05 0.125],'Color','k','Linewidth',2)
        annotation('textbox',[0.18 0.8 0.1 0.1],'String','$A_{k 1}$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.57 0.45 0.05 0.2],'Color','k','Linewidth',2)
        annotation('textbox',[0.5 0.55 0.1 0.1],'String','$A_{k 2}$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.4 0.8 0.18 0.125],'Color','k','Linewidth',2)
        annotation('textbox',[0.6 0.8 0.1 0.1],'String','$A_{k 3}$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    colormap gray; colorbar off    
    plot_para('Filename',para_muller', 'Maximize',true)
%%
figure
    imagesc(abs(T_12)./sqrt(T_11.*T_22))
    Plotsetting_GOM2([0 1],1)
    annotation('rectangle',[0.3 0.45 0.04 0.475],'Color','k','Linewidth',2)
        annotation('textbox',[0.25 0.6 0.1 0.1],'String','$A_{\gamma 1}$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.4 0.45 0.18 0.475],'Color','k','Linewidth',2)
        annotation('textbox',[0.6 0.6 0.1 0.1],'String','$A_{\gamma 2}$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    plot_para('Filename','para_corelation12', 'Maximize',true)
%%
figure
    imagesc((T_11-T_22)./(T_11+T_22))
    Plotsetting_GOM2([0 1],1)
    annotation('rectangle',[0.44 0.3 0.13 0.325],'Color','k','Linewidth',2)
        annotation('textbox',[0.6 0.5 0.1 0.1],'String','$A_{m 1}$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    plot_para('Filename','para_moisture', 'Maximize',true)
%%
figure
    imagesc((T_22-T_33)./(T_22+T_33))
    Plotsetting_GOM2([0 1],1)
    annotation('rectangle',[0.45 0.17 0.14 0.275],'Color','k','Linewidth',2)
        annotation('textbox',[0.4 0.2 0.1 0.1],'String','$A_1$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    annotation('rectangle',[0.4 0.65 0.15 0.275],'Color','k','Linewidth',2)
        annotation('textbox',[0.55 0.75 0.1 0.1],'String','$A_2$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    plot_para('Filename','para_roughness', 'Maximize',true)
%%
figure
    %imagesc(atand((T_22+T_33)./T_11))
    imagesc((T_22+T_33)./T_11)
    Plotsetting_GOM2([0 1],1)
    annotation('rectangle',[0.4 0.25 0.18 0.675],'Color','k','Linewidth',2)
        annotation('textbox',[0.6 0.5 0.1 0.1],'String','$A_{\epsilon 1}$','Linestyle','none','Fontsize',40,'Interpreter', 'latex')
    plot_para('Filename','para_dielectric', 'Maximize',true)
%%
figure
    imagesc(sqrt(atand((T_22+T_33)./T_11)))
    Plotsetting_GOM2([0 10],1)
    plot_para('Filename','para_braggalpha', 'Maximize',true)