% This file implements the detection method of oil-spill by using PolSAR.
% Target area is in the Aso, Japan, mission 1
clear 
clc
%% read data
disp('loading data...')
[hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',1);
[N_az, N_ra] = size(hh_hh);
size_N = numel(hh_hh);
span = hh_hh+vv_vv+2*hv_hv;
% 4-component decomposition
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
    %Pauli_decomp(2*T_22, T_33, 2*T_11, 'Filename','Pauli_decomp','saibu',false)
    Pauli_decomp((hh_hh+vv_vv-hh_vv-conj(hh_vv)), 2*hv_hv,(hh_hh+vv_vv+hh_vv+conj(hh_vv)),'Filename','Pauli_decomp')
close all
%% Eigen-decomposition
[H, alpha_bar] = Eigen_decomp('FileName', 'eigen');
close all