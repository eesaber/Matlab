% This file implements the detection method of oil-spill by using PolSAR.
% Target area is in the Gulf of Mexico, mission 9
clear 
clc
chk_pw('output_NorthSea1/')
%% read data
disp('loading data...')
sub_map = 0;
if sub_map == 0
    [hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',11,'ReadNewFile',true,'Type','mlc');
else
    [hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',11,'Test','area1');
end
[N_az, N_ra] = size(hh_hh);
pow_range = [-35 5];
size_N = numel(hh_hh);
span = hh_hh+vv_vv+2*hv_hv;

%% 4-component decomposition
plot_f = @ Plotsetting_NorthSea1;
FourComp_decomp(hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv, '4decomp','Powrange',[-35 5], 'Plotsetting', plot_f)
%% Span
figure
    imagesc(10*log10(span))
    Plotsetting_NorthSea1(pow_range, 1,'Colorbar_unit',"(dB)")
    colormap gray
    plot_para('Filename','span', 'Maximize',true)
    %Section_GOM2(10*log10(span./(trimmean(span,3,2)*ones(1,8000))),x,y,'Span($\bar{\bar{S}}) / \mu$ (dB)' ,'span')
%%
f_name = 'S_hh';
figure
    imagesc(10*log10(hh_hh))
    Plotsetting_NorthSea1(pow_range, 1, 'Colorbar_unit',"(dB)")
    colormap gray
    
    %Section_GOM2(10*log10(hh_hh./(trimmean(hh_hh,3,2)*ones(1,8000))),x,y,'$|S_{hh}|^2/ \mu $ (dB)','hh_hh')
    plot_para('Filename',f_name, 'Maximize',true, 'Ratio',[4 3 1])
%%
f_name = 'S_vv';
figure
    imagesc(10*log10(vv_vv))
    Plotsetting_NorthSea1(pow_range, 1,'Colorbar_unit',"(dB)")
    colormap gray
    plot_para('Maximize',true,'Filename',f_name,'Ratio',[4 3 1])
    
%%
f_name = 'S_hv';
figure
    imagesc(10*log10(hv_hv))
    Plotsetting_NorthSea1(pow_range, 1,'Colorbar_unit',"(dB)")
    colormap gray
    plot_para('Filename',f_name, 'Maximize',true, 'Ratio',[4 3 1])

%% Pauli decomposition
figure
Pauli_decomp((hh_hh+vv_vv-hh_vv-conj(hh_vv)), 2*hv_hv,(hh_hh+vv_vv+hh_vv+conj(hh_vv)),'Filename','Pauli_decomp','saibu',false)
%% Use moving Average T_mn 
mask = ones(9,9)/81;
T_11 = conv2((hh_hh+vv_vv+hh_vv+conj(hh_vv))/2, mask, 'same');
T_22 = conv2((hh_hh+vv_vv-hh_vv-conj(hh_vv))/2, mask, 'same');
T_33 = conv2(2*hv_hv, mask, 'same');
T_12 = conv2((hh_hh-vv_vv-hh_vv+conj(hh_vv))/2, mask, 'same');
T_13 = conv2(hh_hv + conj(hv_vv), mask, 'same');
T_23 = conv2(hh_hv - conj(hv_vv), mask, 'same');
%clear  hh_hh hv_hv vv_vv hh_hv hh_vv hv_vv
%%
temp_T = cat(1,cat(2, reshape(T_11,[1,1,size_N]), reshape(T_12,[1,1,size_N]), reshape(T_13,[1,1,size_N])), ......
             cat(2,reshape(conj(T_12),[1,1,size_N]), reshape(T_22,[1,1,size_N]), reshape(T_23,[1,1,size_N])),.......
             cat(2,reshape(conj(T_13),[1,1,size_N]), reshape(conj(T_23),[1,1,size_N]), reshape(T_33,[1,1,size_N])));
clear T_11 T_22 T_33 T_12 T_13 T_23 
%% Eigen-decomposition
read = 0;
f_name = 'eigen_avg2';
plot_set = @ Plotsetting_NorthSea1;
if read
    [~,~,A] = Eigen_decomp('Filename', f_name, 'Plotsetting', plot_set);
else
    [~,~,A] = Eigen_decomp('T',temp_T,'Calculate',true,'Filename',f_name, 'Plotsetting', plot_set);
end
%% Obtain terrain slope in azimuth and induced angle by terrain slope. 
%get_polangle(temp_T);
%clear temp_T
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
    Plotsetting_NorthSea1([0 3])
    xlabel('Azimuth (km)')
    ylabel('Range (km)') 
    plot_para('Filename','para_F', 'Maximize',true)
%}
%%
close all
f_name = 'para_confomty';
mu = 2*(-hv_hv+real(hh_vv))./(hh_hh + 2*hv_hv + vv_vv);
figure
%imagesc(mu)
imagesc(conv2(mu,ones(3,3)/9,'same'))
Plotsetting_NorthSea1([-0.4 1],1)
Section_GOM2(mu./(trimmean(mu,3,2)*ones(1,8000)),x,y,'$\mu_c$', f_name)
plot_para('Ratio',[4 3 1])

%delete(findall(gcf,'type','annotation'))
y_ano = [0.9 0.85; 0.5 0.55];
annotation('arrow',0.425*ones(1,2), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','k','Linewidth',3)
annotation('arrow',0.49*ones(1,2),  y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','r','Linewidth',3)    
annotation('arrow',0.521*ones(1,2), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','b','Linewidth',3)

annotation('arrow',0.619*ones(2,1), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','k','Linewidth',3)
annotation('arrow',0.679*ones(2,1), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','r','Linewidth',3)
annotation('arrow',0.686*ones(2,1), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','b','Linewidth',3)
plot_para('Filename', [f_name '_am'], 'Maximize',true,'Ratio',[4 3 1])
clear mu
%%
figure
    imagesc(-(hv_hv>abs(real(hh_vv))))
    Plotsetting_NorthSea1([-1 0],1)
    colormap gray; colorbar off
    
    plot_para('Filename', 'para_muller', 'Maximize',true)
%%
close all
f_name = 'para_corelation12';
figure
    imagesc(abs(T_12)./sqrt(T_11.*T_22))
    Plotsetting_NorthSea1([0 1],1)
    Section_GOM2(abs(T_12)./sqrt(T_11.*T_22),x,y,'$\gamma$', f_name)
    plot_para('Maximize',true,'Filename',[f_name '_am'],'Ratio',[4 3 1])
%%
f_name = 'para_moisture';
figure
    imagesc((T_11-T_22)./(T_11+T_22))
    Plotsetting_NorthSea1([0 1],1)
    Section_GOM2((T_11-T_22)./(T_11+T_22),x,y,f_name)   
%%
close all
f_name = 'para_roughness';
figure
    imagesc(conv2((T_22-T_33)./(T_22+T_33), ones(3,3)/9,'same'))
    Plotsetting_NorthSea1([0 1],1)
    Section_GOM2((T_22-T_33)./(T_22+T_33),x,y, '$\gamma_{r r_\perp}$',f_name)
    plot_para('Maximize',true,'Filename',[f_name '_am'],'Ratio',[4 3 1])
%%
f_name = 'para_dielectric';
figure
    %imagesc(atand((T_22+T_33)./T_11))
    imagesc((T_22+T_33)./T_11)
    Plotsetting_NorthSea1([0 1],1)
    Section_GOM2((T_22+T_33)./T_11,x,y, '$m (\phi, \epsilon_r)$',f_name)
    plot_para('Ratio',[4 3 1])
    
    y_ano = [0.9 0.85; 0.55 0.5; 0.33 0.28];
    annotation('arrow',0.425*ones(1,2), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','k','Linewidth',3)
    annotation('arrow',0.49*ones(1,2),  y_ano(2,:),'Headwidth',10,'Headstyle','vback2','color','r','Linewidth',3)    
    annotation('arrow',0.521*ones(1,2), y_ano(3,:),'Headwidth',10,'Headstyle','vback2','color','b','Linewidth',3)
    
    annotation('arrow',0.619*ones(2,1), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','k','Linewidth',3)
    annotation('arrow',0.679*ones(2,1), y_ano(2,:),'Headwidth',10,'Headstyle','vback2','color','r','Linewidth',3)
    annotation('arrow',0.686*ones(2,1), y_ano(3,:),'Headwidth',10,'Headstyle','vback2','color','b','Linewidth',3)
    
    plot_para('Maximize',true,'Filename',[f_name '_am'],'Ratio',[4 3 1])
%%
figure
    imagesc(sqrt(atand((T_22+T_33)./T_11)))
    Plotsetting_NorthSea1([0 10],1)
    Section_GOM2(sqrt(atand((T_22+T_33)./T_11)),x,y)
    %plot_para('Filename','para_braggalpha', 'Maximize',true)
%%
f_name = 'hvratio';
figure 
    imagesc(10*log10(conv2(vv_vv./hh_hh, ones(9,9)/81, 'same')))
    %imagesc(vv_vv./hh_hh)
    Plotsetting_NorthSea1([0 10],1)
    plot_para('Maximize',true,'Filename',f_name,'Ratio',[4 3 1])
