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

%% 4-component decomposition
FourComp_decomp(hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv, '4decomp')
%% Span
x = [1900, 4200, 5100];
y = [3100, 1600, 500];
pow_range = [-35 5];
figure
    imagesc(10*log10(span))
    Plotsetting_GOM2(pow_range, 1, 'Colorbar_unit',[40 -70])
    colormap gray
    %plot_para('Filename','span', 'Maximize',true)
    Section_GOM2(10*log10(span./(trimmean(span,3,2)*ones(1,8000))),x,y,'Span($\bar{\bar{S}}) / \mu$ (dB)' ,'span')
%%
figure
    imagesc(10*log10(hh_hh))
    Plotsetting_GOM2(pow_range, 1, 'Colorbar_unit',[40 -70])
    colormap gray
    plot_para('Ratio',[4 3 1])
    Section_GOM2(10*log10(hh_hh./(trimmean(hh_hh,3,2)*ones(1,8000))),x,y,'$|S_{hh}|^2/ \mu $ (dB)','hh_hh')
    %plot_para('Filename','S_hh', 'Maximize',true, 'Ratio',[4 3 1])
%%
figure
    imagesc(10*log10(vv_vv))
    Plotsetting_GOM2(pow_range, 1,'Colorbar_unit',[40 -70])
    colormap gray
    plot_para('Ratio',[4 3 1])
    Section_GOM2(10*log10(vv_vv./(trimmean(vv_vv,3,2)*ones(1,8000))),x,y,'$|S_{vv}|^2/ \mu $ (dB)','vv_vv')
    %plot_para('Filename','S_vv', 'Maximize',true, 'Ratio',[4 3 1])
%%
figure
    imagesc(10*log10(hv_hv))
    Plotsetting_GOM2(pow_range, 1,'Colorbar_unit',[40 -70])
    plot_para('Ratio',[4 3 1])
    colormap gray
    Section_GOM2(10*log10(hv_hv./(trimmean(hv_hv,3,2)*ones(1,8000))),x,y,'$|S_{hv}|^2$ (dB)','hv_hv')
    %plot_para('Filename','S_hv', 'Maximize',true, 'Ratio',[4 3 1])

%% Pauli decomposition
figure
Pauli_decomp((hh_hh+vv_vv-hh_vv-conj(hh_vv)), 2*hv_hv,(hh_hh+vv_vv+hh_vv+conj(hh_vv)),'Filename','Pauli_decomp','saibu',false)
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
read = 1;
if read
    f_name = 'eigen_avg';
    [~] = Eigen_decomp('FileName', f_name);
else
    f_name = 'eigen_avg';
    [~] = Eigen_decomp('T',temp_T,'Calculate',true,'FileName',f_name);
end
clear temp_T
close all
%% Obtain terrain slope in azimuth and induced angle by terrain slope. 
get_polangle(temp_T);
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
    Plotsetting_GOM2([0 3])
    xlabel('Azimuth (km)')
    ylabel('Range (km)') 
    plot_para('Filename','para_F', 'Maximize',true)
%}
%%
close all
f_name = 'para_confomty';
figure
    %imagesc(2*(real(hh_vv)-hv_hv)./(hh_hh + 2*hv_hv + vv_vv))
    imagesc(conv2(2*(hv_hv-real(hh_vv))./(hh_hh + 2*hv_hv + vv_vv),ones(3,3)/9,'same'))
    Plotsetting_GOM2([-1 1],1)
    %Section_GOM2(2*(hv_hv-real(hh_vv))./(hh_hh + 2*hv_hv + vv_vv),x,y, f_name)
    Section_GOM2((2*(hv_hv-real(hh_vv))./(hh_hh + 2*hv_hv + vv_vv))./(trimmean( (2*(hv_hv-real(hh_vv))./(hh_hh + 2*hv_hv + vv_vv)) ,3,2)*ones(1,8000)),x,y, f_name)    
%%
figure
    imagesc(-(hv_hv>abs(real(hh_vv))))
    Plotsetting_GOM2([-1 0],1)
    colormap gray; colorbar off
    
    plot_para('Filename', 'para_muller', 'Maximize',true)
%%
close all
f_name = 'para_corelation12';
figure
    imagesc(abs(T_12)./sqrt(T_11.*T_22))
    Plotsetting_GOM2([0 1],1)
    Section_GOM2(abs(T_12)./sqrt(T_11.*T_22),x,y,f_name)
%%
f_name = 'para_moisture';
figure
    imagesc((T_11-T_22)./(T_11+T_22))
    Plotsetting_GOM2([0 1],1)
    Section_GOM2((T_11-T_22)./(T_11+T_22),x,y,f_name)   
%%
close all
f_name = 'para_roughness';
figure
    imagesc((T_22-T_33)./(T_22+T_33))
    Plotsetting_GOM2([0 1],1)
    Section_GOM2((T_22-T_33)./(T_22+T_33),x,y, f_name)
%%

f_name = 'para_dielectric';
figure
    %imagesc(atand((T_22+T_33)./T_11))
    imagesc((T_22+T_33)./T_11)
    Plotsetting_GOM2([0 1],1)
    Section_GOM2((T_22+T_33)./T_11,x,y, '$m (\phi, \epsilon_r)$',f_name)
%%
figure
    imagesc(sqrt(atand((T_22+T_33)./T_11)))
    Plotsetting_GOM2([0 10],1)
    Section_GOM2(sqrt(atand((T_22+T_33)./T_11)),x,y)
    %plot_para('Filename','para_braggalpha', 'Maximize',true)
%%
f_name = 'hvratio';
figure 
    imagesc(10*log10(conv2(vv_vv./hh_hh, ones(3,3), 'same')/9))
    %imagesc(vv_vv./hh_hh)
    Plotsetting_GOM2([0 15],1)
    Section_GOM2(10*log10(conv2(vv_vv./hh_hh, ones(3,3), 'same')/9),x,y,'$|S_{vv}|^2 / |S_{hh}|^2$ (dB)',f_name)