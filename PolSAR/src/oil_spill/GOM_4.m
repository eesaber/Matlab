clear 
clc
chk_pw('output_GulfMexico4/')
%% read data
disp('loading data...')
sub_map = 0;
if sub_map == false
    %[hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',10,'ReadNewFIle',true,'Type','mlc');
    [hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',10);
else
    %[hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',9,'Test','area1');
end
[N_az, N_ra] = size(hh_hh);
size_N = numel(hh_hh);
span = hh_hh+vv_vv+2*hv_hv;

%% 4-component decomposition
FourComp_decomp(hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv, '4decomp')
%% Power 
pow_range = [-35 5];
figure
    imagesc(10*log10(span))
    Plotsetting_GOM4(pow_range, 1, 'Colorbar_unit',[40 -70])
    colormap gray
    plot_para('Filename','span', 'Maximize',true)
%%
f_name = 'hh_hh';
figure
    imagesc(10*log10(hh_hh))
    Plotsetting_GOM4(pow_range, 1, 'Colorbar_unit',[40 -70])
    colormap gray
    plot_para('Filename',f_name, 'Maximize',true, 'Ratio',[4 3 1])
%%
f_name = 'vv_vv';
figure
    imagesc(10*log10(vv_vv))
    Plotsetting_GOM4(pow_range, 1,'Colorbar_unit',[40 -70])
    colormap gray
    plot_para('Filename',f_name, 'Maximize',true, 'Ratio',[4 3 1])
%%
figure
    imagesc(10*log10(hv_hv))
    Plotsetting_GOM4(pow_range, 1,'Colorbar_unit',[40 -70])
    colormap gray
    plot_para('Filename','hv_hv', 'Maximize',true, 'Ratio',[4 3 1])
%% Pauli decomposition
figure
Pauli_decomp((hh_hhh+vv_vv-hh_vv-conj(hh_vv)), 2*hv_hv,(hh_hh+vv_vv+hh_vv+conj(hh_vv)),'Filename','Pauli_decomp','saibu',false)
%%
mask = ones(3,3)/9;
T_11 = conv2((hh_hh+vv_vv+hh_vv+conj(hh_vv))/2, mask, 'same');
T_22 = conv2((hh_hh+vv_vv-hh_vv-conj(hh_vv))/2, mask, 'same');
T_33 = conv2(2*hv_hv, mask, 'same');
T_12 = conv2((hh_hh-vv_vv-hh_vv+conj(hh_vv))/2, mask, 'same');
T_13 = conv2(hh_hv + conj(hv_vv), mask, 'same');
T_23 = conv2(hh_hv - conj(hv_vv), mask, 'same');
clear  hh_hh hv_hv vv_vv hh_hv hh_vv hv_vv
close all
temp_T = cat(1,cat(2, reshape(T_11,[1,1,size_N]), reshape(T_12,[1,1,size_N]), reshape(T_13,[1,1,size_N])), ......
             cat(2,reshape(conj(T_12),[1,1,size_N]), reshape(T_22,[1,1,size_N]), reshape(T_23,[1,1,size_N])),.......
             cat(2,reshape(conj(T_13),[1,1,size_N]), reshape(conj(T_23),[1,1,size_N]), reshape(T_33,[1,1,size_N])));
clear T_11 T_22 T_33 T_12 T_13 T_23 
%% Eigen-decomposition
read = 0;
if read
    f_name = 'eigen_avg';
    [~] = Eigen_decomp('FileName', f_name);
else
    f_name = 'eigen_avg';
    [~] = Eigen_decomp('T',temp_T,'Calculate',true,'FileName',f_name);
end

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