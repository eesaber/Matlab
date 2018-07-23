%% Parameters setting
clear 
clc
chk_pw('output_Beaufort/')
mission_num = 7;
pow_range = [-35 5];
mask = @(x) ones(x,x)/x^2;
mask_size = 9;
plot_set = @ Plotsetting_dummy;
%% read data
disp('loading data...')
sub_map = 0;
if sub_map == 0
    %[hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',7,'ReadNewFile',true,'Type','mlc');
    [hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',mission_num);
else
    [hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',mission_num,'Test','area1');
end
%% Radar-cross section
span = hh_hh+vv_vv+2*hv_hv;
figure
    imagesc(10*log10(span))
    set(gca,'Ydir','normal', 'clim',pow_range)
    colormap jet; colorbar
    plot_para('Maximize',true,'Filename','Span')
figure
    imagesc(10*log10(vv_vv))
    set(gca,'Ydir','normal', 'clim',pow_range)
    colormap jet; colorbar
    plot_para('Maximize',true,'Filename','S_vv')
figure
    imagesc(10*log10(hh_hh))
    set(gca,'Ydir','normal', 'clim',pow_range)
    colormap jet; colorbar
    plot_para('Maximize',true,'Filename','S_hh')
figure
    imagesc(10*log10(hv_hv))
    set(gca,'Ydir','normal', 'clim',[-45 -20])
    colormap jet; colorbar
    plot_para('Maximize',true,'Filename','S_hv')
%% Image Analysis
size_N = numel(hh_hh);
%% Pauli decomposition
figure
Pauli_decomp((hh_hh+vv_vv-hh_vv-conj(hh_vv))/2, 2*hv_hv,(hh_hh+vv_vv+hh_vv+conj(hh_vv))/2,'Filename','Pauli_decomp','saibu',false)
%% 4-component decomposition
FourComp_decomp(hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv, '4decomp','Plotsetting', plot_set, 'Powrange',[-35 5])
%% Log-cumulant
if 0
    temp_C = cat(1,cat(2, reshape(hh_hh,[1,1,size_N]), sqrt(2)*reshape(hh_hv,[1,1,size_N]), reshape(hh_vv,[1,1,size_N])), ......
             cat(2,reshape(sqrt(2)*conj(hh_hv),[1,1,size_N]), reshape(2*hv_hv,[1,1,size_N]), reshape(sqrt(2)*hv_vv,[1,1,size_N])),.......
             cat(2,reshape(conj(hh_vv),[1,1,size_N]), reshape(conj(sqrt(2)*hv_vv),[1,1,size_N]), reshape(vv_vv,[1,1,size_N])));
else
    temp_C = cat(1,cat(2, reshape(hh_hh,[1,1,size_N]), reshape(hh_vv,[1,1,size_N])), ......
             cat(2,reshape(conj(hh_vv),[1,1,size_N]), reshape(vv_vv,[1,1,size_N])));
end
clear hh_hh vv_vv hv_hv hh_hv hh_vv hv_vv
[kai_1, kai_2] = Logcumulant(temp_C);
clear temp_C
[hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',mission_num); % Read image again

figure
image(cat(3,reshape(Stat(kai_1),size(hh_hh)), reshape(Stat(kai_2),size(hh_hh)), zeros(size(hh_hh))))
set(gca,'Ydir','normal')
%% Eigen-decomposition
temp_size = 9;
T_11 = conv2((hh_hh+vv_vv+hh_vv+conj(hh_vv))/2, mask(temp_size), 'same');
T_22 = conv2((hh_hh+vv_vv-hh_vv-conj(hh_vv))/2, mask(temp_size), 'same');
T_33 = conv2(2*hv_hv, mask(temp_size), 'same');
T_12 = conv2((hh_hh-vv_vv-hh_vv+conj(hh_vv))/2, mask(temp_size), 'same');
T_13 = conv2(hh_hv + conj(hv_vv), mask(temp_size), 'same');
T_23 = conv2(hh_hv - conj(hv_vv), mask(temp_size), 'same');
clear hh_hh vv_vv hv_hv hh_hv hh_vv hv_vv
temp_T = cat(1,cat(2, reshape(T_11,[1,1,size_N]), reshape(T_12,[1,1,size_N]), reshape(T_13,[1,1,size_N])), ......
             cat(2,reshape(conj(T_12),[1,1,size_N]), reshape(T_22,[1,1,size_N]), reshape(T_23,[1,1,size_N])),.......
             cat(2,reshape(conj(T_13),[1,1,size_N]), reshape(conj(T_23),[1,1,size_N]), reshape(T_33,[1,1,size_N])));
clear T_11 T_22 T_33 T_12 T_13 T_23 
read = 1;
f_name = 'eigen_9x9avg';
if read
    [~,lambda] = Eigen_decomp('Filename', f_name, 'Plotsetting', plot_set);
else
    [~,lambda] = Eigen_decomp('T',temp_T,'Calculate',true,'Filename',f_name, 'Plotsetting', plot_set);
end
[hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',mission_num); % Read image again
%%
