clear 
clc
chk_pw('output_grandmesa/')
%% read data
disp('loading data...')
sub_map = 0;
if sub_map == 0
    [hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',4,'ReadNewFile',true);
else
    [hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',4,'Test','area1');
end
[N_az, N_ra] = size(hh_hh);
pow_range = [-35 5];
size_N = numel(hh_hh);
span = hh_hh+vv_vv+2*hv_hv;
mask = @(x) ones(x,x)/x^2;
mask_size = 9;
figure
Pauli_decomp((hh_hh+vv_vv-hh_vv-conj(hh_vv))/2, 2*hv_hv,(hh_hh+vv_vv+hh_vv+conj(hh_vv))/2,'Filename','Pauli_decomp','saibu',false)
%% Span
figure
    imagesc(10*log10(span))
    caxis(pow_range)
    colormap gray
%%
figure
    imagesc(10*log10(hh_hh))
    caxis(pow_range)
    colormap gray
    plot_para('Ratio',[4 3 1])
%%
f_name = 'vv_vv';
figure
    imagesc(10*log10(vv_vv))
    caxis(pow_range)
    colormap gray    
%%
figure
    imagesc(10*log10(conv2(hv_hv, mask(mask_size), 'same')))
    caxis(pow_range)
    colormap gray
%%
figure
    imagesc(10*log10(conv2(vv_vv./hh_hh, mask(9), 'same')))
    caxis([0 10])