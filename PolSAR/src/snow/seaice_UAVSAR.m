%% Parameters setting
clear 
clc
chk_pw('output_Beaufort/')
mission_num = 7;
pow_range = [-35 5];
mask = @(x) ones(x,x)/x^2;
mask_size = 9;
plot_set = @ Plotsetting_dummy;
global im_size
%% read data
disp('loading data...')
sub_map = 0;
if sub_map == 0
    %[hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',mission_num,'ReadNewFile',true,'Type','mlc');
    [hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',mission_num);
else
    [hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',mission_num,'Test','area1');
end
%% Radar-cross section
span = hh_hh+vv_vv+2*hv_hv;
figure
    imagesc(10*log10(span))
    Plotsetting_dummy(pow_range)
    colormap gray
    plot_para('Maximize',true,'Filename','Span')
    %%
figure
    imagesc(10*log10(vv_vv))
    Plotsetting_dummy(pow_range)
    plot_para('Maximize',true,'Filename','sigma_vv')
    %%
figure
    imagesc(10*log10(hh_hh))
    Plotsetting_dummy(pow_range)
    plot_para('Maximize',true,'Filename','sigma_hh')
    %%
figure
    imagesc(10*log10(hv_hv))
    Plotsetting_dummy([-45 -20])
    plot_para('Maximize',true,'Filename','sigma_hv')
% Noise equivalent sigam zero (NESZ)
phi = atand((4512.83+linspace(0,22520.48,3300))/12496.8).';
NESZ = 1.9664e-2*phi.^2 -1.5561*phi - 24.0269; % (dB)
figure
    plot(phi, NESZ)
figure
    imagesc(10*log10(hh_hh)-repmat(NESZ,1,im_size(2))>6)
    title('\sigma_{hh}')
figure
    imagesc(10*log10(hv_hv)-repmat(NESZ,1,im_size(2))>6)
    title('\sigma_{hv}')
figure
    imagesc(10*log10(vv_vv)-repmat(NESZ,1,im_size(2))>6)
    title('\sigma_{vv}')
%% Image Analysis
size_N = numel(hh_hh);
%% Pauli decomposition
figure
Pauli_decomp((hh_hh+vv_vv-hh_vv-conj(hh_vv))/2, 2*hv_hv,(hh_hh+vv_vv+hh_vv+conj(hh_vv))/2,'Filename','Pauli_decomp','saibu',false)
%% 4-component decomposition
FourComp_decomp(hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv, '4decomp','Plotsetting', plot_set, 'Powrange',[-35 5])
%% Log cumulant
clear span
if 0
    temp_C = cat(1,cat(2, reshape(hh_hh,[1,1,size_N]), sqrt(2)*reshape(hh_hv,[1,1,size_N]), reshape(hh_vv,[1,1,size_N])), ......
             cat(2,reshape(sqrt(2)*conj(hh_hv),[1,1,size_N]), reshape(2*hv_hv,[1,1,size_N]), reshape(sqrt(2)*hv_vv,[1,1,size_N])),.......
             cat(2,reshape(conj(hh_vv),[1,1,size_N]), reshape(conj(sqrt(2)*hv_vv),[1,1,size_N]), reshape(vv_vv,[1,1,size_N])));
else
    temp_C = cat(1,cat(2, reshape(hh_hh,[1,1,size_N]), reshape(hh_vv,[1,1,size_N])), ......
             cat(2,reshape(conj(hh_vv),[1,1,size_N]), reshape(vv_vv,[1,1,size_N])));
end
clear hh_hh vv_vv hv_hv hh_hv hh_vv hv_vv
[kai_1, kai_2, kai_3] = Logcumulant(temp_C);
clear temp_C
%%
figure
image(cat(3, reshape(Stat(kai_1),im_size), reshape(Stat(kai_2), im_size), reshape(Stat(kai_3), im_size)))
set(gca,'Ydir','normal','xlim',[1 15000])
%%
set(gca,'Ydir','normal','xlim',[15001 30000])
set(gca,'XTickLabel',strsplit(num2str(get(gca, 'Xtick'))))
plot_para('Maximize',true,'Filename','logcumulant_2')
%% Segementation
range = 35001:45000;
temp = cat(3, reshape(Stat(kai_1), im_size), reshape(Stat(kai_2), im_size), reshape(Stat(kai_3), im_size));
temp = temp(:,range,:);
lab_I = rgb2lab(temp);
clear temp
ab = lab_I(:,:,2:3);
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,size(ab,1)* size(ab,2),2);
nColors = 4;
% repeat the clustering 3 times to avoid local minima
[cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean', ...
                                      'Replicates',3);
pixel_labels = uint8(reshape(cluster_idx,nrows,ncols));
clear ab
figure
    imshow(pixel_labels,[],'Border','tight')
%% Statistic analysis
[hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',mission_num); % Read image again
sigma_hh = hh_hh(:,range);
sigma_vv = vv_vv(:,range);
%% Histogram
bin_limit = [-6, 13];
figure
    histogram(10*log10(sigma_vv(pixel_labels==1)./sigma_hh(pixel_labels==1)),'BinLimits',bin_limit,'Normalization','probability')
    hold on 
    for n = 1 : nColors
        histogram(10*log10(sigma_vv(pixel_labels==n)./sigma_hh(pixel_labels==n)),'BinLimits',bin_limit,'Normalization','probability')
    end
    hold off
    legend('label 1','label 2','label 3','label 4')
    xlabel('$\sigma_{vv}/\sigma_{hh}$', 'interpreter', 'latex')
    ylabel('Pr$(x = \sigma_{vv}/\sigma_{hh})$', 'interpreter', 'latex')
    plot_para('FileName','histo_vhratio','Maximize',true)
%%
bin_limit = [1e-5, 0.025];
figure
    histogram(sigma_vv(pixel_labels==1),'BinLimits', bin_limit,'Normalization','probability')
    hold on 
    for n = nColors:-1:2
        histogram(sigma_vv(pixel_labels==n),'BinLimits', bin_limit,'Normalization','probability')
    end
    hold off
    legend('label 1','label 2','label 3','label 4')
    xlabel('$|S_{vv}|^2$', 'interpreter', 'latex')
    ylabel('Pr$(x = |S_{vv}|^2)$', 'interpreter', 'latex')
    plot_para('FileName','histo_vv','Maximize',true)
%%
figure
    imagesc(kai_3)
    set(gca,'Ydir','normal','clim',[-1.5 1.5],'xlim',[15001 30000])
    plot_para('Maximize',true,'Filename','kappa_3rd')
figure
    histogram2(reshape(kai_3,[1 numel(kai_3)]), reshape(kai_2,[1 numel(kai_2)])......
        ,'DisplayStyle','tile','YBinLimits',[0 1],'XBinLimits',[-2 2])
%% Eigen-decomposition
clear kai_1 kai_2 kai_3
read = 1;
f_name = 'eigen_9x9avg';
if read
    [H, alpha] = Eigen_decomp('Filename', f_name, 'Plotsetting', plot_set);
else
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
    [H, alpha] = Eigen_decomp('T',temp_T,'Calculate',true,'Filename',f_name, 'Plotsetting', plot_set);
    %[hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',mission_num); % Read image again
end
%% Histogram of H/alpha
figure
    histogram2(H(:), alpha(:), 'DisplayStyle','tile','YBinLimits',[0 90],...
        'XBinLimits',[0 1],'NumBins',[700 2100]);
    caxis([0 3000])
    colormap jet; colorbar
    set(gca,'Color',[0 0 0]+0.7,'XGrid','off','YGrid','off')
    hold on 
    H_Alpha([],'Seglinestyle','w')
    hold off	
    
%% Surface Analysis
clear H alpha
T_11 = (hh_hh+vv_vv+hh_vv+conj(hh_vv))/2;
T_22 = (hh_hh+vv_vv-hh_vv-conj(hh_vv))/2;
T_33 = 2*hv_hv;
T_12 = (hh_hh-vv_vv-hh_vv+conj(hh_vv))/2;
T_13 = hh_hv + conj(hv_vv);
T_23 = hh_hv - conj(hv_vv);
% sigma_vv/sigma_hh 
figure
imagesc(10*log10(conv2(vv_vv./hh_hh, mask(9), 'same')))
set(gca,'Ydir','normal','clim',[0 10])
colormap jet; colorbar
%% gamma_12
f_name = 'para_gamma12';
figure
    imagesc(abs(T_12)./sqrt(T_11.*T_22))
    Plotsetting_dummy([0 1])
    xlim([range(1) range(end)])
    plot_para('Maximize',true,'Filename',f_name)
%% σ_vv/ σ_hh
f_name = 'para_hvratio_g';
figure
    imagesc(conv2(vv_vv./hh_hh, mask(9), 'same'))
    Plotsetting_dummy([0 10])
    plot_para('Maximize',true,'Filename',f_name)
