%%
clc
clear
%%
plotSetting = @ Plotsetting_dummy;
po = '/media/akb/2026EF9426EF696C/raw_data/20090811/ALPSRP188921500-L1.1/T3';
pi = '/home/akb/Code/Matlab/PolSAR/output';
x = PolSAR_AnalyzeTool([624 4608],'ALOS PALSAR',plotSetting,...
    'inputDataDir', po,  'outputDataDir', pi);
%% Radar-cross section
x.RCS('vv',true,'Span',true,'hh',true,'hv',true)
%% Noise equivalent sigam zero (NESZ)
NESZ = -29; % (dB)
x.isSNR(NESZ,'vv',true,'hh',true,'hv',true)
%% Image Analysis
close all
size_N = numel(hh_hh);
% Pauli decomposition
figure
Pauli_decomp((hh_hh+vv_vv-hh_vv-conj(hh_vv))/2, 2*hv_hv,(hh_hh+vv_vv+hh_vv+conj(hh_vv))/2,'Filename','Pauli_decomp','saibu',false)
% 4-component decomposition
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
[kai_1, kai_2, kai_3] = Logcumulant(temp_T);
%%
figure
image(cat(3, reshape(Stat(kai_1),im_size), reshape(Stat(kai_2), im_size), reshape(Stat(kai_3), im_size)))
set(gca,'Ydir','normal')
%set(gca,,'xlim',[1 15000],'XTickLabel',strsplit(num2str(get(gca, 'Xtick'))))
plot_para('Maximize',true,'Filename','logcumulant')
%%
close all
figure
scatter3(kai_1(:), kai_2(:), kai_3(:),'filled')
xlabel('\kappa_1')
ylabel('\kappa_2')
zlabel('\kappa_3')
figure
scatter(kai_3(:), kai_2(:),'filled')
xlabel('\kappa_3')
ylabel('\kappa_2')

figure
scatter(kai_2(:), kai_1(:),'filled')
xlabel('\kappa_2')
ylabel('\kappa_1')

figure
    histogram2(reshape(kai_3,[1 numel(kai_3)]), reshape(kai_2,[1 numel(kai_2)])......
        ,'DisplayStyle','tile','YBinLimits',[0 1],'XBinLimits',[-2 2])
%{
%% Statistic analysis
sigma_hh = hh_hh(:,range);
sigma_vv = vv_vv(:,range);
% Histogram
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
%}
%% Eigen-decomposition
read = 1;
f_name = 'eigen_9x9avg';
if read
    [H, alpha] = Eigen_decomp('Filename', f_name, 'Plotsetting', plot_set);
else
    temp_T = cat(1,cat(2, reshape(T_11,[1,1,size_N]), reshape(T_12,[1,1,size_N]), reshape(T_13,[1,1,size_N])), ......
             cat(2,reshape(conj(T_12),[1,1,size_N]), reshape(T_22,[1,1,size_N]), reshape(T_23,[1,1,size_N])),.......
             cat(2,reshape(conj(T_13),[1,1,size_N]), reshape(conj(T_23),[1,1,size_N]), reshape(T_33,[1,1,size_N])));
    [H, alpha] = Eigen_decomp('T',temp_T,'Calculate',true,'Filename',f_name, 'Plotsetting', plot_set);
end
%% Histogram of H/alpha
figure
    histogram2(reshape(H, [1 numel(H)]), reshape(alpha, [1 numel(alpha)])......
    ,'DisplayStyle','tile','YBinLimits',[0 90],'XBinLimits',[0 1],'NumBins',[900 2700]);
    colormap jet 
    set(gca,'Color',[0 0 0]+0.7,'XGrid','off','YGrid','off')
    hold on 
    H_Alpha([],'Seglinestyle','w')
    hold off	
    
%% gamma_12
f_name = 'para_gamma12';
figure
    imagesc(conv2(abs(T_12)./sqrt(T_11.*T_22), mask(9), 'same'))
    Plotsetting_dummy([0 1])
    plot_para('Maximize',true,'Filename',f_name)
%% σ_vv/ σ_hh
f_name = 'para_hvratio_g';
figure
    imagesc(10*log10(conv2(vv_vv./hh_hh, mask(9), 'same')))
    Plotsetting_dummy([-3 3])
    plot_para('Maximize',true,'Filename',f_name)
%% meteorology
fid = fopen('/media/akb/2026EF9426EF696C/raw_data/20090811/eng-hourly-08012009-08312009.csv','r');
C = textscan(fid,'%q %q %q %q %q %q %*[^\n]',320,'HeaderLines',16,'Delimiter',',');
fclose(fid);
T = cellfun(@str2double, C{6});
D = strcat(C{4},{'-Aug-'},C{2},{' '},C{5});
close all

dt = datetime(D,'InputFormat','dd-MMM-yyyy HH:mm','TimeZone','America/Cambridge_Bay');
dt = tzoffset(dt);
Day = datetime(D,'InputFormat','dd-MMM-yyyy HH:mm') -dt;
figure
plot(Day,T,'k','LineWidth',2)
ylabel('T $(^\circ)$','interpreter','latex')
grid on
xlim([datetime('Aug 03, 2009, 00:00:00','InputFormat','MMM dd, yyyy, hh:mm:ss') ...
	datetime('Aug 12, 2009, 00:00:00','InputFormat','MMM dd, yyyy, hh:mm:ss')])
xtickformat('dd')
plot_para('Filename','meteorology2','Maximize',true)