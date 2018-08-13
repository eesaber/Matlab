%%
clc; clear;
plotSetting = @ Plotsetting_dummy;
po = '/media/akb/2026EF9426EF696C/raw_data/20090811/ALPSRP188921500-L1.1/T3';
pi = '/home/akb/Code/Matlab/PolSAR/output';
x = PolSAR_AnalyzeTool([624 4608],'ALOS PALSAR',plotSetting,...
    'inputDataDir', po,  'outputDataDir', pi);
%% Radar-cross section
x.RCS('vv',true,'Span',true,'hh',true,'hv',true)
NESZ = -29; % (dB)
x.isSNR(NESZ,'vv',true,'hh',true,'hv',true) % Noise equivalent sigam zero (NESZ)

%% Conventional decomposition
x.pauliDecomposition('Filename','Pauli_decomp') % Pauli decomposition
x.fourComponentDecomposition() % 4-component decomposition
[H, alpha] = x.eigenDecomposition(1,'eigen_9x9avg','SaveResults', false); % Eigen-decomposition
x.HAlphaDiagram(H, alpha) % H/alpha diagram

%% Surface parameters
x.paraRatioVVHH(3, 3) % σ_vv/ σ_hh
x.paraGamma12(3, 3) % gamma_12

%% Log-cumulant
[kai_1, kai_2, kai_3] = x.logCumulant(); % Log cumulant

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