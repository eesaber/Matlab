%% ALOS PALSAR sea ice 
clc; clear; close all;
plotSetting = @ Plotsetting_dummy;
fin = '/media/akb/2026EF9426EF696C/raw_data/20070426/ALPSRP066691480-L1.1/T3';
fout = '/home/akb/Code/Matlab/PolSAR/output/20070426/480';
x = PolSAR_AnalyzeTool([624 4608],'ALOS PALSAR', plotSetting,...
    'inputDataDir', fin,  'outputDataDir', fout); % [#col, #row]

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
x.paraRatioVVHH(2, 4) % σ_vv/ σ_hh
x.paraRoughness1(2, 4) % gamma_12
x.paraMoisture(2, 4) % 
x.getPolAngle([-30 30]); % polarimetric angle
%% Log-cumulant
[kai_1, kai_2, kai_3] = x.logCumulant(); % Log cumulant
%%
kai_1 = kai_1(2:end-1,3:end-2);
kai_2 = kai_2(2:end-1,3:end-2);
kai_3 = kai_3(2:end-1,3:end-2);
%%
figure
histogram2(kai_3(:), kai_2(:), 750,...
        'DisplayStyle','tile','YBinLimits',[0 5],'XBinLimits',[-3 3])
colormap jet
set(gca,'Color',[0 0 0]+0.7,'XGrid','off','YGrid','off')
hold on 
var = linspace(0,4,1000);
plot(psi(3,var), -psi(2,var), 'w' , -psi(3,var), -psi(2,var), 'w--', 'Linewidth',2.5)
set(gca, 'Xtick', -5:5,'xlim',[-3 3], 'ylim', [0 2])
xlabel('$\kappa_3$','interpreter','latex')
ylabel('$\kappa_2$','interpreter','latex')
plot_para('Maximize',true,'Filename','log_cumulant_diagram','Ratio',[4 3 1])
hold off
%% Segementation
temp = cat(3, reshape(Stat(kai_1,'uint16'), size(kai_1)),...
    reshape(Stat(kai_2, 'uint16'), size(kai_1)), ...
    reshape(Stat(kai_3, 'uint16'), size(kai_1)));
%temp = cat(3, kai_1, kai_2, kai_3);
lab_I = rgb2lab(temp);
clear temp
ab = lab_I(:,:,2:3);
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,size(ab,1)* size(ab,2),2);
nColors = 3;
% repeat the clustering 3 times to avoid local minima
[cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','sqeuclidean', ...
                                      'Replicates',5,'MaxIter',500);
pixel_labels = uint8(reshape(cluster_idx,nrows,ncols));
clear ab
%%
figure
%imshow(pixel_labels,[],'Border','tight')
imagesc(pixel_labels)
cmap  = flag(nColors);
colormap(cmap)
L = line(ones(nColors),ones(nColors), 'LineWidth',2);               % generate line 
set(L,{'color'},mat2cell(cmap, ones(1,nColors),3)); % set the colors according to cmap
legend(num2str((1:nColors)'),'Location','bestoutside')
plot_para('Maximize',true,'Filename','label')
%%
labels = padarray(pixel_labels, [1,2], -1,'both');
kai_1 = padarray(kai_1, [1,2], -1,'both');
kai_2 = padarray(kai_2, [1,2], -1,'both');
kai_3 = padarray(kai_3, [1,2], -1,'both');
%%
bin_limit = [1e-5, 0.15];
figure
    hold on 
    for n = 1:nColors
        histogram(x.vv_vv(labels==n),'BinLimits', bin_limit,'Normalization','probability')
    end
    hold off
    legend('label 1','label 2','label 3','label 4', 'label 5')
    xlabel('$|S_{vv}|^2$', 'interpreter', 'latex')
    ylabel('Pr$(x = |S_{vv}|^2)$', 'interpreter', 'latex')
    plot_para('FileName','histo_vv','Maximize',true)
%% 2D 
figure
hold on 
scatter(kai_3(labels==1), kai_2(labels==1),'filled')
scatter(kai_3(labels==2), kai_2(labels==2),'filled')
scatter(kai_3(labels==3), kai_2(labels==3),'filled')
scatter(kai_3(labels==4), kai_2(labels==4),'filled')
hold off
set(gca, 'Xtick', -5:5,'xlim',[-3 3], 'ylim', [0 2])
xlabel('$\kappa_3$','interpreter','latex')
ylabel('$\kappa_2$','interpreter','latex')
%% 3D
figure
hold on 
scatter3(kai_3(labels==1), kai_2(labels==1), kai_1(labels==1), 'filled')
scatter3(kai_3(labels==2), kai_2(labels==2), kai_1(labels==2), 'filled')
scatter3(kai_3(labels==3), kai_2(labels==3), kai_1(labels==3), 'filled')
scatter3(kai_3(labels==4), kai_2(labels==4), kai_1(labels==4), 'filled')
hold off
%set(gca, 'Xtick', -5:5,'xlim',[-3 3], 'ylim', [0 2])
xlabel('$\kappa_3$','interpreter','latex')
ylabel('$\kappa_2$','interpreter','latex')
ylabel('$\kappa_1$','interpreter','latex')
%% Fitting the distribution to 
% Possible distribution: Gamma dist, Inverse Gamma, Fisher, K distribution. 
r = linspace(0,0.2,1000);
figure
hold on 
for n = 1 : nColors
    pd = fitdist(double(x.vv_vv(labels==n)),'Gamma');
    plot(r,pdf(pd,r),'Linewidth',2.5)
end
hold off
%{
%% 2009/08/11 meteorology
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
%}