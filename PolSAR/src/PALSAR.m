
clc; clear;
plotSetting = @ Plotsetting_seaice;
%%
winter = '/media/akb/2026EF9426EF696C/raw_data/20070426/ALPSRP066691460-L1.1/T3'; % 420~540
advance_melt = '/media/akb/2026EF9426EF696C/raw_data/20090811/ALPSRP188921510-L1.1/T3'; % 500~520
fin = advance_melt ;
fout = '/home/akb/Code/Matlab/PolSAR/output/20090811/510';
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
[kai_1, kai_2, kai_3] = x.logCumulant(); % Log cumulant0.00001
x.logCumulantDiagram(kai_2, kai_3)

%% K-means Classfication
nColors = 3;
if 1
    im = cat(3, reshape(intensityMapping(kai_1,'Bit','uint16'), size(kai_1)),...
        reshape(intensityMapping(kai_2,'Bit','uint16'), size(kai_2)), ...
        reshape(intensityMapping(kai_3,'Bit','uint16'), size(kai_3)));
    im = im(2:end-1, 2:end-1);
else
    R = kai_1(2:end-1, 2:end-1);
    G = kai_2(2:end-1, 2:end-1);
    B = kai_3(2:end-1, 2:end-1);
    im = cat(3, reshape(intensityMapping(R,'Bit','uint16'), size(R)),...
        reshape(intensityMapping(G,'Bit','uint16'), size(G)), ...
        reshape(intensityMapping(B,'Bit','uint16'), size(B)));
    clear R G B
end

y_hat = imageKCLuster(x, im, nColors);
%%
labels = padarray(y_hat, [1,1], -1,'both');
showLabels(x, labels, nColors)

%% Make old ice to 3 and new ice to 1
temp = labels;
labels(temp == 1) = 3;
labels(temp == 2) = 1;
labels(temp == 3) = 2;

clear temp
%%
showLabels(x, labels, nColors-1)

%% Hidden Markove random fields 
temp = labels(2:end-1, 2:end-1);
y_hatMRF = uint8(HMRF(x, im, temp , nColors-1));
clear temp;
%%
labels_MRF = padarray(y_hatMRF, [1,1], -1,'both');
showLabels(x, labels_MRF, nColors)
%%
temp = labels_MRF;
labels_MRF(temp == 3) = 1;
labels_MRF(temp == 1) = 3;
%labels_MRF(temp == 1) = 2;
%%
temp = labels_MRF;
temp(labels_MRF==2) = 3;
showLabels(x, temp, nColors)
clear temp

    
%% 2D 
figure
hold on 
for n = 1 : nColors
    scatter(kai_3(labels==n), kai_2(labels==n),'filled')
end
hold off
set(gca, 'Xtick', -5:5,'xlim',[-3 3], 'ylim', [0 2])
xlabel('$\kappa_3$','interpreter','latex')
ylabel('$\kappa_2$','interpreter','latex')
%% 3D
figure
hold on 
for n = 1 : nColors
    scatter3(kai_3(labels==n), kai_2(labels==n), kai_1(labels==n), 'filled')
end
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