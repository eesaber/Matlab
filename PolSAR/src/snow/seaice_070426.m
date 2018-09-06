clc; clear;

%% Data 
plotSetting = @ Plotsetting_dummy;
fin = [repmat('/media/akb/2026EF9426EF696C/raw_data/20070426/',3,1) ...
         ['ALPSRP066691430-L1.1/T3';...
          'ALPSRP066691460-L1.1/T3';...
          'ALPSRP066691480-L1.1/T3']];
fout = [repmat('/home/akb/Code/Matlab/PolSAR/output/20070426/',3,1) ...
        ['430'; '460';'480']];
% data shape is [#col, #row]
x_a = SeaIce([624 4608],'ALOS PALSAR', plotSetting,...
    'inputDataDir', fin(1,:),  'outputDataDir', fout(1,:));
x_b = SeaIce([624 4608],'ALOS PALSAR', plotSetting,...
    'inputDataDir', fin(2,:),  'outputDataDir', fout(2,:));
x_c = SeaIce([624 4608],'ALOS PALSAR', plotSetting,...
    'inputDataDir', fin(3,:),  'outputDataDir', fout(3,:));
%% Reference Data
fin = [repmat('/media/akb/2026EF9426EF696C/raw_data/20070426/ERS2/',5,1) ...
        ['SAR_IMS_1PNESA20070427_061604_00000018A125_00349_62828_0000.tif';...
        'SAR_IMS_1PNESA20070427_061620_00000018A125_00349_62828_0000.tif';...
        'SAR_IMS_1PNESA20070427_061635_00000018A125_00349_62828_0000.tif';...
        'SAR_IMS_1PNESA20070427_061650_00000018A125_00349_62828_0000.tif';...
        'SAR_IMS_1PNESA20070427_061707_00000016A125_00349_62828_0000.tif']];
fout = '/home/akb/Code/Matlab/PolSAR/output/20070426/ers';
% image No.3
y_3 = PolSAR_AnalyzeTool([],'ERS-2', plotSetting,...
    'inputDataDir', fin(3,:),  'outputDataDir', fout);
y_3.vv_vv = flipud(rot90(y_3.vv_vv));
% image No.4
y_b = PolSAR_AnalyzeTool([],'ERS-2', plotSetting,...
    'inputDataDir', fin(4,:),  'outputDataDir', fout);
y_b.vv_vv = flipud(rot90(y_b.vv_vv));
% image No.5
y_c = PolSAR_AnalyzeTool([],'ERS-2', plotSetting,...
    'inputDataDir', fin(5,:),  'outputDataDir', fout);
y_c.vv_vv = flipud(rot90(y_c.vv_vv));
% Region C
height = 675;
temp = [y_b.vv_vv(1:height,2730:2849) y_c.vv_vv(1:height,1:1600)];
temp = imrotate(temp,3,'nearest','crop');
temp = temp(25:end,1:1654);
y_c.vv_vv = imresize(temp, x_c.IMAGE_SIZE,'method','nearest','Antialiasing',true);
% Region B
temp = [y_3.vv_vv(1:600, 2430:2536), y_b.vv_vv(1:600,1:1600)];
temp = imrotate(temp,3.5,'nearest','crop');
temp = [zeros(100, size(temp,2)); temp];
temp = temp(1:610,35:1650);
y_b.vv_vv = imresize(temp, x_b.IMAGE_SIZE,'method','nearest','Antialiasing',true);

clear fin fout temp y_3
%% Compute log-cumulant
x_a.logCumulant();
x_b.logCumulant();
x_c.logCumulant();
%x_430.logCumulantDiagram(kai_2, kai_3)
%x_460.logCumulantDiagram(kai_2, kai_3)
%x_480.logCumulantDiagram(kai_2, kai_3)
%% K-means Classfication
nColors = 3;
if 1
    disp('Appling k-means clustering to x_a')
    x_a.imageKCluster(nColors)
    x_b.showLabels(x_b.y_hat, nColors)
end
if 1
    disp('Appling k-means clustering to x_b')
    x_b.imageKCluster(nColors)
    x_b.showLabels(x_b.y_hat, nColors)
end
if 1
    disp('Appling k-means clustering to x_c')
    x_c.imageKCluster(nColors)
    x_c.showLabels(x_c.y_hat, nColors)
end

%% Make old ice to 3 and new ice to 1
%{
label = x_c.y_hat;
temp = x_c.y_hat;
labels(temp == 1) = 3;
labels(temp == 2) = 1;
labels(temp == 3) = 2;
clear temp
showLabels(x, labels, nColors-1)
%}
%% Hidden Markove random fields 
x_c.HMRF(nColors);
x_c.showLabels(x_c.y_hat_MRF, nColors)
%% Make old ice to 3 and new ice to 1
%{
label = x_c.y_hat_MRF;
temp = x_c.y_hat_MRF;
labels(temp == 1) = 3;
labels(temp == 2) = 1;
labels(temp == 3) = 2;
clear temp
showLabels(x, labels, nColors-1)
%}

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
