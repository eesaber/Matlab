clc; clear;

%% Data 
plotSetting = @ Plotsetting_dummy;
fin = [repmat('/media/akb/2026EF9426EF696C/raw_data/20090811/',3,1) ...
         ['ALPSRP188921500-L1.1/T3';...
          'ALPSRP188921510-L1.1/T3';...
          'ALPSRP188921520-L1.1/T3']];
fout = [repmat('/home/akb/Code/Matlab/PolSAR/output/20090811/',3,1) ...
        ['500'; '510';'520']];
% data shape is [#col, #row]
x_a = SeaIce([624 4608],'ALOS PALSAR', plotSetting,...
    'inputDataDir', fin(1,:),  'outputDataDir', fout(1,:));
x_b = SeaIce([624 4608],'ALOS PALSAR', plotSetting,...
    'inputDataDir', fin(2,:),  'outputDataDir', fout(2,:));
x_c = SeaIce([624 4608],'ALOS PALSAR', plotSetting,...
    'inputDataDir', fin(3,:),  'outputDataDir', fout(3,:));
%% Reference Data
fin = ['/media/akb/2026EF9426EF696C/raw_data/20090811/ers/'...
    'SAR_IMS_1PNESA20090811_194714_00000018A149_00314_74817_0000.tif'];
fout = '/home/akb/Code/Matlab/PolSAR/output/20090811/ers';
y = PolSAR_AnalyzeTool([],'ERS-2', plotSetting,...
    'inputDataDir', fin,  'outputDataDir', fout);
y.vv_vv = flipud(rot90(y.vv_vv));

%%
vq = int16(interp1([1 624],[1167 2958], 1:624,'linear'));
mask = ones(x_c.IMAGE_SIZE);
for n = 1 :624
    mask(n,vq(n):end) = 0;
end
figure
imagesc(10*log10((x_c.vv_vv)).*mask)
colormap gray
set(gca,'Ydir','normal','clim', [-25 -5])
plot_para('Maximize',true)

%%
height = 675;
temp = [y_b.vv_vv(1:height,2730:2849) y_c.vv_vv(1:height,1:1600)];
temp = imrotate(temp,3,'nearest','crop');
temp = temp(25:end,1:1654);
y_c.vv_vv = imresize(temp, x_c.IMAGE_SIZE,'method','nearest','Antialiasing',true);
 

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
if 0
    disp('Appling k-means clustering to x_a')
    x_a.imageKCluster(nColors)
    x_a.showLabels(x_a.y_hat, nColors)
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
if 0
    disp('Appling HMRF to x_a')
    x_a.imageKCluster(nColors)
    x_a.showLabels(x_a.y_hat, nColors)
end
if 1
    disp('Appling HMRF to x_b')
    x_b.HMRF(nColors);
    x_b.showLabels(x_b.y_hat_MRF, nColors)
end
if 1
    disp('Appling HMRF to x_c')
    x_c.HMRF(nColors);
    x_c.showLabels(x_c.y_hat_MRF, nColors)
end

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
