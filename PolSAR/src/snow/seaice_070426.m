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
if 0
    disp('Appling k-means clustering to x_a')
    x_a.imageKCluster(nColors)
    x_a.showLabels(x_a.y_hat, nColors)
end
if 0
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

temp = x_c.y_hat;
x_c.y_hat(temp == 1) = 2;
x_c.y_hat(temp == 2) = 1;

clear temp
x_c.showLabels(x_c.y_hat, nColors)

%% Hidden Markove random fields 
if 0
    disp('Appling HMRF to x_a')
    x_a.HMRF(nColors);
    x_a.showLabels(x_a.y_hat_MRF, nColors)
end
if 1
    disp('Appling HMRF to x_b')
    x_b.HMRF(nColors);
    x_b.showLabels(x_b.y_hat_MRF, nColors-1)
end
if 1
    disp('Appling HMRF to x_c')
    x_c.HMRF(nColors);
    x_c.showLabels(x_c.y_hat_MRF, nColors-1)
end
%% Image B

old = 10*log10(y_b.vv_vv)>-12;
mask = zeros(size(y_b.vv_vv));
mask(2:end-1,2:end-1) = 1;
mask = mask.*(~isinf(10*log10(y_b.vv_vv)));
nPixel = sum(sum(mask));
% Without HMRF
Tt = sum(sum( ((x_b.y_hat==3) + (x_b.y_hat==2)).*(old).*mask))/nPixel;
Ff = sum(sum( (x_b.y_hat==1).*(~old).*mask))/nPixel;
Tf = sum(sum( ((x_b.y_hat==3) + (x_b.y_hat==2)).*(~old).*mask))/nPixel;
tF = sum(sum( (x_b.y_hat==1).*(old).*mask))/nPixel;
fprintf('Before MRF:\n%f, %f\n%f ,%f \n', Tt, Tf, tF, Ff);
% With HMRF
Tt = sum(sum( ((x_b.y_hat_MRF==3) + (x_b.y_hat_MRF==2)).*(old).*mask))/nPixel;
Ff = sum(sum( (x_b.y_hat_MRF==1).*(~old).*mask))/nPixel;
Tf = sum(sum( ((x_b.y_hat_MRF==3) + (x_b.y_hat_MRF==2)).*(~old).*mask))/nPixel;
tF = sum(sum( (x_b.y_hat_MRF==1).*(old).*mask))/nPixel;
fprintf('After MRF:\n%f, %f\n%f ,%f \n', Tt, Tf, tF, Ff);

%% Image C
% Possible frost flower
c = [3613,4376,4492,4573,4608,4608,4527,4117];
r = [1,95,131,124,133,96,66,1];
ROI_frost = ~roipoly(ones(x_c.IMAGE_SIZE),c,r);

old = 10*log10(y_c.vv_vv)>-12;
mask = zeros(size(y_c.vv_vv));
mask(2:end-1,2:end-1) = 1;
mask = mask.*(~isinf(10*log10(y_c.vv_vv)));
nPixel = sum(sum(mask));
% Without HMRF
Tt = sum(sum( ((x_c.y_hat==3) + (x_c.y_hat==2)).*(old).*mask))/nPixel;
Ff = sum(sum( (x_c.y_hat==1).*(~old).*mask))/nPixel;
Tf = sum(sum( ((x_c.y_hat==3) + (x_c.y_hat==2)).*(~old).*mask))/nPixel;
tF = sum(sum( (x_c.y_hat==1).*(old).*mask))/nPixel;
fprintf('Before MRF:\n%f, %f\n%f ,%f \n', Tt, Tf, tF, Ff);
% Without HMRF mask forst flower
mask = ROI_frost.*mask;
nPixel = sum(sum(mask));
Tt = sum(sum( ((x_c.y_hat==3) + (x_c.y_hat==2)).*(old).*mask))/nPixel;
Ff = sum(sum( (x_c.y_hat==1).*(~old).*mask))/nPixel;
Tf = sum(sum( ((x_c.y_hat==3) + (x_c.y_hat==2)).*(~old).*mask))/nPixel;
tF = sum(sum( (x_c.y_hat==1).*(old).*mask))/nPixel;
fprintf('Before MRF:\n%f, %f\n%f ,%f \n', Tt, Tf, tF, Ff);

% With HMRF
Tt = sum(sum( ((x_c.y_hat_MRF==3) ).*(old).*mask))/nPixel;
Ff = sum(sum( (x_c.y_hat_MRF==1).*(~old).*mask))/nPixel;
Tf = sum(sum( ((x_c.y_hat_MRF==3) ).*(~old).*mask))/nPixel;
tF = sum(sum( (x_c.y_hat_MRF==1).*(old).*mask))/nPixel;
fprintf('After MRF:\n%f, %f\n%f ,%f \n', Tt, Tf, tF, Ff);