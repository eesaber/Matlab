%% ERS-2 sea ice 
clc; clear;
plotSetting = @ Plotsetting_dummy;

%{
figure
cbar = colorbar('northoutside');
caxis([-25 -5])
colormap gray
title(cbar,'(dB)','Position', [1150 0])
plot_para('Maximize',true,'Filename', 'colorbar')
%}
%% Winter Radar-cross section
f = [repmat('/media/akb/2026EF9426EF696C/raw_data/20070426/ERS2/',5,1) ...
    ['SAR_IMS_1PNESA20070427_061604_00000018A125_00349_62828_0000.tif';...
    'SAR_IMS_1PNESA20070427_061620_00000018A125_00349_62828_0000.tif';...
    'SAR_IMS_1PNESA20070427_061635_00000018A125_00349_62828_0000.tif';...
    'SAR_IMS_1PNESA20070427_061650_00000018A125_00349_62828_0000.tif';...
    'SAR_IMS_1PNESA20070427_061707_00000016A125_00349_62828_0000.tif']];

%{
for n = 1 : 3
    y = PolSAR_AnalyzeTool([],'ERS-2', plotSetting,...
    'inputDataDir', f(n,:),  'outputDataDir', fout); % [#col, #row]
    figure
    imshow(10*log10(y.vv_vv),[],'Border','tight')
    y.plotSetting(y.POW_RANGE,'Colorbar_unit',"(dB)")
    colormap gray
end
%}
cd /home/akb/Code/Matlab/PolSAR/output/20070426/ers
% image No.2 
y_A = PolSAR_AnalyzeTool([],'ERS-2', plotSetting,...
    'inputDataDir', f(2,:),  'outputDataDir',...
    '/home/akb/Code/Matlab/PolSAR/output/20070426/ers'); % [#col, #row]
y_A.vv_vv = flipud(rot90(y_A.vv_vv));
% image No.3
y_3 = PolSAR_AnalyzeTool([],'ERS-2', plotSetting,...
    'inputDataDir', f(3,:),  'outputDataDir',...
    '/home/akb/Code/Matlab/PolSAR/output/20070426/ers'); % [#col, #row]
y_3.vv_vv = flipud(rot90(y_3.vv_vv));
% image No.4
y_B = PolSAR_AnalyzeTool([],'ERS-2', plotSetting,...
    'inputDataDir', f(4,:),  'outputDataDir',...
    '/home/akb/Code/Matlab/PolSAR/output/20070426/ers'); % [#col, #row]
y_B.vv_vv = flipud(rot90(y_B.vv_vv));
% image No.5
y_C = PolSAR_AnalyzeTool([],'ERS-2', plotSetting,...
    'inputDataDir', f(5,:),  'outputDataDir',...
    '/home/akb/Code/Matlab/PolSAR/output/20070426/ers'); % [#col, #row]
y_C.vv_vv = flipud(rot90(y_C.vv_vv));
%% Region C
%{
figure
imagesc(10*log10(y_C.vv_vv(1:610,1:1550)))
colormap gray
set(gca,'Ydir','normal','clim', [-25 -5])

im_c = [y_B.vv_vv(1:610,2750:2849) y_C.vv_vv(1:610,1:1515)];
im_c = imresize(im_c, [624 4608],'method','nearest','Antialiasing',true);
im_c = imrotate(im_c,1,'nearest','crop');

im_c = [im_c(41:end,:); zeros(40, 4608)];
%}
height = 675;
im_c = [y_B.vv_vv(1:height,2730:2849) y_C.vv_vv(1:height,1:1600)];
im_c = imrotate(im_c,3,'nearest','crop');
im_c = im_c(25:end,1:1654);
im_c = imresize(im_c, x.IMAGE_SIZE,'method','nearest','Antialiasing',true);


figure
imagesc(10*log10(im_c))
colormap gray
set(gca,'Ydir','normal','clim', [-25 -5])
plot_para('Maximize',true,'Filename','sigma_vv_070426_3E')

figure
imagesc(10*log10(im_c)>-12)
colormap gray
set(gca,'Ydir','normal')
plot_para('Maximize',true,'Filename','y_label_070426_3')
%%
old = 10*log10(im_c)>-12;
mask = zeros(size(im_c));
mask(2:end-1,3:end-2) = 1;
mask = mask.*(~isinf(10*log10(im_c)));
nPixel = sum(sum(mask));
% Without HMRF
Tt = sum(sum( ((labels==3) + (labels_MRF==2) ).*(old).*mask))/nPixel;
Ff = sum(sum( (labels==1).*(~old).*mask))/nPixel;
Tf = sum(sum( ((labels==3) + (labels_MRF==2) ).*(~old).*mask))/nPixel;
tF = sum(sum( (labels==1).*(old).*mask))/nPixel;
fprintf('Before MRF:\n%f, %f\n%f ,%f \n', Tt, Tf, tF, Ff);
% With HMRF
Tt = sum(sum( ((labels_MRF==3) ).*(old).*mask))/nPixel;
Ff = sum(sum( (labels_MRF==1).*(~old).*mask))/nPixel;
Tf = sum(sum( ((labels_MRF==3) ).*(~old).*mask))/nPixel;
tF = sum(sum( (labels_MRF==1).*(old).*mask))/nPixel;
fprintf('After MRF:\n%f, %f\n%f ,%f \n', Tt, Tf, tF, Ff);
%%
figure
imagesc((labels==1).*(old).*mask)
figure
imagesc((labels==3).*(~old).*mask)
%% Region B
im_b = [y_3.vv_vv(1:600, 2430:2536), y_B.vv_vv(1:600,1:1600)];
im_b = imrotate(im_b,3.5,'nearest','crop');
im_b = [zeros(100, size(im_b,2)); im_b];
im_b = im_b(1:610,35:1650);
im_b = imresize(im_b, x.IMAGE_SIZE,'method','nearest','Antialiasing',true);

figure
imagesc(10*log10(im_b))
colormap gray
set(gca,'Ydir','normal','clim', [-25 -5])
plot_para('Maximize',true,'Filename','sigma_vv_070426_2E')

figure(3)
imagesc(10*log10(im_b)>-12)
colormap gray
set(gca,'Ydir','normal')
plot_para('Maximize',true,'Filename','y_label_070426_2')
%%
old = 10*log10(im_b)>-12;
mask = zeros(size(im_b));
mask(2:end-1,3:end-2) = 1;
mask = mask.*(~isinf(10*log10(im_b)));
nPixel = sum(sum(mask));
% Without HMRF
Tt = sum(sum( ((labels==3) + (labels==2)).*(old).*mask))/nPixel;
Ff = sum(sum( (labels==1).*(~old).*mask))/nPixel;
Tf = sum(sum( ((labels==3) + (labels==2)).*(~old).*mask))/nPixel;
tF = sum(sum( (labels==1).*(old).*mask))/nPixel;
fprintf('Before MRF:\n%f, %f\n%f ,%f \n', Tt, Tf, tF, Ff);
% With HMRF
Tt = sum(sum( ((labels_MRF==3) + (labels_MRF==2)).*(old).*mask))/nPixel;
Ff = sum(sum( (labels_MRF==1).*(~old).*mask))/nPixel;
Tf = sum(sum( ((labels_MRF==3) + (labels_MRF==2)).*(~old).*mask))/nPixel;
tF = sum(sum( (labels_MRF==1).*(old).*mask))/nPixel;
fprintf('After MRF:\n%f, %f\n%f ,%f \n', Tt, Tf, tF, Ff);
%% Region A
y_A.vv_vv = y_A.vv_vv(780:2329,150:765);
y_A.RCS('vv',true)
movefile sigma_vv.jpg sigma_vv_070426_1E.jpg
figure 
imagesc(10*log10(y_A.vv_vv)>-12)
colormap gray
set(gca,'Ydir','normal')
plot_para('Maximize',true,'Filename','y_label_070426_1')
%% Early melt
