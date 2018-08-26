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
f = [repmat('/media/akb/2026EF9426EF696C/raw_data/20070426/ERS2/',3,1) ...
    [%'SAR_IMS_1PNESA20070427_061604_00000018A125_00349_62828_0000.tif';...
    'SAR_IMS_1PNESA20070427_061620_00000018A125_00349_62828_0000.tif';...
    %'SAR_IMS_1PNESA20070427_061635_00000018A125_00349_62828_0000.tif';...
    'SAR_IMS_1PNESA20070427_061650_00000018A125_00349_62828_0000.tif';...
    'SAR_IMS_1PNESA20070427_061707_00000016A125_00349_62828_0000.tif']];
%%
for n = 1 : 3
    y = PolSAR_AnalyzeTool([],'ERS-2', plotSetting,...
    'inputDataDir', f(n,:),  'outputDataDir', fout); % [#col, #row]
    figure
    imshow(10*log10(y.vv_vv),[],'Border','tight')
    y.plotSetting(y.POW_RANGE,'Colorbar_unit',"(dB)")
    colormap gray
end
%% 480
y_480 = PolSAR_AnalyzeTool([],'ERS-2', plotSetting,...
    'inputDataDir', f(3,:),  'outputDataDir',...
    '/home/akb/Code/Matlab/PolSAR/output/20070426/ers'); % [#col, #row]
y_480.vv_vv = flipud(rot90(y_480.vv_vv(1:1564,1:610)));
y_480.RCS('vv',true)
%%
figure 
imagesc(10*log10(y_480.vv_vv)>-12)
colormap gray
set(gca,'Ydir','normal')
plot_para('Maximize',true,'Filename','y_label_3')
%% 460
y_460 = PolSAR_AnalyzeTool([],'ERS-2', plotSetting,...
    'inputDataDir', f(2,:),  'outputDataDir',...
    '/home/akb/Code/Matlab/PolSAR/output/20070426/ers'); % [#col, #row]
y_460.vv_vv = flipud(rot90(y_460.vv_vv(1:1564,1:610)));
y_460.RCS('vv',true)
figure 
imagesc(10*log10(y_460.vv_vv)>-12)
colormap gray
set(gca,'Ydir','normal')
plot_para('Maximize',true,'Filename','y_label_2')

%% 430
y_430 = PolSAR_AnalyzeTool([],'ERS-2', plotSetting,...
    'inputDataDir', f(1,:),  'outputDataDir',...
    '/home/akb/Code/Matlab/PolSAR/output/20070426/ers'); % [#col, #row]
y_430.vv_vv = flipud(rot90(y_430.vv_vv(779:2420,150:765)));
y_430.RCS('vv',true)
figure 
imagesc(10*log10(y_430.vv_vv)>-12)
colormap gray
set(gca,'Ydir','normal')
plot_para('Maximize',true,'Filename','y_label_1')
%% Early melt
