%% ERS-2 sea ice 
clc; clear;
plotSetting = @ Plotsetting_dummy;
fin = '/media/akb/2026EF9426EF696C/raw_data/20070426/SAR_IMS_1PNESA20070427_061635_00000018A125_00349_62828_0000';
fout = '/home/akb/Code/Matlab/PolSAR/output/20070426/';
x = PolSAR_AnalyzeTool([],'ERS-2', plotSetting,...
    'inputDataDir', fin,  'outputDataDir', fout); % [#col, #row]
%% Radar-cross section
x.RCS('vv',true)
