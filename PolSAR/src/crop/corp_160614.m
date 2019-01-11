mission_num = 13;
plotSetting = @Plotsetting_dummy;
fout = '/home/akb/Code/Matlab/PolSAR/output/California/160614';
fin = '/media/akb/2026EF9426EF696C/raw_data/California_Valley/160614/';
x = PolSAR_AnalyzeTool([3300 26145],...
                        'UAVSAR', plotSetting,...
                        'inputDataDir', fin,...
                        'outputDataDir', fout);
x.POW_RANGE = [-25 15];
%%
x.RCS('span',true,'vv',true,'hv',true,'hh',true);
x.logCumulant();
x.eigenDecomposition(true,'1');
x.fourComponentDecomposition();
