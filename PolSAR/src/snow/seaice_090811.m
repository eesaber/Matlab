clc; clear;
%%
img4roi = H;
bw_testing = logical(zeros([x_b.IMAGE_SIZE, 10]));
%%
figure
[bw_testing(:,:,8), ver_x, ver_y] = roipoly(flipud(img4roi));
%%
bw = zeros(x_b.IMAGE_SIZE);
for it = 1 : size(bw_testing,3), bw = bw+bw_testing(:,:,it); end
bw = bw > 0;
figure
imagesc(bw)

%% Reference Data
%save('/media/akb/2026EF9426EF696C/raw_data/20090811/groundtruth', 'bw', 'bw_testing')
%gt = logical(bw);
%save('/home/akb/Code/PolSAR_ML/data/final_version/mask_090811', 'gt')
load('/media/akb/2026EF9426EF696C/raw_data/20090811/groundtruth')

gt = bw;
figure
imagesc(gt)
set(gca,'Ydir','normal','Visible','off')
colormap(gca, [0,120,0; 255,255,0]/255)
plot_para('Filename',[x_b.OUTPUT_PATH '/label_090811'],'Maximize',true)
%save('/home/akb/Code/PolSAR_ML/data/final_version/image_090811_(2).mat','im')
%save('/home/akb/Code/PolSAR_ML/data/final_version/mask_090811.mat', 'gt')
%% Data 
plotSetting = @ Plotsetting_seaice;
fin = [repmat('/media/akb/2026EF9426EF696C/raw_data/20090811/',3,1) ...
         ['ALPSRP188921500-L1.1/T3';...
          'ALPSRP188921510-L1.1/T3';...
          'ALPSRP188921520-L1.1/T3']];
fout = [repmat('/home/akb/Code/Matlab/PolSAR/output/20090811/',3,1) ...
        ['500'; '510';'520']];
% data shape is [#col, #row]
x_b = SeaIce([624 4608],'ALOS PALSAR', plotSetting,...
    'inputDataDir', fin(2,:),  'outputDataDir', fout(2,:));
x_c = SeaIce([624 4608],'ALOS PALSAR', plotSetting,...
    'inputDataDir', fin(3,:),  'outputDataDir', fout(3,:));
x_b.vv_vv = [x_b.vv_vv, x_c.vv_vv(:,678:end)];
x_b.hh_hh = [x_b.hh_hh, x_c.hh_hh(:,678:end)];
x_b.hv_hv = [x_b.hv_hv, x_c.hv_hv(:,678:end)];
x_b.hh_vv = [x_b.hh_vv, x_c.hh_vv(:,678:end)];
x_b.hh_hv = [x_b.hh_hv, x_c.hh_hv(:,678:end)];
x_b.hv_vv = [x_b.hv_vv, x_c.hv_vv(:,678:end)];
x_b.IMAGE_SIZE = size(x_b.vv_vv);
clear x_c
%% Generate data
input_vector = 2;
[im, texture] = x_b.generateImage4Classification(input_vector,'texture',texture);
if 0
    bai = 100;
    pixel_cov = bai*permute( ...
    reshape( ...
    reshape( ...
    cat(3, x_b.hh_hh, sqrt(2)*x_b.hh_hv, x_b.hh_vv, ...
    sqrt(2)*conj(x_b.hh_hv), 2*x_b.hv_hv, sqrt(2)*x_b.hv_vv, ...
    conj(x_b.hh_vv), sqrt(2)*conj(x_b.hv_vv), x_b.vv_vv), ...
    [], 9).', ...
    3,3,[]), ...
    [2, 1, 3]);
end
%save(['/home/akb/Code/PolSAR_ML/data/image_090811_('+ num2str(input_vector)+').mat'],'im')
%% Test sigle algorithm
algo = 'GFCM';
num_c = 5;
num_subc = 2;
distance_metric = 'squaredeuclidean';
%distance_metric = 'wishart';
drg_m = 2;
drg_n = 2;
tic 
[label_fcm, c, u] = x_b.myFCM(double(im),...
                            'clusterNum', int32(num_c), ...
                            'subClusterNum', int32(num_subc), ...
                            'max_iter', int32(50), ...
                            'drg_m', drg_m, ...
                            'drg_n', drg_n, ...
                            'algo', algo, ...
                            'distance', distance_metric, ...
                            'covariance', pixel_cov);
toc
%[c, u] = fcm(im,num_c,[drg_m, 75, 1e-5, false]);
%[~, label_fcm] = max(u,[],1);
%%
x_b.showLabels(reshape(label_fcm,x_b.IMAGE_SIZE), num_c)
title(sprintf('(squaredeuclidean) algorithm: %s, (I,J,m,n) = (%i,%i,%i,%i)', algo, num_c, num_subc, drg_m, drg_n))
%title([algo,' $(I,J,m,n) = ($', num2str(num_c),',',num2str(num_subc),',',...
%    num2str(drg_m),',',num2str(drg_n),')'],'Interpreter', 'latex')
set(gca,'Ydir','normal','yticklabels',[],'xticklabels',[])

%%
label_svm = x_b.mySVM(double(reshape(im, [x_b.IMAGE_SIZE, size(im,2)])), double(bw));
%%
k_cluster = 5;
[label_kmeans, ~] = kmeans(im, k_cluster,'distance','sqeuclidean', ...
                                          'Replicates',3,'MaxIter',150,...
                                          'Options',statset('Display','final'));
%%
x_b.showLabels(reshape(label_kmeans,x_b.IMAGE_SIZE), k_cluster)
title(sprintf('k-means cluster number: %i',k_cluster))
set(gca,'yticklabels',[],'xticklabels',[])
%% NN 
load('/home/akb/Code/PolSAR_ML/data/output/y_hat_090811.mat')
label_nn = y_hat;
clear y_hat
%%
y_hat = reshape(label_fcm, x_b.IMAGE_SIZE);
y_hat = y_hat == 2;
figure
imagesc(y_hat)
set(gca,'Ydir','normal','Visible','off')
colormap(gca, [0,120,0; 255,255,0]/255)
plot_para('Filename',[x_b.OUTPUT_PATH '/label_090811_svm'],'Maximize',true)
%%
disp('--------------fcm-------------')
fprintf('Acc: %f\n', sum(sum(y_hat == bw))/numel(bw))
% Uppercase: predict, lowercase: actual
Mm = sum(sum((y_hat==1).*(bw==1)))/numel(bw);
Mf = sum(sum((y_hat==1).*(bw==0)))/numel(bw);
Fm = sum(sum((y_hat==0).*(bw==1)))/numel(bw);
Ff = sum(sum((y_hat==0).*(bw==0)))/numel(bw);
fprintf('Acc:\n %f, %f\n %f ,%f \n', Mm, Mf, Fm, Ff);

%%
%%
%x_b.HAlphaDiagram(H, alpha) % H/alpha diagram
%x_b.HAlphaDiagram(H(:,4270:8539), alpha(:,4270:8539)) % H/alpha diagram
x_b.HAlphaDiagram(H(:,1:4269), alpha(:,1:4269)) % H/alpha diagram
%%
temp = 10*log10(x_b.vv_vv./x_b.hh_hh);
figure
    histogram2(H(:), temp(:),'DisplayStyle','tile',...
        'YBinLimits',[-4 4],'XBinLimits',[0 1],'NumBins',[750,750],...
        'Normalization','pdf');
    xlabel('$H$','Interpreter', 'latex')
    ylabel('$\sigma_{vv}/\sigma_{hh}$','Interpreter', 'latex')
    colormap jet
    c = colorbar;
    set(gca,'clim',[0,0.8],'Color',[0 0 0]+0.7,'XGrid','off','YGrid','off')
    plot_para('Maximize',true, 'Filename',[x_b.OUTPUT_PATH '/H_Copol_090811'])
    %c.Ticks = [0.1:0.1:0.3,0.35];
%%
figure
imagesc(H)
x_b.plotSetting([0 1])
set(gca,'Visible','off')
colorbar
plot_para('Maximize',true,'Filename',[x_b.OUTPUT_PATH '/Entropy'])
%%
x_b.RCS('vv',true)
%imagesc(10*log10([x_b.vv_vv, x_c.vv_vv(:,678:end)]))
hold on 
plot(x_b.IMAGE_SIZE(2)*[1,1]/2  ,[-100,800],'b','linewidth',2)
set(gca,'Ydir','normal','yticklabels',[],'xticklabels',[],'clim',[-25 -5],'xColor','none','yColor','none','box','off')
colormap gray
plot_para('Filename','/home/akb/Code/Matlab/PolSAR/output/090811','Maximize',true)
