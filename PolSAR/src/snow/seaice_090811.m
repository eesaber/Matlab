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
%save('/media/akb/2026EF9426EF696C/raw_data/20090811/groundtruth', 'bw', 'bw_testing')
%gt = logical(bw);
%save('/home/akb/Code/PolSAR_ML/data/final_version/mask_090811', 'gt')

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
global summer 
summer = 1;
%% Reference Data
load('/media/akb/2026EF9426EF696C/raw_data/20090811/groundtruth')
x_b.y_hat = flipud(bw~=0);
clear gt bw_testing bw
figure
imagesc(x_b.y_hat)
set(gca,'Ydir','normal','Visible','off')
colormap(gca, [0,120,0; 255,255,0]/255)
plot_para('Filename',[x_b.OUTPUT_PATH '/label_090811'],'Maximize',true)

%% Generate data
input_vector = 4;
if input_vector==4 || input_vector==3
    [~] = x_b.generateImage4Classification(input_vector);
    if input_vector == 4
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
else
    x_b.generateImage4Classification(input_vector,'texture',texture);
    %[~, texture] = x_b.generateImage4Classification(input_vector);
end
%save(['/home/akb/Code/PolSAR_ML/data/image_090811_('+ num2str(input_vector)+').mat'],'im')
%% Test sigle algorithm
algo = 'GFCM';
num_c = 3;
num_subc = 2;
%distance_metric = 'squaredeuclidean';
distance_metric = 'wishart';
drg_m = 2;
drg_n = 2;
[label_fcm, c, u] = x_b.myFCM(x_b.im,...
                            'clusterNum', int32(num_c), ...
                            'subClusterNum', int32(num_subc), ...
                            'max_iter', int32(15), ...
                            'u', 1,...
                            'drg_m', drg_m, ...
                            'drg_n', drg_n, ...
                            'algo', algo, ...
                            'distance', distance_metric, ...
                            'covariance', pixel_cov);

%%
x_b.showLabels(reshape(label_fcm,x_b.IMAGE_SIZE), num_c)
z = sprintf('%s %s, $(I,J,m,n)$ = (%i,%i,%i,%i)', algo, distance_metric, num_c, num_subc, drg_m, drg_n);
title(z,'Interpreter', 'latex')
set(gca,'Ydir','normal','Visible','off')
plot_para('ratio',[16 9 1], 'Filename', [x_b.OUTPUT_PATH '/label_fcm_4_3c'],'maximize',true)
%%
label_svm = x_b.mySVM(double(reshape(x_b.im, [x_b.IMAGE_SIZE, size(x_b.im,2)])), double(x_b.y_hat));
%%
global kpp
kpp = 1;
%%
k_cluster = 3;
distance_metric = 'wishart'; 
if strcmp(distance_metric, 'wishart')
    [label_kmeans, c] = x_b.myKmeans(x_b.im,...
                            'maxClusterNum', int32(k_cluster),...
                            'max_iter', int32(30), ... % test if 10 times will cause warning
                            'replicates',int32(1), ...
                            'distance', distance_metric, ...
                            'covariance', pixel_cov);
    u_pre = 1/(k_cluster-1)*single(ones(size(x_b.im,1),k_cluster));
    for it = 1 : k_cluster
        u_pre(label_kmeans==it,it) = 1;
    end
else
    [label_kmeans, ~] = kmeans(x_b.im, k_cluster,...
                        'distance','sqeuclidean',...
                        'Replicates',3,...
                        'MaxIter',30,...
                        'Options',statset('Display','final'));
end                    
%%
x_b.showLabels(reshape(label_kmeans,x_b.IMAGE_SIZE), k_cluster)
title(sprintf('k-means cluster number: %i',k_cluster))
set(gca,'yticklabels',[],'xticklabels',[])
%% NN 
load('/home/akb/Code/PolSAR_ML/data/output/y_hat_090811.mat')
label_nn = y_hat;
clear y_hat
%%
x_b.y_hat = reshape(label_kmeans, x_b.IMAGE_SIZE);
x_b.y_hat = x_b.y_hat == 2;
figure
imagesc(x_b.y_hat)
set(gca,'Ydir','normal','Visible','off')
colormap(gca, [0,120,0; 255,255,0]/255)
plot_para('Filename',[x_b.OUTPUT_PATH '/label_090811_km_4'],'Maximize',true)

disp('--------------fcm-------------')
fprintf('Acc: %f\n', sum(sum(x_b.y_hat == x_b.y_hat))/numel(x_b.y_hat))
% Uppercase: predict, lowercase: actual
Mm = sum(sum((x_b.y_hat==1).*(x_b.y_hat==1)))/numel(x_b.y_hat);
Mf = sum(sum((x_b.y_hat==1).*(x_b.y_hat==0)))/numel(x_b.y_hat);
Fm = sum(sum((x_b.y_hat==0).*(x_b.y_hat==1)))/numel(x_b.y_hat);
Ff = sum(sum((x_b.y_hat==0).*(x_b.y_hat==0)))/numel(x_b.y_hat);
fprintf('Acc:\n %f, %f\n %f ,%f \n', Mm, Mf, Fm, Ff);
clear Mm Mf Fm Ff
%% CNN
img_h = 256;
img_w = 256;
new_size = [512 8448];
%%
it = 1;
temp_im = imresize(reshape(x_b.im, [x_b.IMAGE_SIZE, size(x_b.im,2)]), new_size, 'method', 'nearest');
x_val = single(zeros(new_size(1)/img_h*new_size(2)/img_w, img_h, img_w, size(x_b.im,2)));
disp(size(x_val))
for it = 1 : size(x_b.im,2)
    tile = im2col(temp_im(:,:,it),[img_h, img_w],'distinct');
    tile = permute(reshape(tile, img_h, img_w, []),[3,1,2]);
    for it_t = 1 : size(tile,1)
        x_val(it_t,:,:,it) = squeeze(tile(it_t,:,:));

    end
    it = it +1;
end
clear temp_im 
y_val = imresize(reshape(x_b.y_hat, x_b.IMAGE_SIZE), new_size, 'method', 'nearest');
y_val = reshape(im2col(y_val,[img_h, img_w],'distinct'), img_h, img_w, []);
y_val = permute(y_val, [3,1,2]);
save(['/home/akb/Code/PolSAR_ML/data_val/x_val_090811_', num2str(input_vector),'.mat'], 'x_val')
save('/home/akb/Code/PolSAR_ML/data_val/y_val_090811.mat', 'y_val')
%%
load(['/home/akb/Code/PolSAR_ML/output/y_hat_090811_', num2str(input_vector),'.mat'])
bw_cnn = imresize(reshape(bw, x_b.IMAGE_SIZE), new_size, 'method', 'nearest');
temp = y_hat;
x_b.y_hat = zeros(size(bw_cnn));
if 1
    it = 1;
    for it_c = 1 : img_w : size(bw_cnn, 2)
        for it_r = 1 : img_h: size(bw_cnn,1)
            x_b.y_hat(it_r: it_r+img_h-1, it_c: it_c+img_w-1) = temp(it,:,:,1);
            it = it+1;
        end
    end
    x_b.y_hat = x_b.y_hat<0.5;
else 
    x_b.y_hat = squeeze(x_b.y_hat(:,:,:,1))<0.5;
    bw_cnn = imresize(x_b.y_hat, [496, 496], 'method','nearest');
end
figure
imagesc(x_b.y_hat)
set(gca,'Ydir','normal','Visible','off')
colormap(gca, [0,120,0; 255,255,0]/255)
plot_para('Filename', [x_b.OUTPUT_PATH '/label_090811_cnn'],'Maximize',true)

disp('--------------cnn-------------')
fprintf('Acc: %f\n', sum(sum(x_b.y_hat == bw_cnn))/sum(sum(bw_cnn>=0)))
% Uppercase: predict, lowercase: actual
Mm = sum(sum((x_b.y_hat==1).*(bw_cnn==1)))/sum(sum(bw_cnn>=0));
Mf = sum(sum((x_b.y_hat==1).*(bw_cnn==0)))/sum(sum(bw_cnn>=0));
Fm = sum(sum((x_b.y_hat==0).*(bw_cnn==1)))/sum(sum(bw_cnn>=0));
Ff = sum(sum((x_b.y_hat==0).*(bw_cnn==0)))/sum(sum(bw_cnn>=0));
fprintf('Acc:\n %f, %f\n %f ,%f \n', Mm, Mf, Fm, Ff);
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
%imagesc(10*log10([x_b.vv_vv, x_b.vv_vv(:,678:end)]))
hold on 
plot(x_b.IMAGE_SIZE(2)*[1,1]/2  ,[-100,800],'b','linewidth',2)
set(gca,'Ydir','normal','yticklabels',[],'xticklabels',[],'clim',[-25 -5],'xColor','none','yColor','none','box','off')
colormap gray
plot_para('Filename','/home/akb/Code/Matlab/PolSAR/output/090811','Maximize',true)
%% Geographical size of image 
%distance(lat1,lon1,lat2,lon2)
disp('ROI 3')
fprintf('Ground range: %f (km)\n',distance('gc',74.30728912353516,-117.95664855815409,74.40733805010396,-117.06706010794709, wgs84Ellipsoid('km')))
fprintf('Along track length: %f (km)\n', distance('gc',74.30728912353516,-117.95664855815409,75.29965659942047,-119.63900528494702,referenceEllipsoid('GRS80','km')))