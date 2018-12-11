%% Data IO
plotSetting = @ Plotsetting_dummy;
fin = [repmat('/media/akb/2026EF9426EF696C/raw_data/20070426/',3,1) ...
         ['ALPSRP066691430-L1.1/T3';...
          'ALPSRP066691460-L1.1/T3';...
          'ALPSRP066691480-L1.1/T3']];
fout = [repmat('/home/akb/Code/Matlab/PolSAR/output/20070426/',3,1) ...
        ['430'; '460';'480']];
x_c = SeaIce([624 4608],'ALOS PALSAR', plotSetting,...
        'inputDataDir', fin(3,:),  'outputDataDir', fout(3,:));
global summer
summer = 0;
%% Ground truth data
load('/media/akb/2026EF9426EF696C/raw_data/20070426/mask_070426_03.mat')
bw = bw_1 + bw_2 + bw_3 + bw_4 + bw_5 + bw_6 + bw_7 + bw_8;
bw = logical(bw);
x_c.groundtruth = bw;
clear bw bw_1  bw_2 bw_3 bw_4 bw_5 bw_6 bw_7 bw_8 bwf_1
figure
imagesc(x_c.groundtruth)
set(gca,'Ydir','normal')
colormap gray
%% Generate data
input_vector = 4;
if input_vector==4 || input_vector==3
    x_c.generateImage4Classification(input_vector);
    if input_vector == 4
        bai = 100;
        pixel_cov = bai*permute( ...
        reshape( ...
        reshape( ...
        cat(3, x_c.hh_hh, sqrt(2)*x_c.hh_hv, x_c.hh_vv, ...
        sqrt(2)*conj(x_c.hh_hv), 2*x_c.hv_hv, sqrt(2)*x_c.hv_vv, ...
        conj(x_c.hh_vv), sqrt(2)*conj(x_c.hv_vv), x_c.vv_vv), ...
        [], 9).', ...
        3,3,[]), ...
        [2, 1, 3]);
    end
else
    x_c.generateImage4Classification(input_vector,'texture',texture);
    %[~, texture] = x_c.generateImage4Classification(input_vector);
end

%% FCM
algo = 'GFCM';
distance_metric = 'wishart';
%{
wishart
squaredeuclidean
%}
num_c = 3;
num_subc = 2;
drg_m = 2;
drg_n = 2;
[label_fcm, c, u] = x_c.myFCM(x_c.im,...
                            'clusterNum', int32(num_c), ...
                            'subClusterNum', int32(num_subc), ...
                            'max_iter', int32(30), ...
                            'drg_m', drg_m, ...
                            'drg_n', drg_n, ...
                            'algo', algo, ...
                            'distance', distance_metric, ...
                            'covariance', pixel_cov);

x_c.showLabels(reshape(label_fcm,x_c.IMAGE_SIZE), num_c)
z = sprintf('%s %s, $(I,J,m,n)$ = (%i,%i,%i,%i)', algo, distance_metric, num_c, num_subc, drg_m, drg_n);
title(z,'Interpreter', 'latex')
%%
figure
scatter3(c(:,1),c(:,2),c(:,3), 32, linspecer(num_c,'qualitative'),'LineWidth',3)
set(gca,'Xlim',[0,1],'Ylim',[0,1],'Zlim',[0,1])
xlabel('$\tilde \kappa_1$','Interpreter','latex')
ylabel('$\tilde \kappa_2$','Interpreter','latex')
zlabel('$\tilde f_2$','Interpreter','latex')
%%

cmap = linspecer(num_c,'qualitative');
figure
hold on
for it = 1 : num_c
label = x_c.im(label_fcm==it,:);
%scatter(label(:,1),label(:,2), 16, cmap(5,:),'filled','LineWidth',3, 'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
histogram2(label(:,1),label(:,2),'FaceColor',cmap(it,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
end
scatter(c(:,1),c(:,2), 128, cmap, '^','filled','LineWidth',3)
hold off
%scatter(label(:,1),label(:,2), 16, cmap(4,:),'filled','LineWidth',3, 'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
set(gca,'Xlim',[0,1],'Ylim',[0,1])
xlabel('$\tilde \kappa_1$','Interpreter','latex')
ylabel('$\tilde \kappa_2$','Interpreter','latex')
hold off
%% SVM
[label_svm, model_svm] = x_c.mySVM(double(reshape(x_c.im, [x_c.IMAGE_SIZE, size(x_c.im,2)])), double(x_c.groundtruth));
save([x_c.OUTPUT_PATH, '/model_svm_(' num2str(input_vector) ')'], 'model_svm')


%% k-means
k_cluster = 3;
distance_metric = '1wishart'; 
if strcmp(distance_metric, 'wishart')
    [label_kmeans, c] = x_c.myKmeans(x_c.im,...
                            'maxClusterNum', int32(k_cluster),...
                            'max_iter', int32(11), ... % test if 10 times will cause warning
                            'replicates',int32(1), ...
                            'distance', distance_metric, ...
                            'covariance', pixel_cov);
else
    [label_kmeans, c] = kmeans(x_c.im, k_cluster,...
                        'distance','sqeuclidean',...
                        'Replicates',3,...
                        'MaxIter',50,...
                        'Options',statset('Display','final'));
end
%%
x_c.showLabels(reshape(label_kmeans,x_c.IMAGE_SIZE), k_cluster)
title(sprintf('k-means cluster number: %i',k_cluster))
set(gca,'yticklabels',[],'xticklabels',[])

%%
temp = label_kmeans;
label_kmeans(temp==4) = 2;
label_kmeans(temp==2) = 4;
%%
x_c.y_hat = reshape(label_kmeans, x_c.IMAGE_SIZE);
x_c.y_hat = x_c.y_hat ==2;

figure
imagesc(x_c.y_hat)
set(gca,'Ydir','normal','Visible','off')
colormap(gca, [0,120,0;180,100,50]/255)
plot_para('Filename',[x_c.OUTPUT_PATH '/label_070426_3_km'],'Maximize',true)
%plot_para('Filename','/home/akb/Code/Matlab/PolSAR/output/20070426/label_070426_3','Maximize',true)
%%
disp('--------------fcm-------------')
fprintf('Acc: %f\n', sum(sum(x_c.y_hat == x_c.groundtruth))/numel(x_c.groundtruth))
% Uppercase: predict, lowercase: actual
Mm = sum(sum((x_c.y_hat==1).*(x_c.groundtruth==1)))/numel(x_c.groundtruth);
Mf = sum(sum((x_c.y_hat==1).*(x_c.groundtruth==0)))/numel(x_c.groundtruth);
Fm = sum(sum((x_c.y_hat==0).*(x_c.groundtruth==1)))/numel(x_c.groundtruth);
Ff = sum(sum((x_c.y_hat==0).*(x_c.groundtruth==0)))/numel(x_c.groundtruth);
fprintf('Acc:\n %f, %f\n %f ,%f \n', Mm, Mf, Fm, Ff);
clear Mm Mf Fm Ff
%% Data augmentation
tic
[x_train, y] = x_c.myDataAugmentation(reshape(x_c.im, [x_c.IMAGE_SIZE, size(x_c.im,2)]), x_c.groundtruth,...
    'data_size', 2500, ...
    'shear_range', 0.2, ...
    'rotation_range', 45, ...
    'zoom_range', 0.35, ...
    'horizontal_flip', true, ...
    'vertical_flip', true, ...
    'filename', ['070426_3_(', num2str(input_vector), ')']);
toc

if input_vector == 4
    aa = x_train;
    x_train = aa(1:600,:,:,:);
    save(['/home/akb/Code/PolSAR_ML/data_aug/x_train_070426_3_(',num2str(input_vector), ')_1.mat'], 'x_train')
    x_train = aa(601:1200,:,:,:);
    save(['/home/akb/Code/PolSAR_ML/data_aug/x_train_070426_3_(',num2str(input_vector), ')_2.mat'], 'x_train')
    x_train = aa(1201:1800,:,:,:);
    save(['/home/akb/Code/PolSAR_ML/data_aug/x_train_070426_3_(',num2str(input_vector), ')_3.mat'], 'x_train')
    x_train = aa(1801:end,:,:,:);
    save(['/home/akb/Code/PolSAR_ML/data_aug/x_train_070426_3_(',num2str(input_vector), ')_4.mat'], 'x_train')
    x_train = aa;
    clear aa
end
%% CNN
img_h = 256;
img_w = 256;
new_size = [512 4608];
%%
it = 1;
temp_im = imresize(reshape(x_c.im, [x_c.IMAGE_SIZE, size(x_c.im,2)]), new_size, 'method', 'nearest');
x_val = single(zeros(new_size(1)/img_h*new_size(2)/img_w, img_h, img_w, size(x_c.im,2)));
disp(size(x_val))
for it = 1 : size(x_c.im,2)
    tile = im2col(temp_im(:,:,it),[img_h, img_w],'distinct');
    tile = permute(reshape(tile, img_h, img_w, []),[3,1,2]);
    for it_t = 1 : size(tile,1)
        x_val(it_t,:,:,it) = squeeze(tile(it_t,:,:));

    end
    it = it +1;
end
clear temp_im
y_val = imresize(reshape(x_c.groundtruth, x_c.IMAGE_SIZE), new_size, 'method', 'nearest');
y_val = reshape(im2col(y_val,[img_h, img_w],'distinct'), img_h, img_w, []);
y_val = permute(y_val, [3,1,2]);
save(['/home/akb/Code/PolSAR_ML/data_val/x_val_070426_3_', num2str(input_vector),'.mat'], 'x_val')
save(['/home/akb/Code/PolSAR_ML/data_val/y_val_070426_3.mat'], 'y_val')
%%
for it = 1 : 36
    figure
    imagesc(flipud(squeeze(y_val(it,:,:))))
end
%%
load(['/home/akb/Code/PolSAR_ML/output/x_c.y_hat_070426_3_', num2str(input_vector),'.mat'])
bw_cnn = imresize(reshape(x_c.groundtruth, x_c.IMAGE_SIZE), new_size, 'method', 'nearest');
temp = x_c.y_hat;
x_c.y_hat = zeros(size(bw_cnn));
if 1
    it = 1;
    for it_c = 1 : img_w : size(bw_cnn, 2)
        for it_r = 1 : img_h: size(bw_cnn,1)
            x_c.y_hat(it_r: it_r+img_h-1, it_c: it_c+img_w-1) = temp(it,:,:,1);
            it = it+1;
        end
    end
    x_c.y_hat = x_c.y_hat<0.5;
else 
    x_c.y_hat = squeeze(x_c.y_hat(:,:,:,1))<0.5;
    temp_label = imresize(x_c.groundtruth, [496, 496], 'method','nearest');
end

figure
imagesc(x_c.y_hat)
set(gca,'Ydir','normal','Visible','off')
colormap(gca, [0,120,0;180,100,50]/255)
plot_para('Filename',[x_c.OUTPUT_PATH '/label_070426_3_cnn'],'Maximize',true)

disp('--------------cnn-------------')
fprintf('Acc: %f\n', sum(sum(x_c.y_hat == bw_cnn))/numel(bw_cnn))
% Uppercase: predict, lowercase: actual
Mm = sum(sum((x_c.y_hat==1).*(bw_cnn==1)))/numel(bw_cnn);
Mf = sum(sum((x_c.y_hat==1).*(bw_cnn==0)))/numel(bw_cnn);
Fm = sum(sum((x_c.y_hat==0).*(bw_cnn==1)))/numel(bw_cnn);
Ff = sum(sum((x_c.y_hat==0).*(bw_cnn==0)))/numel(bw_cnn);
fprintf('Acc:\n %f, %f\n %f ,%f \n', Mm, Mf, Fm, Ff);