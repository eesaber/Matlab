%% Ground truth data
load('/media/akb/2026EF9426EF696C/raw_data/20070426/mask_070426_03.mat')
bw = bw_1 + bw_2 + bw_3 + bw_4 + bw_5 + bw_6 + bw_7 + bw_8;
clear bw_1  bw_2 bw_3 bw_4 bw_5 bw_6 bw_7 bw_8 bwf_1
bw = logical(bw);
figure
imagesc(bw)
set(gca,'Ydir','normal')
colormap gray
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
%% Generate data
input_vector = 3;
%[im, ~] = x_c.generateImage4Classification(input_vector,'texture',texture);
[im, ~] = x_c.generateImage4Classification(input_vector);
if 0
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
%save(['/home/akb/Code/PolSAR_ML/data/image_070426_3_(', num2str(input_vector),').mat'],'im')

%%
qq = imresize(reshape(im,624,4608,[]), [256, 256], 'method','bilinear');
figure
imagesc(qq(:,:,1))
set(gca,'Ydir','normal')
%%
qq = reshape(qq,1,496,496,size(im,2));
save('/home/akb/Code/PolSAR_ML/data/image_070426_3_(9).mat','qq')
%{
qq = cat(3,qq, zeros(496,496));

%}
%%
gt = imresize(bw, [496, 496], 'method','nearest');
figure
imagesc(gt)
set(gca,'Ydir','normal','clim',[0 1])

%%
gt = gt(:);
save('/home/akb/Code/PolSAR_ML/data/mask_070426_3_(2).mat', 'gt')
%% FCM
algo = 'GFCM';
distance_metric = 'squaredeuclidean'; 
%{
wishart
squaredeuclidean
%}
num_c = 4;
num_subc = 2;
drg_m = 2;
drg_n = 2;
tic 
[label_fcm, c, u] = x_c.myFCM(double(im),...
                            'clusterNum', int32(num_c), ...
                            'subClusterNum', int32(num_subc), ...
                            'max_iter', int32(75), ...
                            'drg_m', drg_m, ...
                            'drg_n', drg_n, ...
                            'algo', algo, ...
                            'distance', distance_metric, ...
                            'covariance', pixel_cov);
toc

x_c.showLabels(reshape(label_fcm,x_c.IMAGE_SIZE), num_c)
title(sprintf('(squaredeuclidean) algorithm: %s, (I,J,m,n) = (%i,%i,%i,%i)', algo, num_c, num_subc, drg_m, drg_n))
%title([algo,' $(I,J,m,n) = ($', num2str(num_c),',',num2str(num_sufrom skimage.transform import resizebc),',',...
%    num2str(drg_m),',',num2str(drg_n),')'],'Interpreter', 'latex')
set(gca,'Ydir','normal','yticklabels',[],'xticklabels',[])
%% SVM
[label_svm, model_svm] = x_c.mySVM(double(reshape(im, [x_c.IMAGE_SIZE, size(im,2)])), double(bw));
%%
save([x_c.OUTPUT_PATH, '/model_svm_(' num2str(input_vector) ')'], 'model_svm')
%% k-means
k_cluster = 4;
[label_kmeans, ~] = kmeans(im, k_cluster,...
                    'distance','sqeuclidean',...
                    'Replicates',3,...
                    'MaxIter',150,...
                    'Options',statset('Display','final'));
x_c.showLabels(reshape(label_kmeans,x_c.IMAGE_SIZE), k_cluster)
title(sprintf('k-means cluster number: %i',k_cluster))
set(gca,'yticklabels',[],'xticklabels',[])

%%
temp = label_kmeans;
label_kmeans(temp==4) = 2;
label_kmeans(temp==2) = 4;
%%
y_hat = reshape(label_svm, x_c.IMAGE_SIZE);
%y_hat = y_hat < 3;

figure
imagesc(y_hat)
set(gca,'Ydir','normal','Visible','off')
%colormap(gca, [0,120,0;180,100,50]/255)
%plot_para('Filename','/home/akb/Code/Matlab/PolSAR/output/20070426/label_070426_3_svm','Maximize',true)
%plot_para('Filename','/home/akb/Code/Matlab/PolSAR/output/20070426/label_070426_3','Maximize',true)
%%
disp('--------------fcm-------------')
fprintf('Acc: %f\n', sum(sum(y_hat == bw))/numel(bw))
% Uppercase: predict, lowercase: actual
Mm = sum(sum((y_hat==1).*(bw==1)))/numel(bw);
Mf = sum(sum((y_hat==1).*(bw==0)))/numel(bw);
Fm = sum(sum((y_hat==0).*(bw==1)))/numel(bw);
Ff = sum(sum((y_hat==0).*(bw==0)))/numel(bw);
fprintf('Acc:\n %f, %f\n %f ,%f \n', Mm, Mf, Fm, Ff);

%% Data augmentation
tic
[x_train, y_train] = x_c.myDataAugmentation(reshape(im, [x_c.IMAGE_SIZE, size(im,2)]), bw,...
    'data_size', 2500, ...
    'shear_range', 0.3, ...
    'rotation_range', 45, ...
    'zoom_range', 0.9, ...
    'horizontal_flip', true, ...
    'vertical_flip', true, ...
    'filename', ['070428_3_(', num2str(input_vector), ')']);
toc
%%
%aa = x_train;
x_train = aa(1:800,:,:,:);
save(['/home/akb/Code/PolSAR_ML/data_aug/x_train_070426_3_(',num2str(input_vector), ')_1.mat'], 'x_train')
x_train = aa(801:1600,:,:,:);
save(['/home/akb/Code/PolSAR_ML/data_aug/x_train_070426_3_(',num2str(input_vector), ')_2.mat'], 'x_train')
x_train = aa(1601:end,:,:,:);
save(['/home/akb/Code/PolSAR_ML/data_aug/x_train_070426_3_(',num2str(input_vector), ')_3.mat'], 'x_train')
%% CNN
cut_h = 512;
cut_w = 512;
img_h = 256;
img_w = 256;
new_size = [512 4608];
it = 1;
temp_im = imresize(reshape(im, [x_c.IMAGE_SIZE, size(im,2)]), new_size,'method', 'nearest');
qq = single(zeros(size(temp_im,1)/cut_h*size(temp_im,2)/cut_w, img_h, img_w, size(im,2)));
disp(size(qq))
for it = 1 : size(im,2)
    tile = im2col(temp_im(:,:,it),[cut_h, cut_w],'distinct');
    tile = permute(reshape(tile, cut_h, cut_w, []),[3,1,2]);
    for it_t = 1 : size(tile,1)
        qq(it_t,:,:,it) = imresize(squeeze(tile(it_t,:,:)), [img_h, img_w], 'method', 'nearest');
    end
    it = it +1;
end
%{
for it_t = 1 : size(qq,1)
    figure
    imagesc(flipud(squeeze(qq(it_t,:,:,1))))
end
%}
clear temp_im 
save(['/home/akb/Code/PolSAR_ML/data_result/tile_070426_3_', num2str(input_vector),'.mat'], 'qq')
%%
load('/home/akb/Code/PolSAR_ML/output/y_hat_070426_3.mat')
new_size = new_size/2;
if 1
    temp = y_hat(:,:,:,1);
    temp_label = imresize(reshape(bw, x_c.IMAGE_SIZE), new_size, 'method', 'nearest');
    y_hat = zeros(new_size);
    
    %y_hat = permute(y_hat,[2,3,1]);
    %y_hat = reshape(y_hat, size(temp_label));
    it = 1;
    for it_r = 1 : img_h: size(temp_label,1)
        for it_c = 1 : img_w : size(temp_label, 2)
            y_hat(it_r: it_r+img_h-1, it_c: it_c+img_w-1) = temp(it,:,:);
            it = it+1;
        end
    end
    y_hat = y_hat<0.5;
else 
    y_hat = squeeze(y_hat(:,:,:,1))<0.5;
    temp_label = imresize(bw, [496, 496], 'method','nearest');
end
figure
imagesc(y_hat)
set(gca, 'Ydir', 'normal')
plot_para('Filename', [x_c.OUTPUT_PATH, '/1'], 'Maximize', true)
%%
disp('--------------cnn-------------')
fprintf('Acc: %f\n', sum(sum(y_hat == temp_label))/numel(temp_label))
% Uppercase: predict, lowercase: actual
Mm = sum(sum((y_hat==1).*(temp_label==1)))/numel(temp_label);
Mf = sum(sum((y_hat==1).*(temp_label==0)))/numel(temp_label);
Fm = sum(sum((y_hat==0).*(temp_label==1)))/numel(temp_label);
Ff = sum(sum((y_hat==0).*(temp_label==0)))/numel(temp_label);
fprintf('Acc:\n %f, %f\n %f ,%f \n', Mm, Mf, Fm, Ff);
