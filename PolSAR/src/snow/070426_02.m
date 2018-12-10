%% Ground truth data
img4roi = 10*log10(x_b.vv_vv);
img4roi(img4roi>-5) = -5;
img4roi(img4roi<-25) = -25;
img4roi = flipud(rescale(img4roi));
%bw_testing = logical(zeros([x_b.IMAGE_SIZE, 10]));
%%
imageSegmenter(img4roi)
%%
figure
[bw_testing(:,:,4), ver_x, ver_y] = roipoly(flipud(img4roi));
%%
bw = zeros(x_b.IMAGE_SIZE);
for it = 1 : 3, bw = bw+bw_testing(:,:,it); end
bw = bw.*bw_testing(:,:,10) - (~bw_testing(:,:,10));
bw = int8(flipud(bw));
figure
imagesc(bw)
set(gca,'Ydir','normal','Visible','off')
colormap(gca, [0,0,0;0,120,0;180,100,50]/255)
plot_para('Filename',[x_b.OUTPUT_PATH, '/label_070426_2'],'Maximize',true)
save('/media/akb/2026EF9426EF696C/raw_data/20070426/mask_070426_2', 'bw', 'bw_testing')
%% load groundtruth 
load('/media/akb/2026EF9426EF696C/raw_data/20070426/mask_070426_2')

clear  bw_testing
%%
figure
imagesc(bw)
set(gca,'Ydir','normal','Visible','off')
colormap(gca, [0,0,0;0,120,0;180,100,50]/255)
%% Data IO
plotSetting = @ Plotsetting_dummy;
fin = [repmat('/media/akb/2026EF9426EF696C/raw_data/20070426/',3,1) ...
         ['ALPSRP066691430-L1.1/T3';...
          'ALPSRP066691460-L1.1/T3';...
          'ALPSRP066691480-L1.1/T3']];
fout = [repmat('/home/akb/Code/Matlab/PolSAR/output/20070426/',3,1) ...
        ['430'; '460';'480']];
x_b = SeaIce([624 4608],'ALOS PALSAR', plotSetting,...
        'inputDataDir', fin(2,:),  'outputDataDir', fout(2,:));
%% Generate data
input_vector = 1;
%[im, ~] = x_b.generateImage4Classification(input_vector,'texture',texture);
[im, ~] = x_b.generateImage4Classification(input_vector);
%[im, texture] = x_b.generateImage4Classification(input_vector);

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
%save(['/home/akb/Code/PolSAR_ML/data/image_070426_2_(', num2str(input_vector),').mat'],'im')
%%
qq = imresize(reshape(im,624,4608,[]), [496, 496], 'method','nearest');
save('/home/akb/Code/PolSAR_ML/data/image_070426_2_(9).mat','qq')
qq = cat(3,qq, zeros(496,496));

figure
image(qq)
set(gca,'Ydir','normal')

%% FCM
algo = 'GFCM';
distance_metric = 'squaredeuclidean'; 
%{
wishart
squaredeuclidean
%}
num_c = 2;
num_subc = 2;
drg_m = 2;
drg_n = 2;
tic 
[label_fcm, c, u] = x_b.myFCM(im,...
                            'clusterNum', int32(num_c), ...
                            'subClusterNum', int32(num_subc), ...
                            'max_iter', int32(40), ...
                            'drg_m', drg_m, ...
                            'drg_n', drg_n, ...
                            'algo', algo, ...
                            'distance', distance_metric,...
                            'covariance', pixel_cov);
toc

x_b.showLabels(reshape(label_fcm,x_b.IMAGE_SIZE), num_c)
z = sprintf('%s %s, $(I,J,m,n)$ = (%i,%i,%i,%i)', algo, distance_metric, num_c, num_subc, drg_m, drg_n);
title(z,'Interpreter', 'latex')
%set(gca,'Ydir','normal','Visible','off')
%%
figure
scatter3(c(:,1),c(:,2),c(:,3), 32, linspecer(num_c,'qualitative'),'LineWidth',3)
set(gca,'Xlim',[0,1],'Ylim',[0,1],'Zlim',[0,1])
xlabel('$\tilde \kappa_1$','Interpreter','latex')
ylabel('$\tilde \kappa_2$','Interpreter','latex')
zlabel('$\tilde f_2$','Interpreter','latex')
%% NN
load('/home/akb/Code/PolSAR_ML/output/y_hat_070426_2_4.mat')
label_nn = y_test_hat;
clear y_hat_test_hat
%% SVM
load('/home/akb/Code/Matlab/PolSAR/output/20070426/480/model_svm_(1)')
%label_svm = x_b.mySVM(double(reshape(im, [x_b.IMAGE_SIZE, 3])), double(bw));
disp(datetime)
tic
label_svm = svmpredict(zeros(numel(x_b.vv_vv),1), double(im), model_svm);
toc
%%
k_cluster = 2;
distance_metric = 1;'wishart'; 
if strcmp(distance_metric, 'wishart')
    label_kmeans = x_b.myKmeans(im,...
                            'maxClusterNum', int32(k_cluster),...
                            'max_iter', int32(7), ... % test if 10 times will cause warning
                            'replicates',int32(1), ...
                            'distance', distance_metric, ...
                            'covariance', pixel_cov);
else
    [label_kmeans, ~] = kmeans(im, k_cluster,...
                        'distance','sqeuclidean',...
                        'Replicates',3,...
                        'MaxIter',30,...
                        'Options',statset('Display','final'));
end
%%
x_b.showLabels(reshape(label_kmeans,x_b.IMAGE_SIZE), k_cluster)
title(sprintf('k-means cluster number: %i',k_cluster))
set(gca,'Visible','off')
%%
temp = label_kmeans;
label_kmeans(temp==2) = 4;
label_kmeans(temp==4) = 2;
clear temp

%%
y_hat = reshape(label_kmeans, x_b.IMAGE_SIZE);
%y_hat = ~((y_hat==1)+(y_hat==3));
y_hat = (y_hat==2);

figure
imagesc(y_hat)
set(gca,'Ydir','normal','Visible','off')

colormap(gca, [0,120,0;180,100,50]/255)
plot_para('Filename',[x_b.OUTPUT_PATH, '/label_070426_2_kmeans'],'Maximize',true)
%%
disp('--------------kmeans-------------')
fprintf('Acc: %f\n', sum(sum(y_hat == bw))/sum(sum(bw>=0)))
% Uppercase: predict, lowercase: actual
Mm = sum(sum((y_hat==1).*(bw==1)))/sum(sum(bw>=0));
Mf = sum(sum((y_hat==1).*(bw==0)))/sum(sum(bw>=0));
Fm = sum(sum((y_hat==0).*(bw==1)))/sum(sum(bw>=0));
Ff = sum(sum((y_hat==0).*(bw==0)))/sum(sum(bw>=0));
fprintf('Acc:\n %f, %f\n %f ,%f \n', Mm, Mf, Fm, Ff);

%% CNN
img_h = 256;
img_w = 256;
new_size = [512 4608];
it = 1;
%%
temp_im = imresize(reshape(im, [x_b.IMAGE_SIZE, size(im,2)]), new_size, 'method', 'nearest');
x_val = single(zeros(new_size(1)/img_h*new_size(2)/img_w, img_h, img_w, size(im,2)));

disp(size(x_val))
for it = 1 : size(im,2)
    tile = im2col(temp_im(:,:,it),[img_h, img_w],'distinct');
    tile = permute(reshape(tile, img_h, img_w, []),[3,1,2]);
    for it_t = 1 : size(tile,1)
        x_val(it_t,:,:,it) = squeeze(tile(it_t,:,:));

    end
    it = it +1;
end
clear temp_im 
save(['/home/akb/Code/PolSAR_ML/data_val/x_val_070426_2_', num2str(input_vector),'.mat'], 'x_val')
%save(['/home/akb/Code/PolSAR_ML/data_val/x_val_070426_2_6.mat'], 'x_val')
%%
load(['/home/akb/Code/PolSAR_ML/output/y_hat_070426_2_', num2str(input_vector),'.mat'])
bw_cnn = imresize(reshape(bw, x_b.IMAGE_SIZE), new_size, 'method', 'nearest');
temp = y_hat;
y_hat = zeros(size(bw_cnn));
if 1
    it = 1;
    for it_c = 1 : img_w : size(bw_cnn, 2)
        for it_r = 1 : img_h: size(bw_cnn,1)
            y_hat(it_r: it_r+img_h-1, it_c: it_c+img_w-1) = temp(it,:,:,1);
            it = it+1;
        end
    end
    y_hat = y_hat<0.5;
else 
    y_hat = squeeze(y_hat(:,:,:,1))<0.5;
    bw_cnn = imresize(bw, [496, 496], 'method','nearest');
end
figure
imagesc(y_hat)
set(gca,'Ydir','normal','Visible','off')
colormap(gca, [0,120,0;180,100,50]/255)
plot_para('Filename',[x_b.OUTPUT_PATH '/label_070426_2_cnn'],'Maximize',true)

disp('--------------cnn-------------')
fprintf('Acc: %f\n', sum(sum(y_hat == bw_cnn))/sum(sum(bw_cnn>=0)))
% Uppercase: predict, lowercase: actual
Mm = sum(sum((y_hat==1).*(bw_cnn==1)))/sum(sum(bw_cnn>=0));
Mf = sum(sum((y_hat==1).*(bw_cnn==0)))/sum(sum(bw_cnn>=0));
Fm = sum(sum((y_hat==0).*(bw_cnn==1)))/sum(sum(bw_cnn>=0));
Ff = sum(sum((y_hat==0).*(bw_cnn==0)))/sum(sum(bw_cnn>=0));
fprintf('Acc:\n %f, %f\n %f ,%f \n', Mm, Mf, Fm, Ff);