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
[im, texture] = x_c.generateImage4Classification(input_vector);
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
save(['/home/akb/Code/PolSAR_ML/data/image_070426_3_(', num2str(input_vector),').mat'],'im')

%%
qq = imresize(reshape(im,624,4608,[]), [496, 496], 'method','nearest');
save('/home/akb/Code/PolSAR_ML/data/image_070426_3_(5).mat','qq')
qq = cat(3,qq, zeros(496,496));

figure
image(qq)
set(gca,'Ydir','normal')

%%
gt = imresize(bw, [496, 496], 'method','nearest');
figure
image(gt)
set(gca,'Ydir','normal')

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
label_agg = clusterdata(randi(5,10,5)+1j*randi(5,10,5),...
    'Criterion', 'inconsistent',...
    'Cutoff', 1,...
    'Linkage','ward',...
    'Maxclust',5,...
    'Distance', 'squaredeuclidean')
%%
D2 = single(zeros(size(x,3), size(c,3)));
    for it_c = 1 : size(c,3)
        ln_det_c = real(log(det(c(:,:,it_c))));
        inv_c = inv(c(:,:,it_c));
        temp = repmat(inv_c,1,1,size(D2,1)).*x;
        temp = temp(1,1,:)+temp(2,2,:)+temp(3,3,:);
        D2(:,it_c) = ln_det_c+squeeze(temp);
    end

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
x_c.myDataAugmentation(im, bw,...
    'data_size', 1000, ...
    'shear_range', 0.1, ...
    'rotation_range', 45, ...
    'zoom_range', 0.3, ...
    'horizontal_flip', true, ...
    'vertical_flip', true, ...
    'filename', '070428_3');