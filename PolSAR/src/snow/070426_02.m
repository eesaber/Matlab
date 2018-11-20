%% Ground truth data
img4roi = 10*log10(x_b.vv_vv);
img4roi(img4roi>-5) = -5;
img4roi(img4roi<-25) = -25;
img4roi = rescale(img4roi);
%bw_testing = logical(zeros([x_b.IMAGE_SIZE, 10]));
%%

figure
[bw_testing(:,:,4), ver_x, ver_y] = roipoly(flipud(img4roi));
%%
bw = zeros(x_b.IMAGE_SIZE);
for it = 1 : 3, bw = bw+bw_testing(:,:,it); end
bw = bw.*bw_testing(:,:,10) - (~bw_testing(:,:,10));
bw = flipud(bw);
figure
imagesc(bw)
set(gca,'Ydir','normal','Visible','off')
colormap(gca, [0,0,0;0,120,0;180,100,50]/255)
%plot_para('Filename',[x_b.OUTPUT_PATH, '/label_070426_2'],'Maximize',true)
%%
gt = int8(bw);
save('/media/akb/2026EF9426EF696C/raw_data/20070426/mask_070426_2', 'gt', 'bw_testing')
%% load groundtruth 
load('/media/akb/2026EF9426EF696C/raw_data/20070426/mask_070426_2')
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
input_vector = 3;
%[im, ~] = x_b.generateImage4Classification(input_vector,'texture',texture);
[im, ~] = x_b.generateImage4Classification(input_vector);
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
save(['/home/akb/Code/PolSAR_ML/data/image_070426_2_(', num2str(input_vector),').mat'],'im')
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
num_c = 5;
num_subc = 2;
drg_m = 2;
drg_n = 2;
tic 
[label_fcm, c, u] = x_b.myFCM(double(im),...
                            'clusterNum', int32(num_c), ...
                            'subClusterNum', int32(num_subc), ...
                            'max_iter', int32(75), ...
                            'drg_m', drg_m, ...
                            'drg_n', drg_n, ...
                            'algo', algo, ...
                            'distance', distance_metric);
                            %'covariance', pixel_cov);
toc
%%
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
load('/home/akb/Code/PolSAR_ML/output/y_hat_070426_2.mat')
label_nn = y_test_hat;
clear y_hat_test_hat
%% SVM
load('/home/akb/Code/Matlab/PolSAR/output/20070426/480/model_svm_(4)')
%label_svm = x_b.mySVM(double(reshape(im, [x_b.IMAGE_SIZE, 3])), double(bw));
label_svm = svmpredict(zeros(numel(gt),1), double(im), model_svm);

%%
k_cluster = 3;
[label_kmeans, ~] = kmeans(im, k_cluster,'distance','sqeuclidean', ...
                                          'Replicates',3,'MaxIter',150,...
                                          'Options',statset('Display','final'));

x_b.showLabels(reshape(label_kmeans,x_b.IMAGE_SIZE), k_cluster)
title(sprintf('k-means cluster number: %i',k_cluster))
set(gca,'Visible','off')
%%
temp = label_fcm;
label_fcm(temp==1) = 4;
label_fcm(temp==4) = 1;
clear temp
%%
y_hat = reshape(label_fcm, x_b.IMAGE_SIZE);
y_hat = (y_hat ~= 4);

figure
imagesc(y_hat)
set(gca,'Ydir','normal','Visible','off')
colormap(gca, [0,120,0;180,100,50]/255)
%plot_para('Filename',[x_b.OUTPUT_PATH, '/label_070426_2_svm'],'Maximize',true)
%%
disp('--------------fcm-------------')
fprintf('Acc: %f\n', sum(sum(y_hat == bw))/sum(sum(bw>=0)))
% Uppercase: predict, lowercase: actual
Mm = sum(sum((y_hat==1).*(bw==1)))/sum(sum(bw>=0));
Mf = sum(sum((y_hat==1).*(bw==0)))/sum(sum(bw>=0));
Fm = sum(sum((y_hat==0).*(bw==1)))/sum(sum(bw>=0));
Ff = sum(sum((y_hat==0).*(bw==0)))/sum(sum(bw>=0));
fprintf('Acc:\n %f, %f\n %f ,%f \n', Mm, Mf, Fm, Ff);

