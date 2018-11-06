%% Ground truth data
load('/media/akb/2026EF9426EF696C/raw_data/20070426/MYI_mask.mat')
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
x_c.logCumulant();
x_c.cov2coh()
[H, alpha] = x_c.eigenDecomposition(1,'eigen_9x9avg','SaveResults', false); % Eigen-decomposition
%% Generate data
im = x_c.generateImage4Classification(texture);
%[im, texture] = x_c.generateImage4Classification();
%im = reshape(cat(3, rescale(x_c.kai_1), rescale(x_c.kai_2)), [],2);
%im = reshape(cat(3, x_c.hh_hh, x_c.hv_hv, x_c.vv_vv,...
%    real(x_c.hh_hv), imag(x_c.hh_hv),...
%    real(x_c.hh_vv), imag(x_c.hh_vv),...
%    real(x_c.hv_vv), imag(x_c.hv_vv)),[],9);
temp = 10*log10(x_c.vv_vv./x_c.hh_hh);
temp(temp>3) = 3;
temp(temp<-3) = -3;
temp = rescale(temp);
%%
im = [im, H(:), temp(:)];
clear temp
%% Test sigle algorithm
algo = 'GFCM';
distance_metric = 'squaredeuclidean';
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
                            'distance', distance_metric);
toc
x_c.showLabels(reshape(label_fcm,x_c.IMAGE_SIZE), num_c)
title(sprintf('(squaredeuclidean) algorithm: %s, (I,J,m,n) = (%i,%i,%i,%i)', algo, num_c, num_subc, drg_m, drg_n))
%title([algo,' $(I,J,m,n) = ($', num2str(num_c),',',num2str(num_subc),',',...
%    num2str(drg_m),',',num2str(drg_n),')'],'Interpreter', 'latex')
set(gca,'Ydir','normal','yticklabels',[],'xticklabels',[])

%% Grid searching fuzzy parameters
algo = 'GHFCM';
for drg_m = [1.5, 3, 4]
    for drg_n = [1.5, 3, 4]
        [labels, c, u] = x_c.myFCM(double(im),5,50,10,algo,drg_m,drg_n);
        figure
        imagesc(reshape(labels,x_c.IMAGE_SIZE))
        set(gca,'Ydir','normal')
        print([num2str(drg_m), num2str(drg_n)],'-djpeg','-r300','-noui')
    end
end

%% Grid searching number of cluster center 
max_iter = 1000;
im_size = x_c.IMAGE_SIZE;
parfor num_c = 3 : 6

    for num_subc = 2 : 2
        g_name = ['(',num2str(num_c),',',num2str(num_subc),')'];
    %{
        tic 
            labels = x_c.myFCM(double(im),num_c,max_iter,num_subc,'HFCM',2,2);
            figure('name',['hfcm', g_name])
            imagesc(reshape(labels,im_size))
            set(gca,'Ydir','normal')
            print(['hfcm_',g_name],'-djpeg','-noui','-r300')
        toc
    %} 
        tic
            labels = x_c.myFCM(double(im),num_c,max_iter,num_subc,'GHFCM',2,2);
            figure('name',['ghfcm', g_name])
            imagesc(reshape(labels,im_size))
            set(gca,'Ydir','normal')
            print(['ghfcm_',g_name],'-djpeg','-noui','-r300')
        toc
    end
        g_name = ['(',num2str(num_c),')'];
    %{
        tic
            labels = x_c.myFCM(double(im),num_c,max_iter,0,'CFCM',2,2);
            figure('name',['fcm', g_name])
            imagesc(reshape(labels,im_size))
            set(gca,'Ydir','normal')
            print(['fcm_',g_name],'-djpeg','-noui','-r300')
        toc
    %}
        labels = x_c.myFCM(double(im),num_c,max_iter,0,'GFCM',2,2);
        tic 
            figure('name',['gfcm', g_name])
            imagesc(reshape(labels,im_size))
            set(gca,'Ydir','normal')
            print(['gfcm_',g_name],'-djpeg','-noui','-r300')
        toc    
    
end
%%
label_svm = x_c.mySVM(double(reshape(im, [x_c.IMAGE_SIZE, size(im,2)])), double(bw));
%%
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
%y_hat = y_hat > 2;
%%
figure
imagesc(y_hat)
set(gca,'Ydir','normal','Visible','off')
colormap(gca, [0,120,0;180,100,50]/255)
%plot_para('Filename','/home/akb/Code/Matlab/PolSAR/output/20070426/label_070426_3_svm','Maximize',true)
%plot_para('Filename','/home/akb/Code/Matlab/PolSAR/output/20070426/label_070426_3','Maximize',true)
%%
disp('--------------kmeans-------------')
fprintf('Acc: %f\n', sum(sum(y_hat == bw))/numel(bw))
% Uppercase: predict, lowercase: actual
Mm = sum(sum((y_hat==1).*(bw==1)))/numel(bw);
Mf = sum(sum((y_hat==1).*(bw==0)))/numel(bw);
Fm = sum(sum((y_hat==0).*(bw==1)))/numel(bw);
Ff = sum(sum((y_hat==0).*(bw==0)))/numel(bw);
fprintf('Acc:\n %f, %f\n %f ,%f \n', Mm, Mf, Fm, Ff);
