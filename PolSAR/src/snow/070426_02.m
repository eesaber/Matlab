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
x_b = SeaIce([624 4608],'ALOS PALSAR', plotSetting,...
        'inputDataDir', fin(2,:),  'outputDataDir', fout(2,:));

%% Generate data
%im = x_b.generateImage4Classification(texture);
[im, texture] = x_b.generateImage4Classification();

%% Test sigle algorithm
algo = 'GFCM';
num_c = 4;
num_subc = 2;
drg_m = 2;
drg_n = 2;
tic 
[label_fcm, c, u] = x_b.myFCM(double(im),num_c,75,num_subc,algo,drg_m,drg_n);
toc

x_b.showLabels(reshape(label_fcm,x_b.IMAGE_SIZE), num_c)
title(sprintf('(squaredeuclidean) algorithm: %s, (I,J,m,n) = (%i,%i,%i,%i)', algo, num_c, num_subc, drg_m, drg_n))
%title([algo,' $(I,J,m,n) = ($', num2str(num_c),',',num2str(num_subc),',',...
%    num2str(drg_m),',',num2str(drg_n),')'],'Interpreter', 'latex')
set(gca,'Ydir','normal','yticklabels',[],'xticklabels',[])

%% Grid searching fuzzy parameters
algo = 'GHFCM';
for drg_m = [1.5, 3, 4]
    for drg_n = [1.5, 3, 4]
        [labels, c, u] = x_b.myFCM(double(im),5,50,10,algo,drg_m,drg_n);
        figure
        imagesc(reshape(labels,x_b.IMAGE_SIZE))
        set(gca,'Ydir','normal')
        print([num2str(drg_m), num2str(drg_n)],'-djpeg','-r300','-noui')
    end
end

%% Grid searching number of cluster center 
max_iter = 1000;
im_size = x_b.IMAGE_SIZE;
parfor num_c = 3 : 6

    for num_subc = 2 : 2
        g_name = ['(',num2str(num_c),',',num2str(num_subc),')'];
    %{
        tic 
            labels = x_b.myFCM(double(im),num_c,max_iter,num_subc,'HFCM',2,2);
            figure('name',['hfcm', g_name])
            imagesc(reshape(labels,im_size))
            set(gca,'Ydir','normal')
            print(['hfcm_',g_name],'-djpeg','-noui','-r300')
        toc
    %} 
        tic
            labels = x_b.myFCM(double(im),num_c,max_iter,num_subc,'GHFCM',2,2);
            figure('name',['ghfcm', g_name])
            imagesc(reshape(labels,im_size))
            set(gca,'Ydir','normal')
            print(['ghfcm_',g_name],'-djpeg','-noui','-r300')
        toc
    end
        g_name = ['(',num2str(num_c),')'];
    %{
        tic
            labels = x_b.myFCM(double(im),num_c,max_iter,0,'CFCM',2,2);
            figure('name',['fcm', g_name])
            imagesc(reshape(labels,im_size))
            set(gca,'Ydir','normal')
            print(['fcm_',g_name],'-djpeg','-noui','-r300')
        toc
    %}
        labels = x_b.myFCM(double(im),num_c,max_iter,0,'GFCM',2,2);
        tic 
            figure('name',['gfcm', g_name])
            imagesc(reshape(labels,im_size))
            set(gca,'Ydir','normal')
            print(['gfcm_',g_name],'-djpeg','-noui','-r300')
        toc    
    
end
%%
label_svm = x_b.mySVM(double(reshape(im, [x_b.IMAGE_SIZE, 3])), double(bw));
%%
k_cluster = 10;
[label_kmeans, ~] = kmeans(im, k_cluster,'distance','sqeuclidean', ...
                                          'Replicates',3,'MaxIter',150,...
                                          'Options',statset('Display','final'));
                                      
x_b.showLabels(reshape(label_kmeans,x_b.IMAGE_SIZE), k_cluster)
title(sprintf('k-means cluster number: %i',k_cluster))
set(gca,'yticklabels',[],'xticklabels',[])
%%
y_hat = reshape(label_fcm, x_b.IMAGE_SIZE);
y_hat = y_hat < 3;
figure
imagesc(y_hat)
set(gca,'Ydir','normal','yticklabels',[],'xticklabels',[])
colormap(gca, [0,120,0;180,100,50]/255)
plot_para('Filename','label_070428_3_fcm.jpg','Maximize',true)
%%
disp('--------------???-------------')
fprintf('Acc: %f\n', sum(sum(y_hat == bw))/numel(bw))
% Uppercase: predict, lowercase: actual
Mm = sum(sum((y_hat==1).*(bw==1)))/numel(bw);
Mf = sum(sum((y_hat==1).*(bw==0)))/numel(bw);
Fm = sum(sum((y_hat==0).*(bw==1)))/numel(bw);
Ff = sum(sum((y_hat==0).*(bw==0)))/numel(bw);
fprintf('Acc:\n %f, %f\n %f ,%f \n', Mm, Mf, Fm, Ff);
