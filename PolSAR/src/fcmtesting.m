mu = [0 5];
sigma = [2 0; 0 2];
rng(1)  % For reproducibility
R_1 = mvnrnd(mu,sigma,900);
mu = [0 -5];
sigma = [2 0; 0 2];
R_2 = mvnrnd(mu,sigma,900);
R_3 = mvnrnd([-8, 0],0.5*sigma,400);
R_4 = 20*[rand(200,1), rand(200,1)]-10;
figure
plot(R_1(:,1),R_1(:,2),'bx',...
	R_2(:,1),R_2(:,2),'rx',...
	R_3(:,1),R_3(:,2),'gx',...
    R_4(:,1),R_4(:,2),'yx')
xlim([-12 6])
ylim([-10 10])
x = [R_1;R_3;R_2];
tt = PolSAR_AnalyzeTool(0,0,0,'null',true);
clusterNum = 3;
max_iter = 300;
subClusterNum = 3;

%%
algo = 'GHFCM';
[label, c] = tt.myFCM(x, clusterNum, max_iter, subClusterNum, algo);
%[center, u] = fcm(x,clusterNum,[2,100,1e-5,true]);
%{
dist = pdist2(x,center,'euclidean').^2;
sum(sum((u.'.^2).*dist))
%}
fprintf('(%f,%f) (%f,%f) \n', c(1,1,1), c(1,2,1), c(1,1,2), c(1,2,2));

l = max(u);
figure
plot(R_1(:,1), R_1(:,2),'bo',...
	R_2(:,1),R_2(:,2),'ro',...
	R_3(:,1), R_3(:,2),'go')
hold on 
%plot(center(1,1),center(1,2), 'xg',...
%	center(2,1),center(2,2), 'xc','MarkerSize',15,'LineWidth',3)
plot(c(1,1,1),c(1,2,1), 'xm', c(1,1,2),c(1,2,2), 'xm',...
    c(2,1,1),c(2,2,1),'xc',c(2,1,2),c(2,2,2), 'xc','MarkerSize',20,'LineWidth',3)
%plot(c(1,1),c(1,2), '*g',...
%	c(2,1),c(2,2), '*c','MarkerSize',15,'LineWidth',3)
scatter(x(label==1,1), x(label==1,2), 10,'b*')
scatter(x(label==2,1), x(label==2,2), 10,'r*')
scatter(x(label==3,1), x(label==3,2), 10,'g*')
%scatter(x(u(1,:)==l,1), x(u(1,:)==l,2), 10,'b*')
%scatter(x(u(2,:)==l,1), x(u(2,:)==l,2), 10,'r*')
hold off
xlim([-10 4])
ylim([-10 10])
%% 
algo = 'GHFCM'; 
im = imread('peppers.png');
t = PolSAR_AnalyzeTool(0,0,0,'null',true);
lab_I = rgb2lab(im);
ab = reshape(lab_I(:,:,2:3), [], 2);
t.IMAGE_SIZE = [size(im,1) size(im,2)];
[labels, c, u] = t.myFCM(ab,5,300,2,algo);
figure
image(im,'AlphaData',reshape(labels,size(im(:,:,1)))/3)
figure('name','testing')
imagesc(reshape(labels,size(im(:,:,1))))
%% Test all algorithm
num_c = 3; 
num_subc = 2;
max_iter = 500;
figure('name','image')
image(im)
tic
for num_c = 2 : 5
    for num_subc = 2 : min(num_c, 4)
        g_name = ['(',num2str(num_c),',',num2str(num_subc),')'];
        labels = t.myFCM(ab,num_c,max_iter,num_subc,'HFCM');
        figure('name',['hfcm', g_name])
        imagesc(reshape(labels,size(im(:,:,1))))
        drawnow 
        labels = t.myFCM(ab,num_c,max_iter,num_subc,'GHFCM');
        figure('name',['ghfcm', g_name])
        imagesc(reshape(labels,size(im(:,:,1))))
        drawnow 
    end
    g_name = ['(',num2str(num_c),')'];
    labels = t.myFCM(ab,num_c,max_iter,num_subc,'CFCM');
    figure('name',['fcm', g_name])
    imagesc(reshape(labels,size(im(:,:,1))))
    drawnow 
    labels = t.myFCM(ab,num_c,max_iter,num_subc,'GFCM');
    figure('name',['gfcm', g_name])
    imagesc(reshape(labels,size(im(:,:,1))))
    drawnow 
end
toc
%% polsar
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
im = reshape(cat(3, rescale(x_c.kai_1), rescale(x_c.kai_2),...
   rescale(x_c.kai_3), rescale(x_c.kai_4)),[],4);
im = reshape(cat(3, x_c.kai_1, x_c.kai_2,x_c.kai_3, x_c.kai_4),[],4);
%%
algo = 'GHFCM';
tic 
[labels, c, u] = x_c.myFCM(double(im),4,300,2,algo);
tac
figure
imagesc(reshape(labels,x_c.IMAGE_SIZE))
set(gca,'Ydir','normal')
%%
    [heart_scale_label, heart_scale_inst] = libsvmread('../heart_scale');

    % Split Data
    train_data = heart_scale_inst(1:150,:);
    train_label = heart_scale_label(1:150,:);
    test_data = heart_scale_inst(151:270,:);
    test_label = heart_scale_label(151:270,:);

    % Linear Kernel
    model_linear = svmtrain(train_label, train_data, '-t 0');
    [predict_label_L, accuracy_L, dec_values_L] = svmpredict(test_label, test_data, model_linear);

    % Precomputed Kernel
    model_precomputed = svmtrain(train_label, [(1:150)', train_data*train_data'], '-t 4');
    [predict_label_P, accuracy_P, dec_values_P] = svmpredict(test_label, [(1:120)', test_data*train_data'], model_precomputed);

    accuracy_L % Display the accuracy using linear kernel
    accuracy_P % Display the accuracy using precomputed kernel