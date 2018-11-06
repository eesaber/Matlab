clc; clear;
%%
img4roi = H;
bw_testing = logical(zeros([x_b.IMAGE_SIZE, 10]));
figure
[bw_testing(:,:,1), ver_x, ver_y] = roipoly(flipud(img4roi));
%%
bw = zeros(x_b.IMAGE_SIZE);
for it = 1 : size(bw_testing,3), bw = bw+bw_testing(:,:,it); end
figure
imagesc(bw)

%% Reference Data
fin = ['/media/akb/2026EF9426EF696C/raw_data/20090811/ers/'...
    'SAR_IMS_1PNESA20090811_194714_00000018A149_00314_74817_0000.tif'];
fout = '/home/akb/Code/Matlab/PolSAR/output/20090811/ers';
y = PolSAR_AnalyzeTool([],'ERS-2', plotSetting,...
    'inputDataDir', fin,  'outputDataDir', fout);
y.vv_vv = flipud(rot90(y.vv_vv));
%%
im = imrotate(y.vv_vv,125,'nearest');
im = im(723:1353, 78:2955);
figure(1)
imagesc(10*log10(im))
colormap gray
set(gca,'Ydir','normal','clim', [-25 -5])
plot_para('Maximize',true)

%%
vq = int16(interp1([1 624],[1167 2958], 1:624,'linear'));
mask = ones(x_c.IMAGE_SIZE);
for n = 1 :624
    mask(n,vq(n):end) = 0;
end
figure
imagesc(10*log10((x_c.vv_vv)).*mask)
colormap gray
set(gca,'Ydir','normal','clim', [-25 -5])
plot_para('Maximize',true)

%%
height = 675;
temp = [y_b.vv_vv(1:height,2730:2849) y_c.vv_vv(1:height,1:1600)];
temp = imrotate(temp,3,'nearest','crop');
temp = temp(25:end,1:1654);
y_c.vv_vv = imresize(temp, x_c.IMAGE_SIZE,'method','nearest','Antialiasing',true);
 
clear fin fout temp y_3
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
%%
x_b.vv_vv = [x_b.vv_vv, x_c.vv_vv(:,678:end)];
x_b.hh_hh = [x_b.hh_hh, x_c.hh_hh(:,678:end)];
x_b.hv_hv = [x_b.hv_hv, x_c.hv_hv(:,678:end)];
x_b.hh_vv = [x_b.hh_vv, x_c.hh_vv(:,678:end)];
x_b.hh_hv = [x_b.hh_hv, x_c.hh_hv(:,678:end)];
x_b.hv_vv = [x_b.hv_vv, x_c.hv_vv(:,678:end)];
x_b.IMAGE_SIZE = size(x_b.vv_vv);
%%
x_b.logCumulant();
x_b.cov2coh()
[H, alpha] = x_b.eigenDecomposition(1,'eigen_9x9avg','SaveResults', false); % Eigen-decomposition
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


%% Generate data
%im = x_b.generateImage4Classification(texture);
[im, texture] = x_b.generateImage4Classification();
%%
hv_ratio = 10*log10(x_b.vv_vv./x_b.hh_hh);
hv_ratio(hv_ratio>3) = 3;
hv_ratio(hv_ratio<-3) = -3;
hv_ratio = rescale(hv_ratio);
im = [im, hv_ratio(:), H(:)];

%% Test sigle algorithm
algo = 'GFCM';
num_c = 4;
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
toc
%[c, u] = fcm(im,num_c,[drg_m, 75, 1e-5, false]);
%[~, label_fcm] = max(u,[],1);

x_b.showLabels(reshape(label_fcm,x_b.IMAGE_SIZE), num_c)
title(sprintf('(squaredeuclidean) algorithm: %s, (I,J,m,n) = (%i,%i,%i,%i)', algo, num_c, num_subc, drg_m, drg_n))
%title([algo,' $(I,J,m,n) = ($', num2str(num_c),',',num2str(num_subc),',',...
%    num2str(drg_m),',',num2str(drg_n),')'],'Interpreter', 'latex')
set(gca,'Ydir','normal','yticklabels',[],'xticklabels',[])

%%
label_svm = x_b.mySVM(double(reshape(im, [x_b.IMAGE_SIZE, size(im,3)])), double(bw));
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
