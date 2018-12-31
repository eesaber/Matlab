%% k-means
global kpp
kpp = true;
k_cluster = 3;
max_iter = 20;
distance_metric = 'wishart'; 
%%
[label_kpp, c_kpp] = x_b.myKmeans(x_b.im,...
                        'maxClusterNum', int32(k_cluster),...
                        'max_iter', int32(max_iter), ... % test if 10 times will cause warning
                        'replicates',int32(1), ...
                        'distance', distance_metric, ...
                        'covariance', pixel_cov);
%%
x_b.showLabels(reshape(label_kpp,x_b.IMAGE_SIZE), k_cluster)
title(sprintf('k-means cluster number: %i',k_cluster))
%%
kpp = false;
[label_Nokpp, c_Nokpp] = x_b.myKmeans(x_b.im,...
                        'maxClusterNum', int32(k_cluster),...
                        'max_iter', int32(max_iter), ... % test if 10 times will cause warning
                        'replicates',int32(1), ...
                        'distance', distance_metric, ...
                        'covariance', pixel_cov);
u_pre = 2/3/(k_cluster-1)*single(ones(size(x_b.im,1),k_cluster));
for it = 1 : k_cluster
    u_pre(label_Nokpp==it,it) = 1;
end
u_pre = u_pre./sum(u_pre,2);
%%
x_b.showLabels(reshape(label_Nokpp,x_b.IMAGE_SIZE), k_cluster)
title(sprintf('k-means cluster number: %i',k_cluster))


%% FCM
algo = 'GFCM';
num_c = k_cluster;
num_subc = 0;
%distance_metric = 'squaredeuclidean';
distance_metric = 'wishart';
drg_m = 2;
drg_n = 2;
max_iter = 25;
%% Intitial randomly
[label_fcm_r, c_r, u_r] = x_b.myFCM(x_b.im,...
                            'clusterNum', int32(num_c), ...
                            'subClusterNum', int32(num_subc), ...
                            'max_iter', int32(max_iter), ...
                            'drg_m', drg_m, ...
                            'drg_n', drg_n, ...
                            'algo', algo, ...
                            'distance', distance_metric, ...
                            'covariance', pixel_cov);
%%
x_b.showLabels(reshape(label_fcm_r,x_b.IMAGE_SIZE), num_c)
z = sprintf('%s %s, $(I,J,m,n)$ = (%i,%i,%i,%i)', algo, distance_metric, num_c, num_subc, drg_m, drg_n);
title(z,'Interpreter', 'latex')
plot_para('ratio',[16 9 1], 'Filename', [x_b.OUTPUT_PATH '/label_fcm_4'],'maximize',true)
%% Intitial by k-means results
[label_fcm_k, c_k, u_k] = x_b.myFCM(x_b.im,...
                            'clusterNum', int32(num_c), ...
                            'subClusterNum', int32(num_subc), ...
                            'max_iter', int32(max_iter), ...
                            'u', u_pre,...
                            'drg_m', drg_m, ...
                            'drg_n', drg_n, ...
                            'algo', algo, ...
                            'distance', distance_metric, ...
                            'covariance', pixel_cov);
%%
x_b.showLabels(reshape(label_fcm_k,x_b.IMAGE_SIZE), num_c)
z = sprintf('%s %s, $(I,J,m,n)$ = (%i,%i,%i,%i)', algo, distance_metric, num_c, num_subc, drg_m, drg_n);
title(z,'Interpreter', 'latex')
plot_para('ratio',[16 9 1], 'Filename', [x_b.OUTPUT_PATH '/label_fcm_4'],'maximize',true)
%%
x_b.y_hat = reshape(double(label_fcm_k), x_b.IMAGE_SIZE);
x_b.y_hat(x_b.y_hat== 2) = 0;
x_b.y_hat(x_b.y_hat== 3) = 1;

figure
imagesc(x_b.y_hat)
set(gca,'Ydir','normal','Visible','off','clim',[0 2])
colorbar
colormap(gca, [255,255,0; 0,120,0; 123 1 2;]/255)
%plot_para('Filename',[x_b.OUTPUT_PATH '/label_090811_fcm_4a'],'Maximize',true)
%%

fprintf('\n  accuracy rates: %.4f\n', sum(sum(x_b.y_hat == x_b.groundtruth))/numel(x_b.groundtruth))
% Uppercase: predict, lowercase: actual
confmat = zeros(3,3);
for it_i = 0 : 2
    for it_j = 0 : 2
       confmat(it_i+1, it_j+1) =  sum(sum((x_b.y_hat==it_i).*(x_b.groundtruth==it_j)));
    end
end
confmat = confmat/numel(x_b.groundtruth);
fprintf('%s\n',['+' repmat([repmat('-',1,7) '+' ], 1, size(confmat,1))])
fprintf(repmat(['| ' repmat('%.3f | ',1,size(confmat,1)) '\n'], 1, size(confmat,1)), confmat);
fprintf('%s\n\n',['+' repmat([repmat('-',1,7) '+' ], 1, size(confmat,1))])
%%
x_b.HAlphaDiagram(x_b.H(:,4271:end), x_b.alpha_bar(:,4271:end),'Seglinestyle','k--')
%%
figure
        histogram2(x_b.H(:), x_b.alpha_bar(:), ...
            'DisplayStyle','bar3', ...
            'YBinLimits', [0 90], ...
            'XBinLimits', [0 1], ...
            'NumBins', [150 350], ...
            'Normalization','pdf');
%xlim([0.7 1])
%ylim([20 70])
