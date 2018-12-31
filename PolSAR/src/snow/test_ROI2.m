%% k-means
global kpp
kpp = true;
k_cluster = 3;
max_iter = 15;
distance_metric = 'wishart'; 
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
drg_m = 1.5;
drg_n = 1.5;
max_iter = 15;
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
set(gca,'Ydir','normal','Visible','off')
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
set(gca,'Ydir','normal','Visible','off')
plot_para('ratio',[16 9 1], 'Filename', [x_b.OUTPUT_PATH '/label_fcm_4'],'maximize',true)