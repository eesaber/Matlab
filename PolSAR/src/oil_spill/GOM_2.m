% This file implements the detection method of oil-spill by using PolSAR.
% Target area is in the Gulf of Mexico, mission 9
clear 
clc
chk_pw('output_GulfMexico2/area1/')
%% read data
disp('loading data...')
sub_map = 1;
if sub_map == 0
    [hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',9);
else
    [hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO('MissionNum',9,'Test','area1');
end
[N_az, N_ra] = size(hh_hh);
x = [1900, 4200, 5100];
y = [3100, 1600, 500];
pow_range = [-35 5];
size_N = numel(hh_hh);
span = hh_hh+vv_vv+2*hv_hv;
plot_set = @ Plotsetting_GOM2;
mask = @(x) ones(x,x)/x^2;
mask_size = 9;
%% 4-component decomposition
FourComp_decomp(hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv, '4decomp','Plotsetting', plot_set, 'Powrange',[-35 5])
%FourComp_decomp(temp(hh_hh), temp(hv_hv), temp(vv_vv), temp(hh_hv), temp(hh_vv), temp(hv_vv), '4decomp','Plotsetting', f_plot, 'Powrange',[-35 5])
%% Span
figure
    imagesc(10*log10(span))
    Plotsetting_GOM2(pow_range, 1,'Colorbar_unit',"(dB)")
    colormap gray
    plot_para('Filename','span', 'Maximize',true)
    %Section_GOM2(10*log10(span./(trimmean(span,3,2)*ones(1,8000))),x,y,'Span($\bar{\bar{S}}) / \mu$ (dB)' ,'span')
%%
figure
    imagesc(10*log10(hh_hh))
    Plotsetting_GOM2(pow_range, 1, 'Colorbar_unit',"(dB)")
    colormap gray
    plot_para('Ratio',[4 3 1])
    Section_GOM2(10*log10(hh_hh./(trimmean(hh_hh,3,2)*ones(1,8000))),x,y,'$|S_{hh}|^2/ \mu $ (dB)','hh_hh')
    %plot_para('Filename','S_hh', 'Maximize',true, 'Ratio',[4 3 1])
%%
f_name = 'vv_vv';
figure
    imagesc(10*log10(vv_vv))
    Plotsetting_GOM2(pow_range, 1,'Colorbar_unit',"(dB)")
    colormap gray
    plot_para('Ratio',[4 3 1])
    Section_GOM2(10*log10(vv_vv./(trimmean(vv_vv,3,2)*ones(1,8000))),x,y,'$|S_{vv}|^2/ \mu $ (dB)',f_name)
    plot_para('Ratio',[4 3 1])
    
    y_ano = [0.85 0.8; 0.5 0.55];
    annotation('arrow',0.425*ones(1,2), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','k','Linewidth',3)
    annotation('arrow',0.49*ones(1,2),  y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','r','Linewidth',3)    
    annotation('arrow',0.521*ones(1,2), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','b','Linewidth',3)    
    annotation('arrow',0.619*ones(2,1), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','k','Linewidth',3)
    annotation('arrow',0.679*ones(2,1), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','r','Linewidth',3)
    annotation('arrow',0.686*ones(2,1), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','b','Linewidth',3)
    annotation('arrow',0.279*ones(2,1), y_ano(2,:),'Headwidth',10,'Headstyle','vback2','color','k','Linewidth',3)
    annotation('arrow',0.287*ones(2,1), y_ano(2,:),'Headwidth',10,'Headstyle','vback2','color','r','Linewidth',3)
    plot_para('Maximize',true,'Filename',[f_name '_am'],'Ratio',[4 3 1])
    
%%
figure
    imagesc(10*log10(conv2(hv_hv, mask(mask_size), 'same')))
    Plotsetting_GOM2([-35 -25], 1,'Colorbar_unit',"(dB)")
    plot_para('Ratio',[4 3 1])
    colormap gray
    Section_GOM2(10*log10(hv_hv./(trimmean(hv_hv,3,2)*ones(1,8000))),x,y,'$|S_{hv}|^2$ (dB)','hv_hv')
    %plot_para('Filename','S_hv', 'Maximize',true, 'Ratio',[4 3 1])
%% Find boardar
close all
boarder = conv2(filloutliers(vv_vv(100:end,:), 0,2),mask(mask_size),'same');
%boarder = conv2(filloutliers(vv_vv, 0,2),mask(mask_size),'same');
boarder = filloutliers(boarder,0,1);
boarder = normalize(boarder,2)>0;
%BW = edge(boarder,'zerocross',0,'horizontal');
figure
imagesc(boarder);
set(gca, 'clim',[0 1], 'Ydir', 'normal')
colormap gray
figure
imagesc(10*log10(vv_vv.*~boarder))
set(gca, 'clim',[-35 5], 'Ydir', 'normal')
colormap gray

figure
imagesc(10*log10(vv_vv.*boarder))
set(gca, 'clim',[-35 5], 'Ydir', 'normal')
colormap gray
%%
temp = ones(size(vv_vv(100:end,:)));
temp(1:700,:) = 0;
figure
histogram(vv_vv(100:end,:).*boarder,'BinLimits',[1e-5,1],'Normalization','probability')
hold on
histogram(vv_vv(100:end,:).*~boarder.*temp,'BinLimits',[1e-5,1],'Normalization','probability')
hold off
xlim([0 0.03])
xlabel('$|S_{vv}|^2$', 'interpreter', 'latex')
ylabel('Pr$(x = |S_{vv}|^2)$', 'interpreter', 'latex')
legend('sea surface','oil slick')
plot_para('FileName','histo_vv','Maximize',true)
figure
histogram(hh_hh(100:end,:).*boarder,'BinLimits',[1e-5,.01],'Normalization','probability')
hold on
histogram(hh_hh(100:end,:).*~boarder.*temp,'BinLimits',[1e-5,.01],'Normalization','probability')
hold off
xlim([0 0.003])
xlabel('$|S_{hh}|^2$', 'interpreter', 'latex')
ylabel('Pr$(x = |S_{hh}|^2)$', 'interpreter', 'latex')
legend('sea surface','oil slick')
plot_para('FileName','histo_hh','Maximize',true)


temp = ones(size(vv_vv(100:end,:)));
temp2 = temp;
temp(1:700,:) = 0;
temp2(1:1000,:) = 0; 
figure
histogram(vv_vv(100:end,:)./hh_hh(100:end,:).*boarder,'BinLimits',[1e-5,50],'Normalization','probability','edgealpha',0)
hold on
histogram(vv_vv(100:end,:)./hh_hh(100:end,:).*boarder.*temp2,'BinLimits',[1e-5,50],'Normalization','probability','edgealpha',0)
histogram(vv_vv(100:end,:)./hh_hh(100:end,:).*boarder.*~temp2,'BinLimits',[1e-5,50],'Normalization','probability','edgealpha',0)
histogram(vv_vv(100:end,:)./hh_hh(100:end,:).*~boarder.*temp,'BinLimits',[1e-5,50],'Normalization','probability','edgealpha',0)
hold off
xlabel('$|S_{vv}|^2/|S_{hh}|^2$', 'interpreter', 'latex')
ylabel('$p$', 'interpreter', 'latex')
legend('sea surface','upper sea surface', 'lower sea surface', 'oil slick')
xlim([0 40])
plot_para('FileName','histo_vhratio','Maximize',true)
%%
figure
histogram(10*log10(vv_vv(100:end,:)./hh_hh(100:end,:)).*boarder,'BinLimits',[1e-5,20],'Normalization','probability','edgealpha',0)
hold on
histogram(10*log10(vv_vv(100:end,:)./hh_hh(100:end,:)).*~boarder.*temp,'BinLimits',[1e-5,20],'Normalization','probability','edgealpha',0)
hold off
xlabel('$|S_{vv}|^2/|S_{hh}|^2$', 'interpreter', 'latex')
ylabel('p', 'interpreter', 'latex')
legend('sea surface','oil slick')

%% Pauli decomposition
figure
Pauli_decomp((hh_hh+vv_vv-hh_vv-conj(hh_vv)), 2*hv_hv,(hh_hh+vv_vv+hh_vv+conj(hh_vv)),'Filename','Pauli_decomp','saibu',false)
%% Use moving Average T_mn 
temp_size = 9;
T_11 = conv2((hh_hh+vv_vv+hh_vv+conj(hh_vv))/2, mask(temp_size), 'same');
T_22 = conv2((hh_hh+vv_vv-hh_vv-conj(hh_vv))/2, mask(temp_size), 'same');
T_33 = conv2(2*hv_hv, mask(temp_size), 'same');
T_12 = conv2((hh_hh-vv_vv-hh_vv+conj(hh_vv))/2, mask(temp_size), 'same');
T_13 = conv2(hh_hv + conj(hv_vv), mask(temp_size), 'same');
T_23 = conv2(hh_hv - conj(hv_vv), mask(temp_size), 'same');
temp_T = cat(1,cat(2, reshape(T_11,[1,1,size_N]), reshape(T_12,[1,1,size_N]), reshape(T_13,[1,1,size_N])), ......
             cat(2,reshape(conj(T_12),[1,1,size_N]), reshape(T_22,[1,1,size_N]), reshape(T_23,[1,1,size_N])),.......
             cat(2,reshape(conj(T_13),[1,1,size_N]), reshape(conj(T_23),[1,1,size_N]), reshape(T_33,[1,1,size_N])));
clear T_11 T_22 T_33 T_12 T_13 T_23 
% Obtain terrain slope in azimuth and induced angle by terrain slope. 
%get_polangle(temp_T, plot_set, [-10 10]);
%% Obtain log cumulent of T
temp_C = cat(1,cat(2, reshape(hh_hh,[1,1,size_N]), sqrt(2)*reshape(hh_hv,[1,1,size_N]), reshape(hh_vv,[1,1,size_N])), ......
             cat(2,reshape(sqrt(2)*conj(hh_hv),[1,1,size_N]), reshape(2*hv_hv,[1,1,size_N]), reshape(sqrt(2)*hv_vv,[1,1,size_N])),.......
             cat(2,reshape(conj(hh_vv),[1,1,size_N]), reshape(conj(sqrt(2)*hv_vv),[1,1,size_N]), reshape(vv_vv,[1,1,size_N])));
[kai_1, kai_2] = logCumulent(temp_C);
%%
temp = zeros(size(hh_hh)) == 1;
temp(1500:2000, 4500:5000) = true;
figure
scatter(reshape(temp.*kai_1,[numel(kai_1),1]), reshape(temp.*kai_2,[numel(kai_2),1]),10,'filled','b')
set(gca,'xlim',[-27, -20], 'ylim',[-0.1 2])
hold on 
temp = zeros(size(hh_hh)) == 1;
temp(2001:2500, 1500:2000) = true;
scatter(reshape(temp.*kai_1,[numel(kai_1),1]), reshape(temp.*kai_2,[numel(kai_2),1]),10,'filled','r')
hold off
xlabel('$\kappa_1\{\bar{\bar{C}}\}$','interpreter','latex')
ylabel('$\kappa_2\{\bar{\bar{C}}\}$','interpreter','latex')
set(gca,'box','on')
plot_para('Maximize',true, 'Filename','log_cumulant')
%% Eigen-decomposition
read = 1;
f_name = 'eigen_9x9avg';
if read
    [~,lambda] = Eigen_decomp('Filename', f_name, 'Plotsetting', plot_set);
else
    [~,lambda] = Eigen_decomp('T',temp_T,'Calculate',true,'Filename',f_name, 'Plotsetting', plot_set);
end
clear temp_T
T_11 = (hh_hh+vv_vv+hh_vv+conj(hh_vv))/2;
T_22 = (hh_hh+vv_vv-hh_vv-conj(hh_vv))/2;
T_33 = 2*hv_hv;
T_12 = (hh_hh-vv_vv-hh_vv+conj(hh_vv))/2;
T_13 = hh_hv + conj(hv_vv);
T_23 = hh_hv - conj(hv_vv);

%% Indicator

%{
figure
    imagesc(H+alpha_bar/180*pi+A_1+abs(hh_vv)./sqrt(hh_hh.*vv_vv))
    Plotsetting_GOM2([0 3])
    xlabel('Azimuth (km)')
    ylabel('Range (km)') 
    plot_para('Filename','para_F', 'Maximize',true)
%}
%%
close all
f_name = 'para_confomty';
mu = 2*(-hv_hv+real(hh_vv))./(hh_hh + 2*hv_hv + vv_vv);
figure
%imagesc(mu)
imagesc(conv2(mu,ones(3,3)/9,'same'))
Plotsetting_GOM2([-0.4 1],1)
Section_GOM2(mu./(trimmean(mu,3,2)*ones(1,8000)),x,y,'$\mu_c$', f_name)
plot_para('Ratio',[4 3 1])

%delete(findall(gcf,'type','annotation'))
y_ano = [0.9 0.85; 0.5 0.55];
annotation('arrow',0.425*ones(1,2), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','k','Linewidth',3)
annotation('arrow',0.49*ones(1,2),  y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','r','Linewidth',3)    
annotation('arrow',0.521*ones(1,2), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','b','Linewidth',3)

annotation('arrow',0.619*ones(2,1), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','k','Linewidth',3)
annotation('arrow',0.679*ones(2,1), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','r','Linewidth',3)
annotation('arrow',0.686*ones(2,1), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','b','Linewidth',3)
plot_para('Filename', [f_name '_am'], 'Maximize',true,'Ratio',[4 3 1])
clear mu
%%
figure
    imagesc(-(hv_hv>abs(real(hh_vv))))
    Plotsetting_GOM2([-1 0],1)
    colormap gray; colorbar off
    
    plot_para('Filename', 'para_muller', 'Maximize',true)
%%
close all
f_name = 'para_corelation12';
figure
    imagesc(abs(T_12)./sqrt(T_11.*T_22))
    Plotsetting_GOM2([0 1],1)
    Section_GOM2(abs(T_12)./sqrt(T_11.*T_22),x,y,'$\gamma$', f_name)
    plot_para('Maximize',true,'Filename',[f_name '_am'],'Ratio',[4 3 1])
%%
f_name = 'para_moisture';
figure
    imagesc((T_11-T_22)./(T_11+T_22))
    Plotsetting_GOM2([0 1],1)
    Section_GOM2((T_11-T_22)./(T_11+T_22),x,y,f_name)   
%%
close all
f_name = 'para_roughness';
figure
    imagesc((T_22-T_33)./(T_22+T_33))
    Plotsetting_GOM2([0 1],1)
    Section_GOM2((T_22-T_33)./(T_22+T_33),x,y, '$\gamma_{r r_\perp}$',f_name)
    plot_para('Maximize',true,'Filename',[f_name '_am'],'Ratio',[4 3 1])
%%
f_name = 'para_dielectric';
figure
    %imagesc(atand((T_22+T_33)./T_11))
    imagesc((T_22+T_33)./T_11)
    Plotsetting_GOM2([0 1],1)
    Section_GOM2((T_22+T_33)./T_11,x,y, '$m (\phi, \epsilon_r)$',f_name)
    plot_para('Ratio',[4 3 1])
    
    y_ano = [0.9 0.85; 0.55 0.5; 0.33 0.28];
    annotation('arrow',0.425*ones(1,2), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','k','Linewidth',3)
    annotation('arrow',0.49*ones(1,2),  y_ano(2,:),'Headwidth',10,'Headstyle','vback2','color','r','Linewidth',3)    
    annotation('arrow',0.521*ones(1,2), y_ano(3,:),'Headwidth',10,'Headstyle','vback2','color','b','Linewidth',3)
    
    annotation('arrow',0.619*ones(2,1), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','k','Linewidth',3)
    annotation('arrow',0.679*ones(2,1), y_ano(2,:),'Headwidth',10,'Headstyle','vback2','color','r','Linewidth',3)
    annotation('arrow',0.686*ones(2,1), y_ano(3,:),'Headwidth',10,'Headstyle','vback2','color','b','Linewidth',3)
    
    plot_para('Maximize',true,'Filename',[f_name '_am'],'Ratio',[4 3 1])
%%
figure
    imagesc(sqrt(atand((T_22+T_33)./T_11)))
    Plotsetting_GOM2([0 10],1)
    Section_GOM2(sqrt(atand((T_22+T_33)./T_11)),x,y)
    %plot_para('Filename','para_braggalpha', 'Maximize',true)
%%
f_name = 'hvratio';
figure 
    imagesc(10*log10(conv2(vv_vv./hh_hh, mask(3), 'same')/9))
    %imagesc(vv_vv./hh_hh)
    Plotsetting_GOM2([0 15],1)
    Section_GOM2(10*log10(conv2(vv_vv./hh_hh, mask(3), 'same')/9),x,y,'$\sigma_{vv} / \sigma_{hh}$ (dB)',f_name)
    y_ano = [0.9 0.85];
    annotation('arrow',0.425*ones(1,2), y_ano ,'Headwidth',10,'Headstyle','vback2','color','k','Linewidth',3)
    annotation('arrow',0.49*ones(1,2),  y_ano,'Headwidth',10,'Headstyle','vback2','color','r','Linewidth',3)    
    annotation('arrow',0.521*ones(1,2), y_ano,'Headwidth',10,'Headstyle','vback2','color','b','Linewidth',3)
    
    annotation('arrow',0.619*ones(2,1), y_ano,'Headwidth',10,'Headstyle','vback2','color','k','Linewidth',3)
    annotation('arrow',0.679*ones(2,1), y_ano,'Headwidth',10,'Headstyle','vback2','color','r','Linewidth',3)
    annotation('arrow',0.686*ones(2,1), y_ano,'Headwidth',10,'Headstyle','vback2','color','b','Linewidth',3)
    annotation('arrow',0.279*ones(2,1), y_ano,'Headwidth',10,'Headstyle','vback2','color','k','Linewidth',3)
    annotation('arrow',0.287*ones(2,1), y_ano,'Headwidth',10,'Headstyle','vback2','color','r','Linewidth',3)
    plot_para('Maximize',true,'Filename',[f_name '_am'],'Ratio',[4 3 1])


%% Non-negative matrix factorization (NMF)
opt = statset('MaxIter', 100, 'Display', 'final', 'TolX', 1e-6, 'TolFun', 1e-6,'UseParallel', false);
rng('shuffle') % random seed
num = numel(hh_hh);
hh_hh = conv2(hh_hh, mask(mask_size), 'same');
vv_vv = conv2(vv_vv, mask(mask_size), 'same');
hv_hv = conv2(hv_hv, mask(mask_size), 'same');
hh_vv = conv2(hh_vv, mask(mask_size), 'same');
hh_hv = conv2(hh_hv, mask(mask_size), 'same');
hv_vv = conv2(hv_vv, mask(mask_size), 'same');

M = [reshape(hh_hh, [1 num]); reshape(hv_hv, [1 num]); reshape(vv_vv, [1 num]); reshape(abs(hh_hv), [1 num]); reshape(angle(hh_hv), [1 num]);......
    reshape(abs(hh_vv), [1 num]); reshape(angle(hh_vv), [1 num]); reshape(abs(hv_vv), [1 num]); reshape(angle(hv_vv), [1 num])];
%M = M./repmat(sum(M,1),[9,1]);
%%
[A, X] = nnmf(M, 3, 'algorithm', 'mult', 'options', opt);
%%
global im_size
X = reshape(X, [3 im_size]);
figure
    imagesc(10*log10(X(1,:,:)))
    set(gca,'Ydir','normal','clim',[-100 -50 ])