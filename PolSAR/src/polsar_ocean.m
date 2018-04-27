%% The code implements PolSAR decomposition by using non-negative matrix factorization
clear 
clc
%% read data
disp('Using real data...')
area = 'area_3';
[hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv] = Data_IO();
[N_az, N_ra] = size(hh_hh);
size_N = numel(hh_hh);
span = hh_hh+vv_vv+2*hv_hv;
% 4-component decomposition
FourComp_decomp(hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv, '4decomp')
% Span
figure
    imagesc(10*log10(span))
    Plotsetting_1([-40 0])
    xlabel('azimuth')
    ylabel('range')  
    plot_para('Filename','output/span', 'Maximize',true)
T_11 = (hh_hh+vv_vv+hh_vv+conj(hh_vv))/2;
T_22 = (hh_hh+vv_vv-hh_vv-conj(hh_vv))/2;
T_33 = 2*hv_hv;
T_12 = (hh_hh-vv_vv-hh_vv+conj(hh_vv))/2;
T_13 = hh_hv + conj(hv_vv);
T_23 = hh_hv - conj(hv_vv);
Y = [reshape(T_11, [1, size_N]); reshape(T_22, [1, size_N]); reshape(T_33, [1, size_N]); ......
    reshape(real(T_12), [1, size_N])*sqrt(2); reshape(imag(T_12), [1, size_N])*sqrt(2);
    reshape(real(T_13), [1, size_N])*sqrt(2); reshape(imag(T_13), [1, size_N])*sqrt(2);
    reshape(real(T_23), [1, size_N])*sqrt(2); reshape(imag(T_23), [1, size_N])*sqrt(2)];
clear  hh_hh hv_hv vv_vv hh_hv hh_vv hv_vv
close all
temp_T = cat(1,cat(2, reshape(T_11,[1,1,size_N]), reshape(T_12,[1,1,size_N]), reshape(T_13,[1,1,size_N])), ......
             cat(2,reshape(conj(T_12),[1,1,size_N]), reshape(T_22,[1,1,size_N]), reshape(T_23,[1,1,size_N])),.......
             cat(2,reshape(conj(T_13),[1,1,size_N]), reshape(conj(T_23),[1,1,size_N]), reshape(T_33,[1,1,size_N])));
%% Pauli decomposition
figure
Pauli_decomp(reshape(Y(2,:),[N_az, N_ra]), reshape(Y(3,:),[N_az, N_ra]),......
    reshape(Y(1,:),[N_az, N_ra]),'Filename','Pauli_decomp','saibu',false)
%Pauli_decomp((hh_hh+vv_vv-hh_vv-conj(hh_vv))/2, 2*hv_hv,(hh_hh+vv_vv+hh_vv+conj(hh_vv))/2,'Filename','Pauli_decomp')
close all
%% Eigen-decomposition
Eigen_decomp(temp_T, span)

%% Corelation
imagesc(T_12./sqrt(T_11.*T_22))
Plotsetting_1([0 1])
xlabel('azimuth')
ylabel('range') 

%% Obtain terrain slope in azimuth and induced angle by terrain slope. 
Find_angle(temp_T);
clear T_11 T_22 T_33 T_12 T_13 T_23 
%clear temp_T 
thm_angle = resize(thm_angle, [1,size_N]);

 