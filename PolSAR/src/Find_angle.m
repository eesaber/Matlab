function [theta] = Find_angle(T)
% FIND_ANGLE implements the algorithm [1] that obtain the oreintatation angle 
% from covaraince matrix C
% 
% [1] "On the Estimation of Radar Polarization Orientation Shifts Induced by Terrain Slopes"
 
    %% Read incident file and elevation data
    dir = '/media/akb/2026EF9426EF696C/raw_data/HaywardFault_CA/';
	task = 'Haywrd_23501_17114_003_171019_L090';
    if(0)
        disp('Parsing input file...')
        im_size = [21798, 13827];
        fid = fopen([dir task '_CX_01.hgt'],'r','ieee-le'); 
        hgt = single(rot90(fread(fid, im_size,'real*4')));
        fid = fopen([dir task '_CX_01.inc'],'r','ieee-le'); 
        incd = single(rot90(fread(fid, im_size,'real*4')));
        im_size = [21798*2, 13827];
        fid = fopen([dir task '_CX_01.slope'],'r','ieee-le'); 
        slp_n = single(fread(fid, im_size,'real*4'));
        slp_e = rot90(slp_n(1:2:end, :));
        slp_n = rot90(slp_n(2:2:end, :));
        
        cut = [8600,9100,14400,15200];
        hgt = hgt(cut(1):cut(2), cut(3): cut(4));
        incd = incd(cut(1):cut(2), cut(3): cut(4));
        slp_e = slp_e(cut(1):cut(2), cut(3): cut(4));
        slp_n = slp_n(cut(1):cut(2), cut(3): cut(4));
        save([dir task 'dem_1.mat'],'-v6', 'hgt','incd', 'slp_e', 'slp_n')
    else
        load([dir 'dem_1.mat']);
    end
    %% plot hgt and incident angle
    %{
    Pauli_decomp(reshape(Y(2,:),[N_az, N_ra]), reshape(Y(3,:),[N_az, N_ra]),......
    reshape(Y(1,:),[N_az, N_ra]), 'hgt',0, 'Contour', hgt)
    
    figure 
    imagesc(incd*180/pi)
    set(gca,'Ydir','normal')
    figure
    imagesc(hgt)
    set(gca,'Ydir','normal')
    caxis([-100 600])
    colormap jet
    colorbar
    plot_para('Filename','hgt_dem', 'Maximize',true)
    movefile('hgt_dem.jpg', 'output/')
    %}
    %% Calculate the range terrain slope from DEM
    %{
    temp2 = hgt +10000;
    [found, ind]= max( temp2~=0, [], 2 );
    figure
    plot(found.*ind)
    %}
    [row_, col_] = size(hgt);
    rng_slp = zeros(size(hgt));
    route = repmat([0, ones(1, 7)],[1,200]);
    % row+1 -> col-1 or col 
    % From the boundary of image, the discrete column changes are in
    % "route".
    tilt = 6.177*sqrt(2);
    % bottom
    for n = 1 : col_
        col_num = n;
        for m = 1 : row_-1
            if ~(col_num-route(m) > 0)
                break
            end
            rng_slp(m, col_num) = (hgt(m, col_num) - hgt(m+1, col_num-route(m)))/(6.177*(route(m)==0)+tilt*(route(m)==1));
            col_num = col_num - route(m);
        end
    end
    % right
    for n = 3 : row_
        col_num = col_;
        for m = n : row_-1
            if ~(col_num-route(m) > 0)
                break
            end
            rng_slp(m, col_num) = (hgt(m, col_num) - hgt(m+1, col_num-route(m-1)))/(6.177*(route(m-1)==0)+tilt*(route(m-1)==1));
            
            col_num = col_num - route(m-1);
        end
    end
    %% Calculate the azimuth terrain slope from DEM
    azi_slp = zeros(size(hgt));
    route = repmat([1,2,2,2],[1,800]);
    tilt = 6.177*sqrt(3);
    % bottom
    for n = 1 : col_
        col_num = n;
        for m = 1 : row_-1
            if ~(col_num + route(m) < col_)
                break
            end
            azi_slp(m, col_num) = (hgt(m, col_num) - hgt(m+1, col_num+route(m)))/......
                    (6.177*sqrt(3)*(route(m)==2)+6.177*sqrt(2)*(route(m)==1));
            col_num = col_num + route(m);
        end
    end
    % left 
    for n = 2 : row_
        col_num = 1;
        for m = n : row_-1
            if~(col_num + route(m) < col_)
                break
            end
            azi_slp(m, col_num) = (hgt(m, col_num) - hgt(m+1,col_num+route(m)))/......
                (6.177*sqrt(3)*(route(m)==2)+6.177*sqrt(2)*(route(m)==1));
            col_num = col_num+route(m);
        end
    end    
    %% Calculate the azimuth terrain slope
    phi = 20*pi/180; % Look angle
    theta = squeeze(1/4*(atan(-4*real(T(2,3,:))./2./(-T(2,2,:)+4*T(3,3,:)))+pi)).';
    theta(theta > pi/4) = theta(theta > pi/4) - pi/2;
    rng_slp = reshape(rng_slp, [1, numel(rng_slp)]);
    azi_slp_est = atan(tan(theta).*(-tan(rng_slp)*cos(phi) + sin(phi)));
    azi_slp_est = reshape(azi_slp_est, size(hgt));
    figure 
    imagesc(azi_slp - azi_slp_est)
    set(gca,'Ydir','normal')
    caxis([-0.5,0.5])
    colormap jet
    colorbar
    plot_para('Filename','AziSlpEst', 'Maximize',true)
    movefile('AziSlpEst.jpg', 'output/')
    figure
    imagesc(azi_slp_est)
    set(gca,'Ydir','normal') 
end