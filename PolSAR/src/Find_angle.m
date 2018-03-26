function Find_angle(T)
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
    %%
    method = {'bicubic','bilinear','nearest'}; % Interpolation method 
    method = char(method(1));
    disp(['Interpolation by using ' method])
    ang = -35.6506;
    %ang = -49.5;
    [r_, c_] = size(hgt);
    r_hgt = imrotate(hgt,-ang,method);
    
    refvec = [1/5.556e-05, 37.9686, -121.929];
    [~,~,gN, gE] = gradientm(hgt,refvec,referenceEllipsoid('earth'));
    [~,~,r_azi_slp,r_rng_slp] = gradientm(r_hgt,refvec,referenceEllipsoid('earth'));
    % --------How to find the boundary--------
    % 1. Create array "A" one(size(hgt))
    % 2. rotate "A" with the same angle as "hgt", then rotate "\bar{A}" inversely.
    % 3. by using "imagesc(\bar{A} >= 0)", find the boundary.
    % ----------------------------------------
    rng_slp = imrotate(r_rng_slp, ang, 'nearest');
    azi_slp = imrotate(r_azi_slp, ang, 'nearest');
    rng_slp = -rng_slp(380:880, 238:1038);
    azi_slp = -azi_slp(380:880, 238:1038);
    

    figure
        imagesc(azi_slp)
        set(gca,'Ydir','normal')
        title('azimuth gradient','Interpreter', 'latex')
        caxis([-1,1])
        colormap jet
        colorbar
        plot_para('Filename','a', 'Maximize',true)
        movefile('a.jpg', 'output/')
    figure
        imagesc(rng_slp)
        set(gca,'Ydir','normal')
        title('range gradient','Interpreter', 'latex')
        caxis([-1,1])
        colormap jet
        colorbar
        plot_para('Filename','r', 'Maximize',true)
        movefile('r.jpg', 'output/')
    %{
    figure
        imagesc(r_hgt)
        set(gca,'Ydir','normal')
        title('rotate DEM','Interpreter', 'latex')
        colormap jet
        colorbar
       plot_para('Filename','r', 'Maximize',true)
    figure 
        imagesc(gE)
        set(gca,'Ydir','normal')
        caxis([-1,1])
        colormap jet
        colorbar
        plot_para('Filename','a3', 'Maximize',true)
        movefile('gE.jpg', 'output/')
    figure 
        imagesc(gN)
        set(gca,'Ydir','normal')
        caxis([-1,1])
        colormap jet
        colorbar
        plot_para('Filename','r3', 'Maximize',true)
        movefile('gN.jpg', 'output/')
    figure 
        imagesc(slp_e)
        set(gca,'Ydir','normal')
        caxis([-1,1])
        colormap jet
        colorbar
        plot_para('Filename','a2', 'Maximize',true)
        movefile('slpE.jpg', 'output/')
    figure 
        imagesc(slp_n)
        set(gca,'Ydir','normal')
        caxis([-1,1])
        colormap jet
        colorbar
        plot_para('Filename','r2', 'Maximize',true)
        movefile('slpN.jpg', 'output/')
    figure
        Pauli_decomp(reshape(T_22,[N_az, N_ra]), reshape(T_33,[N_az, N_ra]),......
            reshape(T_11,[N_az, N_ra]), 'hgt',0, 'Contour', hgt)    
    figure
        imagesc(hgt)
        set(gca,'Ydir','normal')
        caxis([-100 600])
        colormap jet
        colorbar
        xlabel('east $\rightarrow$', 'Interpreter', 'latex')
        ylabel('north', 'Interpreter', 'latex')
        plot_para('Filename','hgt_dem', 'Maximize',true)
        movefile('hgt_dem.jpg', 'output/')
    %}
    %% Calculate the azimuth terrain slope
    phi = 20*pi/180; % Look angle
    theta = squeeze(1/4*(atan(-2*real(T(2,3,:))./(-T(2,2,:)+T(3,3,:)))+pi)).';
    a = theta > pi/4;
    theta(a) = theta(a) - pi/2;
    clear a
    figure
    imagesc(reshape(theta,size(hgt)))
    caxis([-pi/2, pi/2])
    colorbar
    plot_para('Filename','1', 'Maximize',true,'Ratio',[4 3 1])
    figure
    imagesc(-atan(-azi_slp./(-rng_slp.*cos(phi)+sin(phi))))
    caxis([-pi/2, pi/2])
    colorbar
    plot_para('Filename','2', 'Maximize',true,'Ratio',[4 3 1])
    
    rng_slp = reshape(rng_slp, [1, numel(rng_slp)]);
    azi_slp_est = -tan(theta).*(-rng_slp.*cos(phi) + sin(phi));   
    azi_slp_est = reshape(azi_slp_est, size(hgt));
    
    %%
    figure
    Pauli_decomp(reshape(T(2,2,:), size(hgt)), reshape(T(3,3,:), size(hgt)),......
    reshape(T(1,1,:), size(hgt)), 'test','contour',hgt)

    hold on
    %figure
        h = imagesc( abs(azi_slp-azi_slp_est));
        title('error')
        set(gca,'Ydir','normal')
        hold off
        set(h, 'AlphaData', (abs(azi_slp-azi_slp_est)<0.1)/2);
        caxis([-0.5,0.5])
        colormap parula
        colorbar
        plot_para('Filename','AziSlpEst', 'Maximize',true,'Ratio',[4 3 1])
        movefile('AziSlpEsterr.jpg', 'output/')
    
    figure
        imagesc(azi_slp_est)
        set(gca,'Ydir','normal') 
        %title('estimation azi-slope')
        caxis([-1,1])
        colormap jet
        colorbar
        plot_para('Filename','AziSlpEst', 'Maximize',true,'Ratio',[4 3 1])
        movefile('AziSlpEst.jpg', 'output/')
end