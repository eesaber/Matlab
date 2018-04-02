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
    method = char(method(2));
    disp(['Interpolation by using ' method])
    ang = 35.6506;
    %ang = -49.5;
    [r_, c_] = size(hgt);
    r_hgt = imrotate(hgt,-ang,method);
    %% test
    temp = zeros(size(hgt));
    for a = 1 : 3
        temp = temp + 10*log10(reshape(T(a,a,:), size(hgt)));
    end
    figure
        imagesc(imrotate(temp,-ang,method));
        set(gca,'Ydir','normal','Clim',[-100 -20])
        xlabel('range','Interpreter', 'latex')
        ylabel('azimuth','Interpreter', 'latex')
    %%
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
     %{
    figure
        imagesc(rng_slp)
        set(gca,'Ydir','normal','Clim',[-1 1])
        colormap jet; colorbar; arrow('r')
        title('range gradient','Interpreter', 'latex')
        plot_para('Filename','slpRng', 'Maximize',true)
        movefile('slpRng.jpg', 'output/')
    figure
        imagesc(azi_slp)
        set(gca,'Ydir','normal','Clim',[-1 1])
        colormap jet; colorbar; arrow('r')
        title('azimuth gradient','Interpreter', 'latex')
        plot_para('Filename','slpAzi', 'Maximize',true)
        movefile('slpAzi.jpg', 'output/')
    
    
    figure
        imagesc(r_hgt)
        set(gca,'Ydir','normal')
        title('rotate DEM','Interpreter', 'latex')
        xlabel('$\leftarrow$ azimuth','Interpreter', 'latex')
        ylabel('$\leftarrow$ range','Interpreter', 'latex')
        colormap jet
        colorbar
        plot_para('Filename','hgt_rot', 'Maximize',true)
        movefile('hgt_rot.jpg', 'output/')
   
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
    % atan2(Y,X), returns values in the closed interval [-pi,pi]
    % atan(X), returns values in the closed interval [-pi/2,pi/2]
    theta = squeeze(1/4*(atan2(-2*real(T(2,3,:)), -T(2,2,:)+T(3,3,:))+pi)).';
    a = theta > pi/4;
    theta(a) = theta(a) - pi/2;
    clear a
    %{
    figure
    imagesc(reshape(theta,size(hgt)))
    caxis([-pi/4, pi/4])
    colormap gray
    plot_para('Filename','1', 'Maximize',true,'Ratio',[4 3 1])
    figure
    imagesc(-atan(-azi_slp./(-rng_slp.*cos(phi)+sin(phi))))
    caxis([-pi/2, pi/2])
    colorbar
    plot_para('Filename','2', 'Maximize',true,'Ratio',[4 3 1])
    %}
    rng_slp = reshape(rng_slp, [1, numel(rng_slp)]);
    azi_slp_est = tan(theta).*(-rng_slp.*cos(phi) + sin(phi));   
    azi_slp_est = reshape(azi_slp_est, size(hgt));
    
    %%
    
    figure
    Pauli_decomp(reshape(T(2,2,:), size(hgt)), reshape(T(3,3,:), size(hgt)),......
    reshape(T(1,1,:), size(hgt)))
    
    hold on
    %figure
        h = imagesc( abs(azi_slp-azi_slp_est));
        %title('error')
        set(gca,'Ydir','normal','Clim',[-1 0.2]) 
        colormap jet
        hold off
        set(h, 'AlphaData', (abs(azi_slp-azi_slp_est)<0.1)/2);
        arrow('y')
        plot_para('Filename','AziSlpEstErr_ovp', 'Maximize',true,'Ratio',[4 3 1])
        movefile('AziSlpEstErr_ovp.jpg', 'output/')
    
    figure
        imagesc(abs(azi_slp-azi_slp_est))
        set(gca,'Ydir','normal','Clim',[0 1])
        %title('estimation azi-slope')
        colormap jet; colorbar
        arrow([250,235,215]/256)
        plot_para('Filename','AziSlpEstErr', 'Maximize',true,'Ratio',[4 3 1])
        movefile('AziSlpEstErr.jpg', 'output/')
end
function arrow(c)
    pos_x = [0.4 0.3];
    pos_y = [0.83 0.9];
    annotation('textarrow',pos_x,pos_y,'String','range','Interpreter', 'latex','Color',c,'Headstyle','none','Linewidth',2.5,'Fontsize',24)
    pos_x = [0.3 0.4];
    pos_y = [0.9 0.83];
    annotation('textarrow',pos_x,pos_y,'String','','Color',c,'Headstyle','plain','Linewidth',2)
    pos_x = [0.3 0.27];
    pos_y = [0.9 0.8];
    annotation('textarrow',pos_x,pos_y,'String','','Color',c,'Headstyle','plain','Linewidth',2)
    pos_x = [0.27 0.3];
    pos_y = [0.8 0.9];
    annotation('textarrow',pos_x,pos_y,'String','azimuth','Interpreter', 'latex','Color',c,'Headstyle','none','Linewidth',2.5,'Fontsize',24)
end