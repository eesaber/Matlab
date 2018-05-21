function Find_angle(T)
% FIND_ANGLE implements the algorithm [1] that obtain the oreintatation angle 
% from covaraince matrix C
% 
% [1] "On the Estimation of Radar Polarization Orientation Shifts Induced by Terrain Slopes"
    
    %% Read incident file and elevation data
    %{
    global dir task im_size;
    if(0)
        disp('Parsing input file...')
        
        fid = fopen([dir task '_CX_01.hgt'],'r','ieee-le'); 
        hgt = rot90(single(fread(fid, im_size,'real*4')));
        %{
        fid = fopen([dir task '_CX_01.inc'],'r','ieee-le'); 
        incd = single(rot90(fread(fid, im_size,'real*4')));
        im_size = [21798*2, 13827];
        fid = fopen([dir task '_CX_01.slope'],'r','ieee-le'); 
        slp_n = single(fread(fid, im_size,'real*4'));
        slp_e = rot90(slp_n(1:2:end, :));
        slp_n = rot90(slp_n(2:2:end, :));
        %}
        %incd = incd(cut(1):cut(2), cut(3): cut(4));
        
        
    else
        prompt = 'What is the filename of .hgt ?';
        x = input(prompt);
        load([dir x '.mat']);
    end

    
    %%
    method = {'bicubic','bilinear','nearest'}; % Interpolation method 
    method = char(method(3));
    disp(['Interpolation by using ' method])
    ang = 35;
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
        ylabel('azimuth','Interpreter', 'l        set(gca,'Ydir','normal','Clim',[-1 1])
        colormap jet; colorbar
        xlabel('east')
        ylabel('north')
        plot_para('Filename','AziSlpEst', 'Maximize',true)
        movefile('AziSlpEst.jpg', 'output/')atex')
    %}
    %%
    %{
    %refvec = [1/5.556e-05, 37.9686, -121.929];
    refvec = [1/0.000055560, 37.04863700, -118.21025500 ];
    [~,~,gN, gE] = gradientm(hgt,refvec,referenceEllipsoid('earth'));
    gN = -gN;
    line = zeros(1,size(gN,1));
    for a=1 : size(gN,1)
        line(a) = find(hgt(a,:)+10000,1);
    end
    %}
     %{
    [~,~,r_azi_slp,r_rng_slp] = gradientm(r_hgt,refvec,referenceEllipsoid('earth'));
    % --------How to find the boundary--------
    % 1. Create array "A" one(size(hgt))
    % 2. rotate "A" with the same angle as "hgt", then rotate "\bar{A}" inversely.
    % 3. by using "imagesc(\bar{A} >= 0)", find the boundary.
    % ----------------------------------------
    rng_slp = imrotate(r_rng_slp, ang, 'nearest');
    azi_slp = imrotate(r_azi_slp, ang, 'nearest');
    rng_slp = -rng_slp(377:877, 237:1037);
    azi_slp = -azi_slp(377:877, 237:1037);
    %azi_slp = -azi_slp(380:880, 238:1038);
    
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
       hold on 
       %plot(sig_angle(1500,:))
       plot(sig_angle(2500,:))
       %plot(sig_angle(3000,:))
       hold off
       xlabel('azimuth')
       xlim([1 31696])
       plot_para('Filename','output/angle_sig_13', 'Maximize',true)
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
        plotset([1000 2000])
        hold on 
        plot(line,1:size(gE,1),'w', 'Linewidth', 1 )
        hold off
        xlabel('east', 'Interpreter', 'latex')
        ylabel('north', 'Interpreter', 'latex')
        plot_para('Filename','output/hgt_dem', 'Maximize',true)
      
    figure 
        imagesc(gE)
        plotset([-1 1])
        hold on 
        plot(line,1:size(gE,1),'k', 'Linewidth', 2 )
        hold off
        xlabel('east')
        ylabel('north')
        %title('range slope')
        plot_para('Filename','output/gE', 'Maximize',true)
    figure 
        imagesc(gN)
        hold on 
        plot(line,1:size(gE,1),'k', 'Linewidth', 2 )
        hold off
        plotset([-1 1])
        %title('azimuth slope')
        xlabel('east')
        ylabel('north')
        plot_para('Filename','output/gN', 'Maximize',true)
    %}
    %% Calculate the azimuth terrain slope
    %{
    [row, ~]= size(hgt);
    phi = repmat(atan((3252+linspace(1,21192,4920))/10023), [row, 1]); % Look angle
    x = input('Corresponding column range:(by a 2x1 array) \n');
    phi = phi(:,x(1):x(2));
    %}
    global im_size
    % atan2(Y,X), returns values in the closed interval [-pi,pi]
    % atan(X), returns values in the closed interval [-pi/2,pi/2]
    theta = squeeze(1/4*(atan2(-4*real(T(2,3,:)), -T(2,2,:)+T(3,3,:))+pi)).';
    a = theta > pi/4;
    theta(a) = theta(a) - pi/2;
    sig_angle = reshape(theta, im_size);
    %rng_slp = gE;
    %azi_slp = gN;
    %thm_angle = -atan(-azi_slp./(-rng_slp.*cos(phi)+sin(phi)));
    clear a theta
    figure
        imagesc(sig_angle/pi*180)
        Plotsetting_1([-15, 15])
        xlabel('Azimuth')
        ylabel('Range')
        colormap jet; colorbar
        plot_para('Filename','output/angle_sig_1', 'Maximize',true)
    %{
    figure
       hold on 
       %plot(sig_angle(1500,:))
       plot(sig_angle(2500,:))
       %plot(sig_angle(3000,:))
       hold off
       xlabel('azimuth')
       xlim([1 31696])
       plot_para('Filename','output/angle_sig_13', 'Maximize',true)
    
    figure
        imagesc(thm_angle/pi*180)
        plotset([-45, 45])
        xlabel('east')
        ylabel('north')
        colormap jet; colorbar
        plot_para('Filename','output/angle_thm', 'Maximize',true)
    close all
    rng_slp = reshape(rng_slp, [1, numel(rng_slp)]);
    phi = reshape(phi, [1, numel(phi)]);
    azi_slp_est = tan(theta).*(-rng_slp.*cos(phi) + sin(phi));   
    azi_slp_est = reshape(azi_slp_est, size(hgt));
    %}
    %%
    %{
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
        %arrow('y')
        plot_para('Filename','AziSlpEstErr_ovp', 'Maximize',true,'Ratio',[4 3 1])
        movefile('AziSlpEstErr_ovp.jpg', 'output/')
    
    figure
        imagesc(abs(azi_slp-azi_slp_est))
        Plotsetting([0 1])
        xlabel('east')
        ylabel('north')
        plot_para('Filename','output/AziSlpEstErr', 'Maximize',true)
    
    figure
        imagesc(azi_slp_est)
        Plotsetting_1([-1 1])
        xlabel('east')
        ylabel('north')
        plot_para('Filename','output/AziSlpEst', 'Maximize',true)
    %}
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