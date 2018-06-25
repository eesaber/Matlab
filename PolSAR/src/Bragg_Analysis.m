%% Incident angle and Bragg wavenumber
global im_size
r_phi = atand(linspace(4602.29004, 4602.29004+22490.8262, im_size(1))/12497);
figure
    plot(r_phi,'k', 'Linewidth',3)
    xlim([1 3300])
    ylim([20 65])
    set(gca,'XTick',0:600:3300, 'Xticklabel',cellstr(int2str((0:600*5:3300*5)'/1000))')
    hold on
    scatter(500,r_phi(500),100,'o','r','filled')
    scatter(1600,r_phi(1600),100,'o','r','filled')
    scatter(3100,r_phi(3100),100,'o','r','filled')
    hold off
    xlabel('Range (km)')
    ylabel('$\phi$ (deg)','Interpreter', 'latex')
    grid on
    pos = get(gca, 'Position');
    yoffset = +0.05;
    pos(2) = pos(2) + yoffset;
    set(gca, 'Position', pos)
    plot_para('Filename','incidence_angle', 'Maximize',true, 'Ratio',[4 3 1])
%% Bragg scattering coefficient 
epsilon_sea = 80-1j*70;
epsilon_oil = 2.3-1j*0.02;
epsilon_oil = 1*epsilon_oil + 0*epsilon_sea;
k = 2*pi/(physconst('LightSpeed')/1.2575e9);
phi = 0:0.5:90;
theta = 7.6;
zeta = theta;
phi_i = acosd(cosd(phi+theta)*cosd(zeta));

B_hh = @(x) (cosd(phi_i)-sqrt(x-sind(phi_i).^2))./(cosd(phi_i)+sqrt(x-sind(phi_i).^2));
B_vv = @(x) ((x-1)*(sind(phi_i).^2 - x*(1 + sind(phi_i).^2)))./(x*cosd(phi_i)+sqrt(x-sind(phi_i).^2)).^2;
%{
B_hh_oil = (cosd(phi_i)-sqrt(epsilon_oil-sind(phi_i).^2))./(cosd(phi_i)+sqrt(epsilon_oil-sind(phi_i).^2));
B_vv_oil = ((epsilon_oil-1)*(sind(phi_i).^2 - epsilon_oil*(1 + sind(phi_i).^2)))./(epsilon_oil*cosd(phi_i)+sqrt(epsilon_oil-sind(phi_i).^2)).^2;
B_hh_sea = (cosd(phi_i)-sqrt(epsilon_sea-sind(phi_i).^2))./(cosd(phi_i)+sqrt(epsilon_sea-sind(phi_i).^2));
B_vv_sea = ((epsilon_sea-1)*(sind(phi_i).^2 - epsilon_sea*(1 + sind(phi_i).^2)))./(epsilon_sea*cosd(phi_i)+sqrt(epsilon_sea-sind(phi_i).^2)).^2;
%}
sigma_hh = @(x) abs((sind(phi+theta)*cosd(zeta)./sind(phi_i)).^2.*B_hh(x) + (sind(zeta)./sind(phi_i)).^2.*B_vv(x)).^2; 
sigma_vv = @(x) abs((sind(phi+theta)*cosd(zeta)./sind(phi_i)).^2.*B_vv(x) + (sind(zeta)./sind(phi_i)).^2.*B_hh(x)).^2;
%{
S_hh_oil = cosd(phi_i).^4.* abs((sind(phi+theta)*cosd(zeta)./sind(phi_i)).^2.*B_hh_oil + (sind(zeta)./sind(phi_i)).^2.*B_vv_oil).^2;
S_vv_oil = cosd(phi_i).^4.* abs((sind(phi+theta)*cosd(zeta)./sind(phi_i)).^2.*B_vv_oil + (sind(zeta)./sind(phi_i)).^2.*B_hh_oil).^2;
S_hh_sea = cosd(phi_i).^4.* abs((sind(phi+theta)*cosd(zeta)./sind(phi_i)).^2.*B_hh_sea + (sind(zeta)./sind(phi_i)).^2.*B_vv_sea).^2;
S_vv_sea = cosd(phi_i).^4.* abs((sind(phi+theta)*cosd(zeta)./sind(phi_i)).^2.*B_vv_sea + (sind(zeta)./sind(phi_i)).^2.*B_hh_sea).^2;
%}
figure
    a = plot(phi, sigma_hh(epsilon_sea),'b--', phi, sigma_vv(epsilon_sea),'b','Linewidth',3);
    hold on 
    b = plot(phi, sigma_hh(epsilon_oil),'k-.', phi, sigma_vv(epsilon_oil),'k:','Linewidth',3);
    grid on
    hold off
    
	set(gca,'Xlim',[0 80],'Ylim',[0,20])
    xlabel('$\phi$ (deg)','Interpreter', 'latex')
    ylabel('$S_{vv}$, $S_{hh}$','Interpreter', 'latex')
    plot_para('Filename','die_angle', 'Maximize',true, 'Ratio',[4 3 1])
%%
close all
point_size = 200;
figure
    plot(phi, 10*log10(abs(sigma_vv(epsilon_sea)./sigma_hh(epsilon_sea))),'k--',......
    phi, 10*log10(abs(sigma_vv(epsilon_oil)./sigma_hh(epsilon_oil))),'k','Linewidth',2.5)
    %plot(phi, sigma_hh(epsilon_sea)./sigma_vv(epsilon_sea),'k--',......
    %phi, sigma_hh(epsilon_oil)./sigma_vv(epsilon_oil),'k','Linewidth',2.5)
    hold on 
    scatter(r_phi(500), 4.2, point_size, 'o','b','filled')
    scatter(r_phi(500), 4.6, point_size, 'd','b','filled')
    scatter(r_phi(1600), 8.45, point_size, 'o','r','filled')
    scatter(r_phi(1600), 9.56, point_size, 'd','r','filled')
    scatter(r_phi(3100), 8.79, point_size, 'o','k','filled')   
    scatter(r_phi(3100), 12.36, point_size, 'd','k','filled')
    hold off
    xlabel('$\phi$ (deg)','Interpreter', 'latex')
    ylabel('$\sigma_{vv} / \sigma_{hh}$ (dB)','Interpreter', 'latex')
    set(gca,'Xlim',[0 80],'Ylim',[0 20],'xGrid','on','yGrid','on')
    %{
    delete(findall(gcf,'type','annotation'))
    y_ano = []
    annotation('arrow',0.425*ones(1,2), y_ano(1,:),'Headwidth',10,'Headstyle','vback2','color','k','Linewidth',3)
    annotation('arrow',0.49*ones(1,2),  y_ano(2,:),'Headwidth',10,'Headstyle','vback2','color','r','Linewidth',3)    
    annotation('arrow',0.521*ones(1,2), y_ano(3,:),'Headwidth',10,'Headstyle','vback2','color','b','Linewidth',3)
    %}
    plot_para('Filename','hv_ratio_thm', 'Maximize',true, 'Ratio',[4 3 1])
%%
figure
    plot(phi,atand(abs((B_hh-B_vv)./(B_hh+B_vv))),'k')
    grid on
    xlabel('angle (deg)')
%%  
figure
    plot(phi,2*pi./((physconst('LightSpeed')/1.2575e9)/2./sind(phi)), 'k', 'Linewidth',3)
    grid on
    xlabel('$\phi$ (deg)','Interpreter', 'latex')
    ylabel('$(2 \pi)/\Lambda$ (rad/m)','Interpreter', 'latex')
    xlim([20, 65])
    plot_para('Filename','braggwavenum', 'Maximize',true, 'Ratio',[4 3 1])
%%
idx = logical(ones(1,8000));
idx(4050:6393) = 0;
10*log10(mean(vv_vv(500, idx)))
10*log10(mean(vv_vv(500, 4050:6393)))