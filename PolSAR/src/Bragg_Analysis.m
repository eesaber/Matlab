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
%epsilon_sea = 80-1j*70;
epsilon_sea = 76-1j*48;
epsilon_water = 86-1j*15;
epsilon_oil = 2.3-1j*0.02;
oil_ratio = 0.8;
epsilon_mix = @(x) x*epsilon_oil + (1-x)*epsilon_sea;
k = 2*pi/(physconst('LightSpeed')/1.2575e9);
phi = 0:0.5:90;
theta = 7.2; zeta = theta;
%theta = 0; zeta = theta;
phi_i = acosd(cosd(phi+theta)*cosd(zeta));

%{
B_hh = @(x) (cosd(phi_i)-sqrt(x-sind(phi_i).^2))./(cosd(phi_i)+sqrt(x-sind(phi_i).^2));
B_vv = @(x) ((x-1)*(sind(phi_i).^2 - x*(1 + sind(phi_i).^2)))./(x*cosd(phi_i)+sqrt(x-sind(phi_i).^2)).^2;
sigma_hh = @(x) abs((sind(phi+theta)*cosd(zeta)./sind(phi_i)).^2.*B_hh(x) + (sind(zeta)./sind(phi_i)).^2.*B_vv(x)).^2; 
sigma_vv = @(x) abs((sind(phi+theta)*cosd(zeta)./sind(phi_i)).^2.*B_vv(x) + (sind(zeta)./sind(phi_i)).^2.*B_hh(x)).^2;
%}
B_hh = @(x) (cosd(phi_i)-sqrt(real(sqrt(x))-sind(phi_i).^2))./(cosd(phi_i)+sqrt(real(sqrt(x))-sind(phi_i).^2));
B_vv = @(x) ((real(sqrt(x))-1)*(sind(phi_i).^2 - real(sqrt(x))*(1 + sind(phi_i).^2)))./(real(sqrt(x))*cosd(phi_i)+sqrt(real(sqrt(x))-sind(phi_i).^2)).^2;
sigma_hh = @(x) abs((sind(phi+theta)*cosd(zeta)./sind(phi_i)).^2.*B_hh(x) + (sind(zeta)./sind(phi_i)).^2.*B_vv(x)).^2; 
sigma_vv = @(x) abs((sind(phi+theta)*cosd(zeta)./sind(phi_i)).^2.*B_vv(x) + (sind(zeta)./sind(phi_i)).^2.*B_hh(x)).^2;

figure
    plot(phi, sigma_hh(epsilon_sea),'b--', phi, sigma_vv(epsilon_sea),'b','Linewidth',3)
    hold on 
    plot(phi, sigma_hh(epsilon_mix(0.85)),'k-.', phi, sigma_vv(epsilon_mix(0.85)),'k:','Linewidth',3)
    plot(phi, sigma_hh(epsilon_mix(0.65)),'k-.', phi, sigma_vv(epsilon_mix(0.65)),'k:','Linewidth',3)
    grid on
    hold off
    
	set(gca,'Xlim',[0 80],'Ylim',[0,20])
    xlabel('$\phi$ (deg)','Interpreter', 'latex')
    ylabel('$S_{vv}$, $S_{hh}$','Interpreter', 'latex')
    plot_para('Filename','die_angle', 'Maximize',true, 'Ratio',[4 3 1])

point_size = 200;
mask = ones(1,9)/9;
figure
    plot(phi, 10*log10(sigma_vv(epsilon_sea)./sigma_hh(epsilon_sea)),'k',......
    phi, 10*log10(sigma_vv(epsilon_mix(0.85))./sigma_hh(epsilon_mix(0.85))),'k-.',.....
    phi, 10*log10(sigma_vv(epsilon_mix(0.65))./sigma_hh(epsilon_mix(0.65))),'k--','Linewidth',2.5)

    %{
    plot(phi, sigma_hh(epsilon_sea)./sigma_vv(epsilon_sea),'k--',......
    phi, sigma_hh(epsilon_mix(0.85))./sigma_vv(epsilon_mix(0.85)),'k',......
    phi, sigma_hh(epsilon_mix(0.65))./sigma_vv(epsilon_mix(0.65)),'k:',.....
    phi, sigma_hh(epsilon_oil)./sigma_vv(epsilon_oil),'g','Linewidth',2.5)
    %}
    hold on
    %{
    plot(r_phi, 10*log10(conv(vv_vv(:,4627)./hh_hh(:,4627),mask,'same')),'Linewidth',3)
    plot(r_phi, 10*log10(conv(vv_vv(:,2727)./hh_hh(:,2727),mask,'same')),'Linewidth',3)
    plot(r_phi, 10*log10(conv(vv_vv(:,4633)./hh_hh(:,4633),mask,'same')),'Linewidth',3)
    %}
    scatter(r_phi(500), 4.2, point_size, 'o','b','filled')
    scatter(r_phi(500), 4.6, point_size, 'd','b','filled')
    scatter(r_phi(1600), 8.45, point_size, 'o','r','filled')
    scatter(r_phi(1600), 9.56, point_size, 'd','r','filled')
    scatter(r_phi(3100), 8.79, point_size, 'o','k','filled')   
    scatter(r_phi(3100), 12.36, point_size, 'd','k','filled')
    
     
    hold off
    xlabel('$\phi$ (deg)','Interpreter', 'latex')
    ylabel('$\sigma_{vv} / \sigma_{hh}$ (dB)','Interpreter', 'latex')
    set(gca,'Xlim',[25 65],'Ylim',[0 15],'xGrid','on','yGrid','on')
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
plot(phi, 10*log10(sigma_vv(epsilon_sea)./sigma_hh(epsilon_sea)),'k','Linewidth',2.5)
xlabel('$\phi$ (deg)','Interpreter', 'latex')
    ylabel('$\sigma_{vv} / \sigma_{hh}$ (dB)','Interpreter', 'latex')
    set(gca,'Xlim',[5 30],'Ylim',[0 5],'xGrid','on','yGrid','on')
plot_para('Filename','hv_ratio_thm', 'Maximize',true, 'Ratio',[4 3 1])
%% Two Layer Scattering Model
psi = asind(sind(phi)*(real(physconst('LightSpeed')/sqrt(epsilon_oil))/physconst('LightSpeed')));
d_oil = 0.1e-3; % oil slick thickness 1mm 
k_oil = 2*pi/(real(physconst('LightSpeed')/sqrt(epsilon_oil))/1.2575e9);
S_hh = @(x) abs((cosd(phi)-sqrt(x-sind(phi).^2))./(cosd(phi)+sqrt(x-sind(phi).^2))).^2;
S_vv = @(x) abs((x*cosd(phi)-sqrt(x-sind(phi).^2))./(x*cosd(phi)+sqrt(x-sind(phi).^2))).^2;
R_vv = abs((S_vv(epsilon_oil)+S_vv(epsilon_sea).*exp(1j*d_oil*2*k_oil*cosd(psi)))./(1 + S_vv(epsilon_oil).*S_vv(epsilon_sea).*exp(1j*2*d_oil*k_oil*cosd(psi)))).^2;
R_hh = abs((S_hh(epsilon_oil)+S_hh(epsilon_sea).*exp(1j*d_oil*2*k_oil*cosd(psi)))./(1 + S_hh(epsilon_oil).*S_hh(epsilon_sea).*exp(1j*2*d_oil*k_oil*cosd(psi)))).^2;

figure 
plot(phi, 10*log10(R_vv./R_hh),'k', 'Linewidth',3)
ylabel('$R_{vv} / R_{hh}$ (dB)','Interpreter', 'latex')
xlabel('$\phi$ (deg)','Interpreter', 'latex')
set(gca,'xGrid','on','yGrid','on', 'xlim',[25 65])
plot_para('Filename','Rvv_Rhh_thm', 'Maximize',true, 'Ratio',[4 3 1])

%% Bragg Wavenumber
figure
    plot(phi,2*pi./((physconst('LightSpeed')/1.2575e9)/2./sind(phi)), 'k', 'Linewidth',3)
    grid on
    xlabel('$\phi$ (deg)','Interpreter', 'latex')
    ylabel('$(2 \pi)/\Lambda$ (rad/m)','Interpreter', 'latex')
    xlim([20, 65])
    plot_para('Filename','braggwavenum', 'Maximize',true, 'Ratio',[4 3 1])




%% Noise equivalent sigam zero (NESZ)
NESZ = 1.9664e-2*r_phi.^2 -1.5561*r_phi - 24.0269; % (dB)
mask = ones(9,1)/9;
figure
plot(r_phi,NESZ, 'r',r_phi, 10*log10(conv(hv_hv(:,4627),mask,'same')),'k',......
    r_phi, 10*log10(conv(hv_hv(:,2727),mask,'same')),'k--','Linewidth',3)
hold on 
plot(r_phi, 10*log10(conv(vv_vv(:,4627),mask,'same')),'g',......
    r_phi, 10*log10(conv(vv_vv(:,2727),mask,'same')),'g--','Linewidth',3)
plot(r_phi,NESZ, 'r',r_phi, 10*log10(conv(hh_hh(:,4627),mask,'same')),'b',......
    r_phi, 10*log10(conv(hh_hh(:,2727),mask,'same')),'b--','Linewidth',3)
hold off
ylabel('(dB)')
xlabel('\phi (deg)')
set(gca, 'Ylim',[-50 -10],'Xlim', [25 inf], 'Xgrid', 'on', 'Ygrid', 'on')
plot_para('Filename','NESZ_vv', 'Maximize',true, 'Ratio',[4 3 1])

figure
plot(r_phi,NESZ, 'r',r_phi, 10*log10(conv(hv_hv(:,4627),mask,'same')),'k',......
    r_phi, 10*log10(conv(hh_hh(:,4627),mask,'same')),'g',......
    r_phi, 10*log10(conv(vv_vv(:,4627),mask,'same')),'b','Linewidth',3)
ylabel('(dB)')
xlabel('\phi (deg)')
set(gca, 'Ylim',[-50 -10],'Xlim', [25 inf], 'Xgrid', 'on', 'Ygrid', 'on')
plot_para('Filename','NESZ_all', 'Maximize',true, 'Ratio',[4 3 1])
%% Eigenvalue change due to system noise
figure;
subplot(2,1,1)
plot(r_phi,lambda(:,4627,1)./sum(lambda(:,4627,:),3),'k','Linewidth',3)
set(gca, 'Xlim', [25 65], 'Xgrid', 'on', 'Ygrid', 'on')
ylabel('$p$','Interpreter','latex')
plot_para('Filename','eigenvalue_incident', 'Maximize',true)
subplot(2,1,2)
plot(r_phi, lambda(:,4627,2)./sum(lambda(:,4627,:),3),'g',......
    r_phi, lambda(:,4627,3)./sum(lambda(:,4627,:),3),'b','Linewidth',3)
ylabel('$p$','Interpreter','latex')
xlabel('\phi (deg)')
set(gca, 'Xlim', [25 65], 'Xgrid', 'on', 'Ygrid', 'on')

plot_para('Filename','eigenvalue_incident', 'Maximize',true)