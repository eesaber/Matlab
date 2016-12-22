over = 51;
step = 1;
x = 0: step: (over - 1) * step ; 
cd ~/Matlab/output
%% v_y = 0 , v_x change
for kk = 1 : over
  [r_PLSR(kk), a_PLSR(kk), rsandb(kk), asandb(kk), t_v_x(kk), t_v_y(kk)] =  SAR((kk - 1) * step, 0);
end

my_plot(t_v_x - x, 0, x, '$\tilde{v}_x - v_x$ (m)','errx_vyfix')
my_plot(r_PLSR, a_PLSR, x, 'PLSR (dB)','plsr_vyfix')
my_plot(rsandb, asandb, x, 'Pulse width (m)','sand_vyfix')

%% v_y change, v_x = 0
for kk = 1 : over
   [r_PLSR(kk), a_PLSR(kk), rsandb(kk), asandb(kk), t_v_x(kk), t_v_y(kk)] = SAR(0, (kk-1) * step);
end
my_plot(real(t_v_y) - x, 0, x, '$\tilde{v}_y - v_y$ (m/s)','erry_vxfix')
my_plot(r_PLSR, a_PLSR, x, 'PLSR (dB)','plsr_vxfix')
my_plot(rsandb, asandb, x, 'Pulse width (m)','sand_vxfix')
%% v_y = v_x 
for kk = 1 : over
   [r_PLSR(kk), a_PLSR(kk), rsandb(kk), asandb(kk), t_v_x(kk), t_v_y(kk)] = SAR((kk-1) * step, (kk-1) * step);
end
my_plot(t_v_x - x, real(t_v_y) - x, x, '$\tilde{v} - v$ (m/s)','errx_')
my_plot(r_PLSR, a_PLSR, x, 'PLSR (dB)','plsr_')
my_plot(rsandb, asandb, x, 'Pulse width (m)','sand_')
