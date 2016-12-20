over = 51;
step = 1;
x = 0: step: (over - 1) * step ; 
cd ~/Matlab/output
%% v_y = 0 , v_x change
for kk = 1 : over
  [r_PLSR(kk), a_PLSR(kk), rsandb(kk), asandb(kk), t_v_x(kk), t_v_y(kk)] =  SAR((kk - 1) * step, 0);
end

my_plot(t_v_x - x,  x, '$\tilde{v}_x - v_x$','1-1')
my_plot(r_PLSR, x, 'Range PLSR','1-3')
my_plot(a_PLSR, x, 'Azimuth PLSR','1-4')
my_plot(rsandb, x, 'Range pulse width','1-5')
my_plot(asandb, x, 'Azimuth pulse width','1-6')

%% v_y change, v_x = 0
for kk = 1 : over
   [r_PLSR(kk), a_PLSR(kk), rsandb(kk), asandb(kk), t_v_x(kk), t_v_y(kk)] = SAR(0, (kk-1) * step);
end
my_plot(real(t_v_y) - x,  x, '$\tilde{v}_y - v_y$','2-2')
my_plot(r_PLSR, x, 'Range PLSR','2-3')
my_plot(a_PLSR, x, 'Azimuth PLSR','2-4')
my_plot(rsandb, x, 'Range pulse width','2-5')
my_plot(asandb, x, 'Azimuth pulse width','2-6')
%% v_y = v_x 
for kk = 1 : over
   [r_PLSR(kk), a_PLSR(kk), rsandb(kk), asandb(kk), t_v_x(kk), t_v_y(kk)] = SAR((kk-1) * step, (kk-1) * step);
end
my_plot(t_v_x - x,  x, '$\tilde{v}_x - v_x$','3-1')
my_plot(real(t_v_y) - x,  x, '$\tilde{v}_y - v_y$','3-2')
my_plot(r_PLSR, x, 'Range PLSR','3-3')
my_plot(a_PLSR, x, 'Azimuth PLSR','3-4')
my_plot(rsandb, x, 'Range pulse width','3-5')
my_plot(asandb, x, 'Azimuth pulse width','3-6')
close all