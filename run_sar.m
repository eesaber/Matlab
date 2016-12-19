over = 6;
step = 3;
x = 0: step: (over - 1) * step ; 
%% v_y = 0 , v_x change
for kk = 1 : over
  [r_PLSR(kk), a_PLSR(kk), rsandb(kk), asandb(kk), t_v_x(kk), t_v_y(kk)] =  SAR((kk - 1) * step, 0);
end
my_plot(t_v_x - x,  x, '$\tilde{v}_x - v_x$')
my_plot(r_PLSR, x, 'Range PLSR')
my_plot(a_PLSR, x, 'Azimuth PLSR')
my_plot(rsandb, x, 'Range pulse width')
my_plot(asandb, x, 'Azimuth pulse width')

%% v_y change, v_x = 0
for kk = 1 : over
   SAR(0, (kk-1) * step)
end
my_plot(t_v_y - x,  x, '$\tilde{v}_y - v_y$')
my_plot(r_PLSR, x, 'Range PLSR')
my_plot(a_PLSR, x, 'Azimuth PLSR')
my_plot(rsandb, x, 'Range pulse width')
my_plot(asandb, x, 'Azimuth pulse width')
%% v_y = v_x 
for kk = 1 : over
   SAR((kk-1) * step, (kk-1) * step)
end
my_plot(t_v_x - x,  x, '$\tilde{v}_x - v_x$')
my_plot(t_v_y - x,  x, '$\tilde{v}_y - v_y$')
my_plot(r_PLSR, x, 'Range PLSR')
my_plot(a_PLSR, x, 'Azimuth PLSR')
my_plot(rsandb, x, 'Range pulse width')
my_plot(asandb, x, 'Azimuth pulse width')