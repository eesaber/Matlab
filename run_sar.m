over = 51;
step = 1;
x = 0: step: (over - 1) * step ; 
%cd ~/Matlab/output
cd D:\Code\Simu
%% v_y = 0 , v_x change
for kk = 1 : over
  [r_PLSR(kk), a_PLSR(kk), rsandb(kk), asandb(kk), t_v_x(kk), t_v_y(kk)] =  SAR((kk - 1) * step, 0);
end
%%
close all
my_plot(t_v_x - x, x, '$\tilde{v}_x - v_x$(m/s)','errx_y0')
my_plotyy(r_PLSR, a_PLSR, x, 'PSLR$_\tau$(dB)', 'PSLR$_\eta$(dB)','pslr_y0')
my_plotyy(rsandb, asandb, x, 'IRW$_\tau$ (m)', 'IRW$_\eta$ (m)','sand_y0')

%% v_y change, v_x = 0
for kk = 1 : over
   [r_PLSR(kk), a_PLSR(kk), rsandb(kk), asandb(kk), t_v_x(kk), t_v_y(kk)] = SAR(0, (kk-1) * step);
end
%%
close all
my_plot(real(t_v_y) - x, x,'$\tilde{v}_y - v_y$(m/s)','erry_x0')
my_plotyy(r_PLSR, a_PLSR, x, 'PSLR$_\tau$(dB)', 'PSLR$_\eta$(dB)','pslr_x0')
my_plotyy(rsandb, asandb, x, 'IRW$_\tau$ (m)', 'IRW$_\eta$ (m)','sand_x0')

%% v_y = v_x 
for kk = 1 : over
   [r_PLSR(kk), a_PLSR(kk), rsandb(kk), asandb(kk), t_v_x(kk), t_v_y(kk)] = SAR((kk-1) * step, (kk-1) * step);
end
%%
close all
my_plotyy(t_v_x - x, real(t_v_y) - x, x, '$\tilde{v_x} - v_x$(m/s)', '$\tilde{v}_y - v_y$(m/s)','errx_')
my_plotyy(r_PLSR, a_PLSR, x, 'PSLR$_\tau$(dB)', 'PSLR$_\eta$(dB)','pslr_')
my_plotyy(rsandb, asandb, x, 'IRW$_\tau$ (m)', 'IRW$_\eta$ (m)','sand_')
