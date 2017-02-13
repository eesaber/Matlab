v_space = linspace(-20,20,41);
if ispc 
	cd D:\Code\Simu
else
	cd ~/Matlab 
end
%% v_y = 0 , v_x change
v_r = zeros(1,41);
v_rt = v_r;
for m = 1 : 41
   [v_r(m), v_rt(m), ~] = SAR_key(v_space(m), 0,0,0);
end
%%
my_plot(v_rt - v_r, v_space, '$\tilde{v}_r - v_r$(m/s)','ErrKey_vx')

%% v_x = 0 , v_y change
v_yt = zeros(1,41);
for m = 1 : 1
   [v_r(m), v_rt(m), v_yt ] = SAR_key(v_space(m), v_space(m), 0, 0);
end
%%
my_plot(v_yt - v_space, v_space, '$\tilde{v}_y - v_y$(m/s)','ErrKey_vy')