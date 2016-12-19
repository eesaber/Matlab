clear
over = 4;
step = 3;
x = linsapce(0, over * step, step); 
%% v_y = 0 , v_x change
for kk = 1 : over
  [r_PLSR(kk), a_PLSR(k), rsandb(kk), asandb(kk), t_v_x(kk), t_v_y(kk)] =  SAR((kk - 1) * step, 0);
end
figure 
subplot
%% v_y change, v_x = 0
for kk = 1 : over
   SAR(0, kk * step)
end
%% v_y = v_x 
for kk = 1 : over
   SAR(kk * step, kk * step)
end