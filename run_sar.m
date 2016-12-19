clear
%% v_y = 0 , v_x change
over = 4;
for kk = 1 : over
  [asandb(kk), PLSR(kk), rsandb(kk), t_v_x(kk), t_v_y(kk)] =  SAR((kk - 1) * 3, 0);
end
%% v_y change, v_x = 0
for kk = 1 : over
   SAR(0, kk * 3)
end
%% v_y = v_x 
for kk = 1 : over
   SAR(kk * 3, kk * 3)
end