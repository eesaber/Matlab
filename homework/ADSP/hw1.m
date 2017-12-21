clear
close all
N = 21;
k = 10;
a = [0 0.18 0.22 0.5];
f = 0:0.001:0.5;
d = 0.0001;
H_d = [zeros(1,201) ones(1,300)];
% 12 exterme point 
F = [0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.23, 0.25, 0.30, 0.35, 0.4, 0.45];
W = (F - 0.22 > 0 ) + (F - 0.18 <0)/2;
w_m = (f - 0.22 > 0) + (f - 0.18 < 0)/2;
H = (F - 0.22 > 0)*2/2;
suck = (mod(1:1:12,2)-0.5)*2;

error = 100000000;
count = 1;
while 1
	s = horzcat( cos(2 * pi * F.' * linspace(0,k,k+1)), (suck./W).' )^(-1) * H.';
	err = (cos(2 * pi * f.' *linspace(0,k,k+1) ) * s(1:k+1) - H_d.' ).' .* w_m ;
	F = f( (find( (err(2:end-1) - err(3:end)).*(err(2:end-1) - err(1:end-2))>0)));
	
	if length(F) == k
		F = [f(1), F, f(end)];
	end
	if length(F) == k+1
		if abs(err(1)) > abs(err(end))
			F = [f(1), F];
		else 
			F = [F,f(end)];
		end
	end
	Err_irt(count) = max(abs(err));
	if error >= Err_irt(count)  && error - Err_irt(count)  < d
		break;
	end
	error = Err_irt(count) ;
	count = count+1;
	W = (F - 0.22 >= 0) + (F - 0.18 <=0)/2;
	H = (F - 0.22 >= 0);
end

fprintf('Error in each iteration') 
Err_irt
figure
plot(f,cos(2 * pi * f.' *linspace(0,k,k+1) ) * s(1:k+1),'k','linewidth',3)
	set(gca,'FontSize',40,'Fontname','CMU Serif Roman')
	grid on 
	pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1); 
	%export_fig Fres.jpg 
figure
stem([s(end-1:-1:2).'/2, s(1), s(2:end-1).'/2],'k','linewidth',3)
	xlim([1,21])
	set(gca,'FontSize',40,'Fontname','CMU Serif Roman')
	grid on 
	pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1); 
	export_fig Tres.jpg 
%%


