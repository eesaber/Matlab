function [phase_] = ambi_re(qq, lamda_, para_)
	% Note that the frequency is f_a > f_b > f_c 
	% fprintf('Note that the frequency is (f_a,f_b,f_c) f_a > f_b > f_c \n')
	
	d_a = para_(1); R_0 = para_(2);
	v_p =para_(3); v_y=para_(4); 
	y_0 = para_(5);
	N = 20;
	T = 90;
	N_down = min(floor(4*pi./lamda_*d_a/2/R_0*(-T*(v_p - v_y) /2 + y_0) / 2 /pi));
	N_up = max(floor(4*pi./lamda_*d_a/2/R_0*(T*(v_p - v_y) /2 + y_0) / 2 /pi));
	
	space = 2*pi*(N_down: N_up);
	k = repmat(qq.',1,length(space)) + lamda_.' * space ;
	k_bar = [k(:,1),  diff(k,1,2)];
	%% Plot the 
	x_l = 6;
	%{
	figure 
	subplot('position',[0.1 0.5 0.8 0.2]);
	%subplot('position',[0.1 0.3375 0.8 0.1125]);
	barh(k_bar,1,'stack','FaceColor','none','Linewidth',1.25,'EdgeColor','k'...
	,'BaseValue', -100, 'ShowBaseLine', 'on');
	set(gca,'TickLength',[0 0])
	xlabel('$(\beta_\alpha \tilde{y}_\alpha + 2 N \pi)\lambda_\alpha$','Interpreter', 'latex')
	yticklabels({'\lambda_4 ', '\lambda_3', '\lambda_2', '\lambda_1'})
	xlim([-x_l,x_l])
	ylim([0.5 numel(k)/length(k)+0.5])
	hold all
	x_l = 2*x_l*(N>0) - x_l;
	%}
	amb_flag = false;
	sa_ = 100 ;
	[temp, ind_] = sort(reshape(k,1,numel(k)));
	row_ = reshape( (1:numel(k)/length(k)).' *ones(1,length(space)), 1, numel(k)); 
	row_ = row_(ind_);
	for i = 1 : numel(k)- numel(k)/length(k)
		if temp(i+numel(k)/length(k)-1) - temp(i) < 0.05 && sum(diff(sort(row_(i: i + numel(k)/length(k) - 1))) == 0) == 0
			if temp(i+numel(k)/length(k)-1) - temp(i) <= sa_
				%if amb_flag
				%	fprintf('DANGER! AMBIGUITY HAPPEN')
				%end
				%assert(amb_flag, 'DANGER AMBIGUITY HAPPEN')
				sa_ = temp(i+numel(k)/length(k)-1) - temp(i) ;
				phase_ = mean(temp(i: i+numel(k)/length(k)-1));
				amb_flag = true;
				for j = 0 : 3
					plot([temp(i+j) temp(i+j)],[row_(i+j)-0.5 row_(i+j)+0.5],'Color','r','Linewidth',2); 
				end
			end
		end
	end
	%annotation('arrow', (0.5+ 0.8*(4*pi*d_a*y_0/2/R_0 / 2 / x_l ))*ones(1,2), [.69 .7],'Color','k')
	%plot(4*pi*d_a*y_0/2/R_0*[1 1],[.5 4.5],'--k', 'Linewidth',2.5)
	plot_para('Maximize',true,'Filename','AmbiGene1')

end