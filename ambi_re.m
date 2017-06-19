function ambi_re(qq, lamda_)
	% Note that the frequency is f_a > f_b > f_c 
	fprintf('Note that the frequency is (f_a,f_b,f_c) f_a > f_b > f_c \n')
	load('para.mat')
	N = 20;
	T = 20;
	N_down = min(floor(4*pi./lamda_*d_a/2/R_0*(-T*(v_p - v_y) /2 + y_0) / 2 /pi));
	N_up = max(floor(4*pi./lamda_*d_a/2/R_0*(T*(v_p - v_y) /2 + y_0) / 2 /pi));
	space = 2*pi*(N_down: N_up);
	k = repmat(qq.',1,length(space)) + lamda_.' * space ;
	k_bar = [k(:,1),  diff(k,1,2)];
	%% Plot the 
	x_l = 1;
	figure 
	subplot('position',[0.1 0.5 0.8 0.16]);
	%subplot('position',[0.1 0.3375 0.8 0.1125]);
	barh(k_bar,1,'stack','FaceColor','none','Linewidth',1.25,'EdgeColor',[204, 158, 86]/256	...
		,'BaseValue', -100, 'ShowBaseLine', 'on');
	set(gca,'TickLength',[0 0])
	ylim([0.5 numel(k)/length(k)+0.5])

	%yticklabels({'\lambda_1 ', '\lambda_2','\lambda_3'})
	yticks([0.75 2 3 4.2])
	yticklabels({'\lambda_1 ', '\lambda_2', '\lambda_3', '\lambda_4'})
	hold all
	x_l = 2*x_l*(N>0) - x_l;
	xlim([-x_l,x_l])
	xlabel(['$\lambda_k ( \beta_k \hat{y}_k + 2 N \pi)$, ', int2str(N_down), '$ \leq N \leq $', int2str(N_up)],'Interpreter', 'latex')

	[temp, ind_] = sort(reshape(k,1,numel(k)));
	row_ = reshape( (1:numel(k)/length(k)).' *ones(1,length(space)), 1, numel(k)); 
	row_ = row_(ind_);
	for i = 1 : numel(k)- numel(k)/length(k)
		if temp(i+numel(k)/length(k)-1) - temp(i) < 0.01 && sum(diff(sort(row_(i: i + numel(k)/length(k) - 1))) == 0) == 0
			for j = 0 : 3
				%animatedline([temp(i+j) temp(i+j)],[row_(i+j)-0.5 row_(i+j)+0.5],'Color',[ 2,65,226]/256,'Linewidth',2); 
				%drawnow
				plot([temp(i+j) temp(i+j)],[row_(i+j)-0.5 row_(i+j)+0.5],'Color',[ 2,65,226]/256,'Linewidth',2); 
			end
		end
	end

	annotation('arrow', (0.5+ 0.8*(4*pi*d_a*y_0/2/R_0 / 2 / x_l ))*ones(1,2), [.64 .66],'Color','k')
	plot(4*pi*d_a*y_0/2/R_0*[1 1],[.5 4.5],'Color', 'k', 'Linewidth',2)
	plot_para('Maximize',true,'Filename','AmbiGene1')

end