function H_Alpha(T)
% H_ALPHA implements the H/alpha diagram 
% T is matrix with size 3x3xn.

	[~,~,num]= size(T);
	figure
	hold on
	for r = 1 : num
		[U, L] = eig(T(:,:,r));
		P = L/sum(L);
		H = -sum(P.*log(P)/log(3));
		alpha = sum(P'.*acosd(abs(U(1,:))));
		scatter(H,alpha,24,'filled','b')
	end
	n = 100;
	x = linspace(0,1,n);
	P = [ones(size(x));x;x]./repmat(2*x+1,[3,1]);
	H = -sum(P.*log(P)/log(3),1);
	H(isnan(H)) = 0;
	plot(H,acosd(0)*(P(2,:)+P(3,:)),'k','Linewidth',2)
	%
	P = [zeros(size(x)); ones(size(x)); 2*x]./repmat(2*x+1,[3,1]);
	H = -sum(P.*log(P)/log(3),1);
	H(isnan(H)) = 0;
	plot(H(1:n/2),x(1:n/2),'k','Linewidth',2)
	%
	P = [2*x-1; ones(size(x)); ones(size(x))]./repmat(2*x+1,[3,1]);
	H = -sum(P.*log(P)/log(3),1);
	plot(H(n/2+1:end), acosd(0)*2./(2*x(n/2+1:end)+1),'k','Linewidth',2)
    %
	set(gca,'Xlim', [0, 1], 'Ylim', [0, 90], 'XTick', 0:0.2:1,'YTick',0:15:90,'Box','On')
    grid on
	hold off
    xlabel('entropy $H$','Interpreter', 'latex')
	ylabel('$\langle \alpha \rangle^\circ$','Interpreter', 'latex')
    plot_para('Maximize',true,'Filename', 'H_alpha_decomp')
    movefile 'H_alpha_decomp.jpg' output
end