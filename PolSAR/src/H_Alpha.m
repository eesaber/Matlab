function H_Alpha(T, varargin)
% H_ALPHA implements the H/alpha diagram 
% T is matrix with size 3x3xn.
    parse_ = inputParser;
	validationFcn_1_ = @(x) validateattributes(x,{'char'    },{}); 
	addParameter(parse_,'Pointcolor','b',validationFcn_1_);
	parse(parse_,varargin{:})
    
    
	[~,~,num]= size(T);
	hold on
	for r = 1 : num
		[U, L] = eig(T(:,:,r));
		P = L/sum(L);
		H = -sum(P.*log(P)/log(3));
		alpha = sum(P'.*acosd(abs(U(1,:))));
		scatter(H,alpha,24,'filled',parse_.Results.Pointcolor)
	end
	n = 100;
	x = linspace(0,1,n);
	P = [ones(size(x));x;x]./repmat(2*x+1,[3,1]);
	H = -sum(P.*log(P)/log(3),1);
	H(isnan(H)) = 0;
	plot(H,acosd(0)*(P(2,:)+P(3,:)),'k','Linewidth',2)
	% Plot the boundary
	P = [zeros(size(x)); ones(size(x)); 2*x]./repmat(2*x+1,[3,1]);
	H = -sum(P.*log(P)/log(3),1);
	H(isnan(H)) = 0;
	plot(H(1:n/2),x(1:n/2),'k','Linewidth',2)
	%
	P = [2*x-1; ones(size(x)); ones(size(x))]./repmat(2*x+1,[3,1]);
	H = -sum(P.*log(P)/log(3),1);
	plot(H(n/2+1:end), acosd(0)*2./(2*x(n/2+1:end)+1),'k','Linewidth',2)
    % Plot segmetation lines
	plot([0.5,0.5],[0, 90],'--','Linewidth', 2)
	plot([0.9, 0.9],[0 90],'--','Linewidth', 2)
	plot([0,0.5],[47.5, 47.5],'--','Linewidth', 2)
	plot([0,0.5],[42.5, 42.5],'--','Linewidth', 2)
	plot([0.5,0.9],[50, 50],'--','Linewidth', 2)
	plot([0.5,1],[40, 40],'--','Linewidth', 2)
	plot([0.9, 1], [57.5, 57.5],'--','Linewidth', 2)
	% Image setting
	set(gca,'Xlim', [0, 1], 'Ylim', [0, 90], 'XTick', 0:0.2:1,'YTick',0:15:90,'Box','On')
	hold off
    xlabel('entropy $H$','Interpreter', 'latex')
	ylabel('$\langle \alpha \rangle^\circ$','Interpreter', 'latex')
    plot_para('Maximize',true,'Filename', 'H_alpha_decomp')
    movefile 'H_alpha_decomp.jpg' output
end

%% Plot Each Region onto H/alpha diagram
% bare soil
fan = [343,461;214,368]; % [row, row; col, col]
size_f = (fan(1,2)-fan(1,1)+1)*(fan(2,2)-fan(2,1)+1);
temp = cat(1,cat(2, reshape(T_11(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f]), reshape(T_12(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f]), reshape(T_13(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f])), ......
             cat(2,reshape(conj(T_12(fan(1,1):fan(1,2), fan(2,1): fan(2,2))),[1,1,size_f]), reshape(T_22(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f]), reshape(T_23(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f])),.......
             cat(2,reshape(conj(T_13(fan(1,1):fan(1,2), fan(2,1): fan(2,2))),[1,1,size_f]), reshape(conj(T_23(fan(1,1):fan(1,2), fan(2,1): fan(2,2))),[1,1,size_f]), reshape(T_33(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f])));
ind_ = randperm(size_f);
figure
H_Alpha(temp(:,:,ind_(1:500)),'Pointcolor','b')
% building 
fan = [176,273; 1,212]; % [row, row; col, col]
size_f = (fan(1,2)-fan(1,1)+1)*(fan(2,2)-fan(2,1)+1);
temp = cat(1,cat(2, reshape(T_11(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f]), reshape(T_12(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f]), reshape(T_13(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f])), ......
             cat(2,reshape(conj(T_12(fan(1,1):fan(1,2), fan(2,1): fan(2,2))),[1,1,size_f]), reshape(T_22(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f]), reshape(T_23(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f])),.......
             cat(2,reshape(conj(T_13(fan(1,1):fan(1,2), fan(2,1): fan(2,2))),[1,1,size_f]), reshape(conj(T_23(fan(1,1):fan(1,2), fan(2,1): fan(2,2))),[1,1,size_f]), reshape(T_33(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f])));
ind_ = randperm(size_f);
H_Alpha(temp(:,:,ind_(1:500)),'Pointcolor','r')
% canopy
fan = [214,285; 661,801]; % [row, row; col, col]
size_f = (fan(1,2)-fan(1,1)+1)*(fan(2,2)-fan(2,1)+1);
temp = cat(1,cat(2, reshape(T_11(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f]), reshape(T_12(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f]), reshape(T_13(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f])), ......
             cat(2,reshape(conj(T_12(fan(1,1):fan(1,2), fan(2,1): fan(2,2))),[1,1,size_f]), reshape(T_22(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f]), reshape(T_23(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f])),.......
             cat(2,reshape(conj(T_13(fan(1,1):fan(1,2), fan(2,1): fan(2,2))),[1,1,size_f]), reshape(conj(T_23(fan(1,1):fan(1,2), fan(2,1): fan(2,2))),[1,1,size_f]), reshape(T_33(fan(1,1):fan(1,2), fan(2,1): fan(2,2)),[1,1,size_f])));
ind_ = randperm(size_f);
H_Alpha(temp(:,:,ind_(1:500)),'Pointcolor','k')
