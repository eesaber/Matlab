function HAlphaDiagram(obj, H, alpha_bar, varargin)
	% HALPHADIAGRAM implements the H/alpha diagram 
    %
    % Syntax: eigenDecomp(H, alpha_bar, NAME, VALUE)
    % Specifies function properties using one or more Name,Value pair arguments. 
    % Name-value pair settings apply to all the lines plotted. 
    %
    % Inputs:
	%    H			- Entropy matrix
	%	 alpha_bar  - Average alpha angle matrix
	% Name-Value Pair Arguments:
	%    ('Seglinestyle', s) - Set the segmentation line style. s is a char array.
	%    ('BackgroundColor', c) - Set the background color. c is a char or 
    % Example: 
    %    1. HALPHADIAGRAM(H, alpha_bar)
    %       
    %    2. HALPHADIAGRAM(H, alpha_bar, 'Seglinestyle', 'b', )
    %
    % Other m-files required: plot_para
    % Subfunctions: none
    % MAT-files required: none
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------
    parse_ = inputParser;
	validationFcn_1_ = @(x) validateattributes(x,{'char'},{}); 
	addParameter(parse_,'Seglinestyle','w--',validationFcn_1_);
	parse(parse_,varargin{:})
	
	figure
		histogram2(H(:), alpha_bar(:),'DisplayStyle','tile',...
			'YBinLimits',[0 90],'XBinLimits',[0 1],'NumBins',[450 1350]);
		colormap jet 
		set(gca,'Color',[0 0 0]+0.7,'XGrid','off','YGrid','off')
    hold on
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
	plot([0.5,0.5],[0, 90],parse_.Results.Seglinestyle,'Linewidth', 2)
	plot([0.9, 0.9],[0 90],parse_.Results.Seglinestyle,'Linewidth', 2)
	plot([0,0.5],[47.5, 47.5],parse_.Results.Seglinestyle,'Linewidth', 2)
	plot([0,0.5],[42.5, 42.5],parse_.Results.Seglinestyle,'Linewidth', 2)
	plot([0.5,0.9],[50, 50],parse_.Results.Seglinestyle,'Linewidth', 2)
	plot([0.5,1],[40, 40],parse_.Results.Seglinestyle,'Linewidth', 2)
	plot([0.9, 1], [57.5, 57.5],parse_.Results.Seglinestyle,'Linewidth', 2)
	% Image setting
	set(gca,'Xlim', [0, 1], 'Ylim', [0, 90], 'XTick', 0:0.2:1,'YTick',0:15:90,'Box','On')
	hold off
	set(gca,'Color', 0.8*ones(3,1));
    xlabel('$H$','Interpreter', 'latex')
	ylabel('$\langle \alpha \rangle$ $(^\circ)$','Interpreter', 'latex')
    plot_para('Maximize',true,'Filename',[obj.OUTPUT_PATH, '/H_alpha'])
    
end