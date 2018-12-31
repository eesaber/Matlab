function HAlphaDiagram(obj, H, alpha_bar, varargin)
    % HALPHADIAGRAM implements the H/alpha diagram 
    %
    % Syntax: eigenDecomp(H, alpha_bar, NAME, VALUE)
    % Specifies function properties using one or more Name,Value pair arguments. 
    % Name-value pair settings apply to all the lines plotted. 
    %
    % Inputs:
    %    H          - Entropy matrix
    %	 alpha_bar  - Average alpha angle matrix
    % Name-Value Pair Arguments:
    %    ('xlim', xlim)
    %    ('ylim', ylim)
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
    p_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'char'},{}); 
    validationFcn_2_ = @(x) validateattributes(x,{'numeric'},{'numel', 2,'increasing'}); 
    addParameter(p_,'Seglinestyle','w--',validationFcn_1_);
    addParameter(p_,'XBinLimits',[0 1],validationFcn_2_);
    addParameter(p_,'YBinLimits',[0 90],validationFcn_2_);
    parse(p_,varargin{:})
    p_ = p_.Results;
    %H = H(H>=p_.XBinLimits(1));
    %H = H(H<=p_.XBinLimits(2));
    %alpha_bar = alpha_bar(H>=p_.XBinLimits(1));
    %alpha_bar = alpha_bar(H<=p_.XBinLimits(2));
    figure
        histogram2(H(:), alpha_bar(:), ...
            'DisplayStyle','tile', ...
            'YBinLimits', p_.YBinLimits, ...
            'XBinLimits', p_.XBinLimits, ...
            'Normalization','pdf',...
            'NumBins', [450 1350]);
        colormap jet
        c = colorbar;
        set(gca,'clim',[0,0.4],'Color',[0 0 0]+0.7)
        set(gca,'XGrid','off','YGrid','off')
    hold on
    %{
    t=-180:1:180;
    x=0.78+0.15*cosd(t+10);
    y=49.37+21*sind(t);
    plot(x,y, 'r','Linewidth',2)
    x=0.3+0.1*cosd(t-10);
    y=14.6+5*sind(t);
    plot(x,y, 'r','Linewidth',2)
    x=0.11+0.1*cosd(t-10);
    y=8.2+5*sind(t);
    plot(x,y, 'r','Linewidth',2)
    
    
    rectangle('Position',[0.65, 20, 0.3, 50],'EdgeColor','r','Linewidth',2)
    rectangle('Position',[0.03, 3.2, 0.13, 9.6170],'EdgeColor','r','Linewidth',2)
    rectangle('Position',[0.18, 8, 0.3, 14.7400],'EdgeColor','r','Linewidth',2)
    %}
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
    plot([0.5,0.5],[0, 90],p_.Seglinestyle,'Linewidth', 2)
    plot([0.9, 0.9],[0 90],p_.Seglinestyle,'Linewidth', 2)
    plot([0,0.5],[47.5, 47.5],p_.Seglinestyle,'Linewidth', 2)
    plot([0,0.5],[42.5, 42.5],p_.Seglinestyle,'Linewidth', 2)
    plot([0.5,0.9],[50, 50],p_.Seglinestyle,'Linewidth', 2)
    plot([0.5,1],[40, 40],p_.Seglinestyle,'Linewidth', 2)
    plot([0.9, 1], [57.5, 57.5],p_.Seglinestyle,'Linewidth', 2)
    % Image setting
    set(gca,'Xlim', p_.XBinLimits, 'Ylim', p_.YBinLimits, ...
        'XTick', linspace(p_.XBinLimits(1), p_.XBinLimits(2),5), ...
        'YTick',linspace(p_.YBinLimits(1), p_.YBinLimits(2),5),'Box','On')
    hold off
    set(gca,'Color', 0.8*ones(3,1));
    xlabel('$H$','Interpreter', 'latex')
    ylabel('$\langle \alpha \rangle$ (deg)','Interpreter', 'latex')
    %set(gca,'Xticklabel',[], 'YTickLabel',[],'TickLength',[1e-10, 1e-10])
    Plotsetting_gray()
    plot_para('Maximize',true,'Filename',[obj.OUTPUT_PATH, '/H_alpha'])
    
end