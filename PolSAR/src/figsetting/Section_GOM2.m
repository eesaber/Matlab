function Section_GOM2(I, x, y, y_label,f_name)
    global im_size;
    linew = 2;

    hold on
    plot([1, im_size(2)], y(1)*ones(2,1),'k','Linewidth', linew)
    plot([1, im_size(2)], y(2)*ones(2,1), 'k','Linewidth', linew)
    plot([1, im_size(2)], y(3)*ones(2,1), 'k','Linewidth', linew)
    text(100, y(1)-150,'$L_1$','FontSize',36,'FontName','CMU Serif Roman', 'Interpreter', 'latex')
    text(100, y(2)-150,'$L_2$','FontSize',36,'FontName','CMU Serif Roman', 'Interpreter', 'latex')
    text(100, y(3)-150,'$L_3$','FontSize',36,'FontName','CMU Serif Roman', 'Interpreter', 'latex')
    %{
    plot(x(1)*ones(2,1), [1, im_size(1)], '-k','Linewidth', linew)
    plot(x(2)*ones(2,1), [1, im_size(1)], '--k','Linewidth', linew)
    plot(x(3)*ones(2,1), [1, im_size(1)], '-.k','Linewidth', linew)
    hold off
    %}
    plot_para('Filename',f_name, 'Maximize',true)
    
    font_size = 40;
    movavg = ones(100,1);
    movavg = movavg/numel(movavg);
    %{
    % Along range
    figure
        plot(conv(I(:, x(1)), movavg, 'same'), 'k','Linewidth', linew)
        xlim([1,im_size(1)])
        set(gca,'XTick',1:600:3300, 'Xticklabel',cellstr(int2str((0:600*5:3300*5)'/1000))')
        xlabel('Range (km)')
        hold on 
        plot(conv(I(:,x(2)), movavg, 'same'), 'r','Linewidth', linew)
        plot(conv(I(:,x(3)), movavg, 'same'), 'b','Linewidth', linew)
        hold off
        plot_para('Maximize',true,'Ratio',[4 3 1],'Filename',[f_name '_rm'])
    figure
        plot(I(:,x(1)), '-k', 'Linewidth', linew)
        xlim([1,im_size(1)])
        set(gca,'XTick',1:600:3300, 'Xticklabel',cellstr(int2str((0:600*5:3300*5)'/1000))')
        xlabel('Range (km)')
        hold on 
        plot(I(:,x(2)), '--k', 'Linewidth', linew)
        plot(I(:,x(3)), '-.k','Linewidth', linew)
        hold off
        plot_para('Maximize',true,'Filename',[f_name '_r'])
    %}
    % Along azimuth
    figure
        if 1
            segmentAzimuth(conv(I(y(1),:), movavg, 'same'), 1, 'k', linew)
            set(gca,'XTick',1:1000:8000, 'Xticklabel',cellstr(int2str((0:1000*7:8000*7)'/1000))','FontSize',font_size,'Fontname','CMU Serif Roman','Linewidth',2)
            %plot(conv(I(y(1),:), movavg, 'same'), 'k','Linewidth', linew)
            segmentAzimuth(conv(I(y(2),:), movavg, 'same'), 2, 'r', linew)
            %plot(conv(I(y(2),:), movavg, 'same'), 'r','Linewidth', linew)
            segmentAzimuth(conv(I(y(3),:), movavg, 'same'), 3, 'b', linew)
            %plot(conv(I(y(3),:), movavg, 'same'), 'b','Linewidth', linew)
            grid on
            hold off
            %set(gca, 'Ylim', [0.5, 0.8])
            xlabel('Azimuth (km)')
            ylabel(y_label,'Interpreter','latex')
        else 
            subplot(3,1,1);
            %plot(conv(I(y(1),:), movavg, 'same'), 'k','Linewidth', linew)
            segmentAzimuth(conv(I(y(1),:), movavg, 'same'), 1, 'k', linew)            
            set(gca,'XTick',1:1000:8000, 'Xticklabel',cellstr(int2str((0:1000*7:8000*7)'/1000))','FontSize',font_size,'Fontname','CMU Serif Roman','Linewidth',2)
            %Ylim([0 1])
            subplot(3,1,2);
            %plot(conv(I(y(2),:), movavg, 'same'), 'r','Linewidth', linew)
            segmentAzimuth(conv(I(y(1),:), movavg, 'same'), 1, 'k', linew)            
            set(gca,'XTick',1:1000:8000, 'Xticklabel',cellstr(int2str((0:1000*7:8000*7)'/1000))','FontSize',font_size,'Fontname','CMU Serif Roman','Linewidth',2)
            %Ylim([0 1])
            subplot(3,1,3);
            %plot(conv(I(y(3),:), movavg, 'same'), 'b','Linewidth', linew)
            segmentAzimuth(conv(I(y(1),:), movavg, 'same'), 1, 'k', linew)            
            xlabel('Azimuth (km)')
            ylabel(y_label,'Interpreter','latex')
            set(gca,'XTick',1:1000:8000, 'Xticklabel',cellstr(int2str((0:1000*7:8000*7)'/1000))')
            %Ylim([0 1])
        end
        %plot_para('Maximize',true,'Filename',[f_name '_am'],'Ratio',[4 3 1])
    %{
    figure
        if 1
            plot(I(y(1),:), 'k','Linewidth', linew)
            hold on 
            plot(I(y(2),:), 'r','Linewidth', linew)
            plot(I(y(3),:), 'b','Linewidth', linew)
            hold off
            xlabel('Azimuth (km)')
            ylabel(y_label,'Interpreter','latex')
        else 
            subplot(3,1,1);
            plot(I(y(1),:), 'k','Linewidth', linew)
            set(gca,'XTick',1:1000:8000, 'Xticklabel',cellstr(int2str((0:1000*7:8000*7)'/1000))','FontSize',font_size,'Fontname','CMU Serif Roman','Linewidth',2)
            %Ylim([0 1])
            subplot(3,1,2);
            plot(I(y(2),:), 'r','Linewidth', linew)
            set(gca,'XTick',1:1000:8000, 'Xticklabel',cellstr(int2str((0:1000*7:8000*7)'/1000))','FontSize',font_size,'Fontname','CMU Serif Roman','Linewidth',2)
            %Ylim([0 1])
            subplot(3,1,3);
            plot(I(y(3),:), 'b','Linewidth', linew)
            xlabel('Azimuth (km)')
            ylabel(y_label,'Interpreter','latex')            
            set(gca,'XTick',1:1000:8000, 'Xticklabel',cellstr(int2str((0:1000*7:8000*7)'/1000))')
            %Ylim([0 1])
        end
        plot_para('Maximize',true,'Filename',[f_name '_a'])
        %}
end
function segmentAzimuth(y, pos, linecolor, linewidth)
    linewidth = linewidth+0.5;
    seg = [2672, 5447; 3597, 6300; 4050, 6393]; % corresponding range bins: 3100,1600,500
    plot(1: seg(pos, 1),y(1: seg(pos, 1)), [linecolor ':'], 'Linewidth',linewidth)
    hold on 
    plot(seg(pos, 1): seg(pos, 2), y(seg(pos, 1): seg(pos, 2)), linecolor, 'Linewidth', linewidth)
    plot(seg(pos, 2): numel(y), y(seg(pos, 2): end), [linecolor ':'], 'Linewidth', linewidth)
end