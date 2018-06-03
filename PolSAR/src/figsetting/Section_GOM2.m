function Section_GOM2(I, x, y, f_name)
    global im_size;
    linew = 2;

    hold on
    plot([1, im_size(2)], y(1)*ones(2,1),'-k','Linewidth', linew)
    plot([1, im_size(2)], y(2)*ones(2,1), '--k','Linewidth', linew)
    plot([1, im_size(2)], y(3)*ones(2,1), '-.k','Linewidth', linew)

    plot(x(1)*ones(2,1), [1, im_size(1)], '-k','Linewidth', linew)
    plot(x(2)*ones(2,1), [1, im_size(1)], '--k','Linewidth', linew)
    plot(x(3)*ones(2,1), [1, im_size(1)], '-.k','Linewidth', linew)
    hold off
    plot_para('Filename',f_name, 'Maximize',true)

    lpFilt = designfilt('lowpassfir','PassbandFrequency',0.2, ...
         'StopbandFrequency',0.25,'PassbandRipple',0.5, ...
         'StopbandAttenuation',65,'SampleRate',4,'DesignMethod','kaiserwin');
    %fvtool(lpFilt)
    linew = 0.5;
    font_size = 40;
    movavg = ones(100,1);
    movavg = movavg/numel(movavg);
    % Range
    figure
        plot(filter(lpFilt, I(:, x(1))), 'k','Linewidth', linew)
        xlim([1,im_size(1)])
        set(gca,'XTick',1:600:3300, 'Xticklabel',cellstr(int2str((0:600*5:3300*5)'/1000))')
        xlabel('Range (km)')
        hold on 
        plot(filter(lpFilt, I(:,x(2))), 'r','Linewidth', linew)
        plot(filter(lpFilt, I(:,x(3))), 'b','Linewidth', linew)
        hold off
        plot_para('Maximize',true,'Ratio',[4 3 1],'Filename',[f_name '_rf'])
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
    
    % Azimuth
    figure
        if 1 
            plot(filter(lpFilt, I(y(1),:)), 'k','Linewidth', linew)
            set(gca,'XTick',1:1000:8000, 'Xticklabel',cellstr(int2str((0:1000*7:8000*7)'/1000))')
            hold on 
            plot(filter(lpFilt, I(y(2),:)), 'r','Linewidth', linew)
            plot(filter(lpFilt, I(y(3),:)), 'b','Linewidth', linew)
            hold off
            xlabel('Azimuth (km)')
            %xlim([3000, 5000])
        else
            subplot(3,1,1);
            plot(filter(lpFilt, I(y(1),:)), 'k','Linewidth', linew)
            set(gca,'XTick',1:1000:8000, 'Xticklabel',cellstr(int2str((0:1000*7:8000*7)'/1000))','FontSize',font_size,'Fontname','CMU Serif Roman','Linewidth',2)
            subplot(3,1,2);
            plot(filter(lpFilt, I(y(2),:)), 'r','Linewidth', linew)
            set(gca,'XTick',1:1000:8000, 'Xticklabel',cellstr(int2str((0:1000*7:8000*7)'/1000))','FontSize',font_size,'Fontname','CMU Serif Roman','Linewidth',2)
            subplot(3,1,3);
            plot(filter(lpFilt, I(y(3),:)), 'b','Linewidth', linew)
            xlabel('Azimuth (km)')
            set(gca,'XTick',1:1000:8000, 'Xticklabel',cellstr(int2str((0:1000*7:8000*7)'/1000))')
        end
        plot_para('Maximize',true,'Filename',[f_name '_af'])
    figure
        if 1
            plot(conv(I(y(1),:), movavg, 'same'), 'k','Linewidth', linew)
            hold on 
            plot(conv(I(y(2),:), movavg, 'same'), 'r','Linewidth', linew)
            plot(conv(I(y(3),:), movavg, 'same'), 'b','Linewidth', linew)
            hold off
            xlabel('Azimuth (km)')
        else 
            subplot(3,1,1);
            plot(conv(I(y(1),:), movavg, 'same'), 'k','Linewidth', linew)
            set(gca,'XTick',1:1000:8000, 'Xticklabel',cellstr(int2str((0:1000*7:8000*7)'/1000))','FontSize',font_size,'Fontname','CMU Serif Roman','Linewidth',2)
            %Ylim([0 1])
            subplot(3,1,2);
            plot(conv(I(y(2),:), movavg, 'same'), 'r','Linewidth', linew)
            set(gca,'XTick',1:1000:8000, 'Xticklabel',cellstr(int2str((0:1000*7:8000*7)'/1000))','FontSize',font_size,'Fontname','CMU Serif Roman','Linewidth',2)
            %Ylim([0 1])
            subplot(3,1,3);
            plot(conv(I(y(3),:), movavg, 'same'), 'b','Linewidth', linew)
            xlabel('Azimuth (km)')
            set(gca,'XTick',1:1000:8000, 'Xticklabel',cellstr(int2str((0:1000*7:8000*7)'/1000))')
            %Ylim([0 1])
        end
        plot_para('Maximize',true,'Filename',[f_name '_am'])
    figure
        if 1
            plot(I(y(1),:), 'k','Linewidth', linew)
            hold on 
            plot(I(y(2),:), 'r','Linewidth', linew)
            plot(I(y(3),:), 'b','Linewidth', linew)
            hold off
            xlabel('Azimuth (km)')
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
            set(gca,'XTick',1:1000:8000, 'Xticklabel',cellstr(int2str((0:1000*7:8000*7)'/1000))')
            %Ylim([0 1])
        end
        plot_para('Maximize',true,'Filename',[f_name '_a'])
end