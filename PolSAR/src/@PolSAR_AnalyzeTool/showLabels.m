function showLabels(obj, labels, nColors)
    if size(labels) ~= obj.IMAGE_SIZE
        error('The dimention of image and labels are different.')
    end
    figure
        imagesc(labels)
        caxis([1 nColors])
        if nColors <= 4
            cmap  = flag(nColors);
        else
            cmap  = pink(nColors);
        end
        cmap = linspecer(nColors,'qualitative');
        colormap(cmap)
        L = line(ones(nColors),ones(nColors), 'LineWidth',2);               % generate line 
        set(L,{'color'},mat2cell(cmap, ones(1,nColors),3)); % set the colors according to cmap
        legend(num2str((1:nColors)'),'Location','bestoutside')
        plot_para('Maximize',true,'Ratio',[4 3 1],'Filename',[obj.OUTPUT_PATH '/label'])
    
    %% Histogram of sigma_vv
    % bin_limit = [0, 0.1];
    bin_limit = [0, 0.1];
    figure
        hold on 
        for n = 1:nColors
            % 4e-4   
            % 1e-3
            histogram(obj.vv_vv(labels==n),'BinLimits', bin_limit, 'BinWidth',4e-4,...
                'LineWidth',1,'Normalization','probability', 'FaceColor',cmap(n,:))
        end
        hold off
        legend(strcat('label', num2str((1:nColors)')))
        ylim([0 0.06])
        set(gca, 'Box','on', 'Xlim', bin_limit, 'Ygrid', 'on', 'Xgrid','on')
        xlabel('$\sigma_{vv}$', 'interpreter', 'latex')
        ylabel('Pr$(x = \sigma_{vv})$', 'interpreter', 'latex')
        plot_para('Maximize',true,'Ratio',[4 3 1], 'FileName',[obj.OUTPUT_PATH '/histo_vv'])
end