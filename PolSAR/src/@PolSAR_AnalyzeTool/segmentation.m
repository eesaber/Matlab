function [labels] = segmentation(obj, kai_1, kai_2, kai_3, nColors)
    % SEGMENTATION Implement the classification by using k_means 
    %
    % Syntax: SEGMENTATION(kai_1,kai_2,kai_3)
    %
    % Inputs:
	%    kai_1  - first order log-cumulant. The definition of log-cumulant is in [1].
	%    kai_2  - second order log-cumulant.
    %    kai_3  - third order log-cumulant.
    %
    % Other m-files required: plot_para, linspecer
    % Subfunctions: kappa23Diagram, 
    % MAT-files required: none
    %
    % Reference 
    % [1] Introduction to second kind statistics Application of log-moments and 
    %     log-cumulants to the analysis of radar image distributions
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------

    %% Crop the image for preventing the boarder's value 
    kai_1 = kai_1(2:end-1,3:end-2);
    kai_2 = kai_2(2:end-1,3:end-2);
    kai_3 = kai_3(2:end-1,3:end-2);

    %% (\kappa_3, \kappa_2) diagram
    obj.logCumulantDiagram(kai_2, kai_3)
    
    %% Segementation
    temp = cat(3, reshape(Stat(kai_1,'uint16'), size(kai_1)),...
        reshape(Stat(kai_2, 'uint16'), size(kai_1)), ...
        reshape(Stat(kai_3, 'uint16'), size(kai_1)));
    %temp = cat(3, kai_1, kai_2, kai_3);
    lab_I = rgb2lab(temp);
    clear temp
    ab = lab_I(:,:,2:3);
    nrows = size(ab,1);
    ncols = size(ab,2);
    ab = reshape(ab,size(ab,1)* size(ab,2),2);
    % repeat the clustering 3 times to avoid local minima
    [cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','cityblock', ...
                                        'Replicates',5,'MaxIter',500);
    pixel_labels = uint8(reshape(cluster_idx,nrows,ncols));
    clear ab
    figure
    %imshow(pixel_labels,[],'Border','tight')
    imagesc(pixel_labels)
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
    plot_para('Maximize',true,'Filename',[obj.OUTPUT_PATH '/label'])

    %% Analysis
    labels = padarray(pixel_labels, [1,2], -1,'both');
    distribution(obj, nColors, cmap, labels)
end

function distribution(obj, nColors, cmap, labels)
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
    ylim([0 0.03])
    set(gca, 'Box','on', 'Xlim', bin_limit, 'Ygrid', 'on', 'Xgrid','on')
    xlabel('$\sigma_{vv}$', 'interpreter', 'latex')
    ylabel('Pr$(x = \sigma_{vv})$', 'interpreter', 'latex')
    plot_para('Maximize',true,'Ratio',[4 3 1],'FileName',[obj.OUTPUT_PATH '/histo_vv'])
end

function scatterPlot(kai_1, kai_2, kai_3, labels)
    kai_1 = padarray(kai_1, [1,2], -1,'both');
    kai_2 = padarray(kai_2, [1,2], -1,'both');
    kai_3 = padarray(kai_3, [1,2], -1,'both');
    %% 2D 
    figure
    hold on 
    for n = 1 : nColors
        scatter(kai_3(labels==n), kai_2(labels==n),'filled')
    end
    hold off
    set(gca, 'Xtick', -5:5,'xlim',[-3 3], 'ylim', [0 2])
    xlabel('$\kappa_3$','interpreter','latex')
    ylabel('$\kappa_2$','interpreter','latex')
    %% 3D
    figure
    hold on 
    for n = 1 : nColors
        scatter3(kai_3(labels==n), kai_2(labels==n), kai_1(labels==n), 'filled')
    end
    hold off
    %set(gca, 'Xtick', -5:5,'xlim',[-3 3], 'ylim', [0 2])
    xlabel('$\kappa_3$','interpreter','latex')
    ylabel('$\kappa_2$','interpreter','latex')
    ylabel('$\kappa_1$','interpreter','latex')
end