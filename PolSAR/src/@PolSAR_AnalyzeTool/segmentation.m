function segmentation(obj)
    %%
    kai_1 = kai_1(2:end-1,3:end-2);
    kai_2 = kai_2(2:end-1,3:end-2);
    kai_3 = kai_3(2:end-1,3:end-2);
    %%
    figure
    histogram2(kai_3(:), kai_2(:), 750,...
            'DisplayStyle','tile','YBinLimits',[0 5],'XBinLimits',[-3 3])
    colormap jet
    set(gca,'Color',[0 0 0]+0.7,'XGrid','off','YGrid','off')
    hold on 
    var = [linspace(0,10,100), linspace(10,2000,4000)];
    plot(psi(2,var), psi(1,var), 'w' , -psi(2,var), psi(1,var), 'w--', 'Linewidth',2.5)
    set(gca, 'Xtick', -5:5,'xlim',[-3 3], 'ylim', [0 2])
    xlabel('$\kappa_3$','interpreter','latex')
    ylabel('$\kappa_2$','interpreter','latex')
    plot_para('Maximize',true,'Filename','log_cumulant_diagram','Ratio',[4 3 1])
    hold off
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
    nColors = 3;
    % repeat the clustering 3 times to avoid local minima
    [cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','cityblock', ...
                                        'Replicates',5,'MaxIter',500);
    pixel_labels = uint8(reshape(cluster_idx,nrows,ncols));
    clear ab
    %%
    figure
    %imshow(pixel_labels,[],'Border','tight')
    imagesc(pixel_labels)
    %imagesc(labels)
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
    plot_para('Maximize',true,'Filename','label')
end