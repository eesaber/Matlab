function targetarea_map_2()
    % Function TARGETAREA_MAP is used to plot the target area in the North Sea.
    % The target area is enclosed by the rectangular, and rigs is marked by 'x'.
    fsize = 24;

    chk_pw('output_NorthSea2/')
    latlim = [59.5, 60.5];
    lonlim = [1.5 3];
    fname = 'rig_local_vad';
    %% map
    basemap(latlim, lonlim, fsize)
    fig_setting(fsize)
    targetarea_map(fsize)
    plot_para('Filename',fname,'Maximize',true)

end

function basemap(latlim, lonlim, fsize)
    oclayer = wmsfind('Oceans/Seas','SearchFields','layertitle');
    imageLength = 3000;
    figure
    [A, R] = wmsread(oclayer(2),'Latlim', latlim, 'Lonlim', lonlim, ...
        'ImageHeight', imageLength, 'ImageWidth', imageLength);
    axesm('MapProjection','mercator','MapLatLimit',latlim,'MapLonLimit', lonlim,...
            'Frame','on','Grid', 'on','MlineLocation',1,'PlineLocation',1,...
            'MeridianLabel','on','ParallelLabel','on','FontName','CMU Serif Roman', 'FontSize', fsize)
    geoshow(A,R)
    tightmap
    axis off
end

function fig_setting(fsize)
    annotation('arrow',[.6 .6],[.7 .75],'HeadStyle','plain','HeadWidth',5,'LineWidth',2)
    annotation('textarrow',[.6 .6],[.75 .7],'string','N','Fontsize',fsize,'FontName','CMU Sans Serif','HeadStyle','plain','HeadWidth',0.00000000001,'LineWidth',0.5)
    %% Move figure downward
    pos = get(gca, 'Position');
    yoffset = -0.05;
    pos(2) = pos(2) + yoffset;
    set(gca, 'Position', pos) 
end

function targetarea_map(fsize, c)
    switch nargin
        case 1
            c = 'k';
    end
    
    cor = [60.19714872, 2.17886710; 59.74766832, 2.07531258;......
    59.73444504, 2.46197538; 60.18514776, 2.57097429; 59.94, 23.2]; % OPQRC
    dy = 0.1*[0.5,  -1,-1, 0, 1];
    dx = 0.1*[-0.5,-1.3, 0, 1,-1];
    %% Truncated Target area
    mark = {'O';'P';'Q';'R';'C'};
    cor_m = [28.9757, -88.8472; 28.5558, -88.3719; 28.7032, -88.2056; 29.1036, -88.5936; 28.8496, -88.5598];
    for n = 1 : 5
        scatterm(cor(n,1),cor(n,2), 90,[0 0 .5],'filled')
        textm(cor_m(n,1),cor_m(n,2),mark(n),'color',c,'FontSize',fsize,'FontName','CMU Serif Roman', 'Interpreter', 'latex')
    end
    for n = 1 : 4
        plotm([cor(n,1),cor(1+mod(n,4),1)],[cor(n,2),cor(1+mod(n,4),2)],'color',c,'Linewidth',1)
    end

end
