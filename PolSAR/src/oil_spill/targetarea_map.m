function targetarea_map()
    % Function TARGETAREA_MAP is used to plot the target area in the Gulf of Mexico.
    % The target area is enclosed by the rectangular, and rigs is marked by 'x'.
    oclayer = wmsfind('Oceans/Seas','SearchFields','layertitle');
    latlim = [28, 31];
    lonlim = [-90,-85.5];
    imageLength = 3000;
    fsize = 24;
    
    figure
    [A, R] = wmsread(oclayer(2),'Latlim', latlim, 'Lonlim', lonlim, ...
        'ImageHeight', imageLength, 'ImageWidth', imageLength);
    axesm('MapProjection','mercator','MapLatLimit',latlim,'MapLonLimit', lonlim,...
            'Frame','on','Grid', 'on','MlineLocation',1,'PlineLocation',1,...
            'MeridianLabel','on','ParallelLabel','on','FontName','CMU Serif Roman', 'FontSize', fsize)
    geoshow(A,R)
    tightmap
    axis off
    annotation('arrow',[.5 .5],[.6 .65],'HeadStyle','plain','HeadWidth',5,'LineWidth',2)
    annotation('textarrow',[.5 .5],[.65 .6],'string','N','Fontsize',fsize,'FontName','CMU Sans Serif','HeadStyle','plain','HeadWidth',0.00000000001,'LineWidth',0.5)
    %% Rigs location
    rig = [28.866,-88.056; 29.108,-87.944; 28.755,-88.267; 28.481,-88.204; 28.695,-87.931; ...
    28.893,-87.985; 28.509,-88.031; 28.445,-88.277; 28.785,-88.089; 28.784,-88.235; 28.859,-88.044];
    for r = 1 : numel(rig)/2
        scatterm(rig(r,1),rig(r,2), 30, [2 35 209]/255, 'filled')
    end
    scatterm(28.738,-88.366, 90, [255 45 45]/255, 'x') % The explosion BP Deep Horizontal
    scatterm(28.962, -88.696, 45, [255 45 45]/255, '^','filled') % The Bouy 
    %% Move figure downward
    pos = get(gca, 'Position');
    yoffset = -0.05;
    pos(2) = pos(2) + yoffset;
    set(gca, 'Position', pos)
    targetarea_map_2(fsize)
end
function targetarea_map_1(fsize)
    
    dy = 0.1; dx = 0.1;
    %% Target area
    scatterm(cor(1,1),cor(1,2), 30,[0 0 .5],'filled'); textm(cor(1,1)-dy,cor(1,2),'O','FontSize',fsize,'FontName','CMU Serif Roman'); % O
    scatterm(cor(2,1),cor(2,2), 30,[0 0 .5],'filled'); textm(cor(2,1)-dy,cor(2,2)+dx,'P','FontSize',fsize,'FontName','CMU Serif Roman') % P
    scatterm(cor(3,1),cor(3,2), 30,[0 0 .5],'filled'); textm(cor(3,1)+dy,cor(3,2)+dx,'Q','FontSize',fsize,'FontName','CMU Serif Roman') % Q
    scatterm(cor(4,1),cor(4,2), 30,[0 0 .5],'filled'); textm(cor(4,1)+dy,cor(4,2)-dx,'R','FontSize',fsize,'FontName','CMU Serif Roman') % R
    scatterm(cor(5,1),cor(5,2), 30,[0 0 .5],'filled'); textm(cor(5,1)+dy,cor(5,2)+dx,'C','FontSize',fsize,'FontName','CMU Serif Roman') % C
    plotm([28.88 29.36],[-88.65 -86.18],'k','LineWidth',1) % R-Q
    plotm([28.36 28.88],[-88.52 -88.65],'k','Linewidth',1) % R-O
    plotm([28.84 28.36],[-86.06 -88.52],'k','Linewidth',1) % O-P
    plotm([29.36 28.84],[-86.18 -86.06],'k','Linewidth',1) % P-Q
    plot_para('Filename','output/rig_local','Maximize',true)
end

function targetarea_map_2(fsize)
    cor = [29.1616, -88.9576; 28.319386246,-88.161090341; 28.446506933,-87.984930547; 29.2897,-88.7810; 28.8043,-88.4712]; % OPQRC
    dy = 0.1*[-1.2,  -1,-1, 0, 1];
    dx = 0.1*[-0.8,-1.3, 0, 1,-1];
    %% target area
    mark = {'O';'P';'Q';'R';'C'};
    for n = 1 : 5
        scatterm(cor(n,1),cor(n,2), 30,[0 0 .5],'filled')
        textm(cor(n,1)+dy(n),cor(n,2)+dx(n),mark(n),'FontSize',fsize,'FontName','CMU Serif Roman', 'Interpreter', 'latex')
    end
    for n = 1 : 4
        plotm([cor(n,1),cor(1+mod(n,4),1)],[cor(n,2),cor(1+mod(n,4),2)],'k','Linewidth',1)
    end
    cor = [30.161421459,-89.903288309;29.1616, -88.9576;29.2897,-88.7810; 30.290794451, -89.726117031]; %[O';O;R;R']
    for n = 1 : 4
        plotm([cor(n,1),cor(1+mod(n,4),1)],[cor(n,2),cor(1+mod(n,4),2)],'k--','Linewidth',1)
        scatterm(cor(n,1),cor(n,2), 30,[0 0 .5],'filled')
    end
    textm(cor(1,1)+dy(1),cor(1,2)+dx(1),'O$^\prime$','FontSize',fsize,'FontName','CMU Serif Roman', 'Interpreter', 'latex')
    textm(cor(4,1)+dy(4),cor(4,2)+dx(4),'R$^\prime$','FontSize',fsize,'FontName','CMU Serif Roman', 'Interpreter', 'latex')
    
    plot_para('Filename','rig_local2','Maximize',true)
end