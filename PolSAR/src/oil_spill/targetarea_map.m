function targetarea_map()
    % Function TARGETAREA_MAP is used to plot the target area in the Gulf of Mexico.
    % The target area is enclosed by the rectangular, and rigs is marked by 'x'.
    fsize = 24;
    loca = 2;
    switch loca
        case 1
            chk_pw('output_GulfMexico1')
            latlim = [28, 31];
            lonlim = [-90,-85.5];
            basemap(latlim, lonlim, fsize)
            targetarea_map_1(fsize)
            fname = 'rig_loc';
            fig_setting(fsize)
            plot_para('Filename',fname,'Maximize',true)
        case 2
            chk_pw('output_GulfMexico2/area1')
            latlim = [28, 29.5];
            lonlim = [-90,-88];
            fname = 'rig_loc2';
            %% map
            basemap(latlim, lonlim, fsize)
            targetarea_map_2(fsize)
            fig_setting(fsize)     
            plot_para('Filename',fname,'Maximize',true)
            %% optic iamge projection onto map
            targetarea_map_2_optics(fsize, latlim, lonlim)
            fig_setting(fsize)
            plot_para('Filename',[fname '_opt'],'Maximize',true)
    end
    
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

%% Metohd 
function rig_local()
    %% Rigs location
    rig = [28.866,-88.056; 29.108,-87.944; 28.755,-88.267; 28.481,-88.204; 28.695,-87.931; ...
    28.893,-87.985; 28.509,-88.031; 28.445,-88.277; 28.785,-88.089; 28.784,-88.235; 28.859,-88.044];
    for r = 1 : numel(rig)/2
        scatterm(rig(r,1),rig(r,2), 90, [2 35 209]/255, 'filled')
    end
    scatterm(28.738,-88.366, 270, [255 45 45]/255, 'v','filled') % The explosion BP Deep Horizontal
    scatterm(28.962, -88.696, 135, [255 45 45]/255, '^','filled') % The Bouy 
end

function fig_setting(fsize)
    annotation('arrow',[.7 .7],[.7 .75],'HeadStyle','plain','HeadWidth',5,'LineWidth',2)
    annotation('textarrow',[.7 .7],[.75 .7],'string','N','Fontsize',fsize,'FontName','CMU Sans Serif','HeadStyle','plain','HeadWidth',0.00000000001,'LineWidth',0.5)
    %% Move figure downward
    pos = get(gca, 'Position');
    yoffset = -0.05;
    pos(2) = pos(2) + yoffset;
    set(gca, 'Position', pos) 
end

function targetarea_map_1(fsize)
    cor = [28.36,-88.52; 28.84,-86.06; 29.36,-86.18; 28.88,-88.65; 28.87, -87.36]; %[O;P;Q;R;C]
    dy = 0.1; dx = 0.1;
    rig_local()
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
end


function targetarea_map_2(fsize, c)
    %% 影像範圍 target area 2 
    %           原始資料                      舊範圍 (PQ 沒變)                     新範圍
    %   O = [30.1614, -89.9033]   ;   O' = [29.1616, -88.9576]     ;   O'' [28.9757, -88.7802]
    %   P = [28.3194, -88.1611]   ;   P' = [28.3194, -88.1611]     ;   P'' [28.5758, -88.4019]
    %   Q = [28.4465, -87.9849]   ;   Q' = [28.4465, -87.9849]     ;   Q'' [28.7032, -88.2256]
    %   R = [30.2908, -89.7261]   ;   R' = [29.2897, -88.7810]     ;   R'' [29.1036, -88.6036]
    % 
    % O = [30.1614, -89.9033];  
    % P = [28.3194, -88.1611];
    % Q = [28.4465, -87.9849];
    % R = [30.2908, -89.7261];
    % Opp = 1/36842*(O*13142 + P*23699);
    % Rpp = 1/36842*(R*13142 + Q*23699);
    % Ppp = 1/36842*(O*5143 + P*31698);
    % Qpp = 1/36842*(R*5143 + Q*31698);
    %%
    switch nargin
        case 1
            c = 'k';
    end
    rig_local()
    cor = [28.9757, -88.7802; 28.5758, -88.4019; 28.7032, -88.2256; 29.1036, -88.6036; 28.8396, -88.5028]; % OPQRC
    dy = 0.1*[0.5,  -1,-1, 0, 1];
    dx = 0.1*[-0.5,-1.3, 0, 1,-1];
    %% Truncated Target area
    mark = {'O';'P';'Q';'R';'C'};
    for n = 1 : 5
        scatterm(cor(n,1),cor(n,2), 90,[0 0 .5],'filled')
        textm(cor(n,1)+dy(n),cor(n,2)+dx(n),mark(n),'color',c,'FontSize',fsize,'FontName','CMU Serif Roman', 'Interpreter', 'latex')
    end
    for n = 1 : 4
        plotm([cor(n,1),cor(1+mod(n,4),1)],[cor(n,2),cor(1+mod(n,4),2)],'color',c,'Linewidth',1)
    end
    %% Original Covered Area 
    %{
    cor = [30.161421459,-89.903288309;29.1616, -88.9576;29.2897,-88.7810; 30.290794451, -89.726117031]; %[O';O;R;R']
    for n = 1 : 4
        plotm([cor(n,1),cor(1+mod(n,4),1)],[cor(n,2),cor(1+mod(n,4),2)],'k--','Linewidth',1)
        scatterm(cor(n,1),cor(n,2), 30,[0 0 .5],'filled')
    end
    textm(cor(1,1)+dy(1),cor(1,2)+dx(1),'O$^\prime$','FontSize',fsize,'FontName','CMU Serif Roman', 'Interpreter', 'latex')
    textm(cor(4,1)+dy(4),cor(4,2)+dx(4),'R$^\prime$','FontSize',fsize,'FontName','CMU Serif Roman', 'Interpreter', 'latex')
    %}
end

function targetarea_map_2_optics(fsize, latlim, lonlim)
    global dir
    [A,R]= geotiffread([dir 'optics']);
    %%
    figure
    axesm('MapProjection','mercator','MapLatLimit',latlim,'MapLonLimit', lonlim,...
                'Frame','on','Grid', 'on','MlineLocation',1,'PlineLocation',1,...
                'MeridianLabel','on','ParallelLabel','on','FontName','CMU Serif Roman', 'FontSize', fsize)
    geoshow(A,R)
    tightmap
    axis off
    targetarea_map_2(fsize, [255, 211 6]/255)
end