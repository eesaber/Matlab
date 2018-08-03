function map_beaufort()    
    figure('Color','w')
    set(gca,'Visible','off')
    %{
    axesm('eqaazim','MapLatLimit',[65 90],'FLatLimit',[-180 -60],......
        'Frame','on','Grid', 'on','MlineLocation',30,'PlineLocation',10,......
        'GLineStyle','.-','MeridianLabel','on','ParallelLabel','on','FontName','CMU Serif Roman','Fontsize',32)
    setm(gca,'MLabelParallel',0,'ParallelLabel','on')
    %}
    % eqdconic flatplrq
    axesm('eqdcylin','MapLatLimit',[73 77],'MapLonLimit',[-159 -156],...
    'Frame','off','Grid', 'on','MlineLocation',1,'PlineLocation',1,...
    'GLineStyle','.-','MeridianLabel','on','ParallelLabel','on',...
    'PlabelMeridian','east','MlabelParallel','south',...
    'LabelFormat','signed','LabelRotation','on',...
    'Fontsize',32)
    %daspectm('m',1)
    %view(-33,18)
    %{
    % eqaazim
    axesm('eqdcylin','origin',[60 30],...%'FLatLimit',[-20 20],...
       'geoid',wgs84Ellipsoid,...
        'Frame','off','Grid', 'on','MlineLocation',5,'PlineLocation',5,...
    'GLineStyle','.-','MeridianLabel','on','ParallelLabel','on',...
    'PlabelMeridian','east','MlabelParallel','south',...
    'LabelFormat','signed','LabelRotation','off',...
    'Fontsize',32)
    %}
    geoshow('landareas.shp', 'FaceColor', [255,218,150]/255)
    geoshow('worldrivers.shp','Color', 'blue')
    tightmap
    hold on 
    roi()
    plane_trajectory()
    hold off
    %view(-90, 90);
    %plot_para('Maximize',true,'Filename','map','Ratio',[5,1 ,1])
    
end
function plane_trajectory()
    lat = double(ncread('/home/akb/Code/Matlab/RDWES1B_CS_OFFL_SIR_SAR_1B_20151005T123541_20151005T124356_C001.nc','lat'));
    lon = double(ncread('/home/akb/Code/Matlab/RDWES1B_CS_OFFL_SIR_SAR_1B_20151005T123541_20151005T124356_C001.nc','lon'));
    elev = double(ncread('/home/akb/Code/Matlab/RDWES1B_CS_OFFL_SIR_SAR_1B_20151005T123541_20151005T124356_C001.nc','elev'));
    plotm(lat, lon, 'LineWidth',3)
    lat = lat(1:3000);
    lon = lon(1:3000);
    elev = elev(1:3000);
    l = 0.01;
    surfacem([lat(:) l+lat(:)], [lon(:) l+lon(:)], [elev(:) elev(:)],...
        'CData',[Stat(elev(:)).' Stat(elev(:)).'],...
        'AlphaData',0.1*ones(size([lat lat])),...
        'AlphaDataMapping','none')
    colormap jet
end
function roi()
    info = geotiffinfo('/media/akb/2026EF9426EF696C/raw_data/Beaufort/beaufo_01105_15148_004_151006_L090_CX_01_hgt.tif');
    fid = fopen('/media/akb/2026EF9426EF696C/raw_data/Beaufort/beaufo_01105_15148_004_151006_L090VVVV_CX_01.grd','r','ieee-le'); 
    vv_vv = (single(fread(fid, [17186 53672],'real*4'))).';
    t = imresize(10*log10(vv_vv), 0.1);
    t = reshape(Stat(t), size(t));
    t(t==0) = 255;
    info.SpatialRef.RasterSize = size(t);
    geoshow(t,info.SpatialRef,'DisplayType','Image')
    colormap gray
end