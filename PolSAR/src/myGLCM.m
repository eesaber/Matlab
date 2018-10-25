function texture = myGLCM(im, f_row, f_col)
    % MYGLCM generate the textrue feature
    %
    % Syntax:
	%  *
	%
    % Inputs:
    %  *
    %
    % Outputs:
    %  *
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------

%%
    FILTER_SIZE = [f_row, f_col];
    FILTER_ROW_SIZE = FILTER_SIZE(1);
    FILTER_COL_SIZE = FILTER_SIZE(2); 
    IM_ROW_SIZE = size(im,1);
    IM_COL_SIZE = size(im,2);
    texture_meanfirst = zeros([size(im)-FILTER_SIZE-1, 4]);

    parfor it_r = 1 : IM_ROW_SIZE-FILTER_ROW_SIZE-1
        offset = [0 1; -1 1; -1 0; -1 -1];
        %offset = [0, 1; -1, 1; -1, 0; -1, -1; 0, -1; 1, -1; 1, 0; 1, 1];
        prop = {'Contrast','Correlation','Energy','Homogeneity'};
        for it_c = 1 : IM_COL_SIZE-FILTER_COL_SIZE-1
            glcm = graycomatrix(im(it_r:it_r+FILTER_ROW_SIZE, it_c:it_c+FILTER_COL_SIZE), ...
            'NumLevels',64, 'Offset', offset);

            glcm = floor(mean(glcm, 3));
            stats = graycoprops(glcm, prop);
            texture_meanfirst(it_r, it_c,:) = cat(3, stats.Contrast, ...
                stats.Correlation, stats.Energy, stats.Homogeneity);
            
        end
    end

    %% Pad the boarder as 
    f = @(x) padarray(x,floor(FILTER_SIZE/2),'replicate','both');
    texture = f(texture_meanfirst);
    %%
    prop = {'Contrast','Correlation','Energy','Homogeneity'};
    figure
    for it = 1 : 4
        subplot(2,2,it)
        imagesc(texture(:,:,it))
        title(prop(it))
        set(gca,'Ydir','normal')
    end
    a = axes;
    t = title(['filter size: ' num2str(FILTER_SIZE)], 'Fontsize',24);
    a.Visible = 'off';
    t.Visible = 'on';
    set(gcf, 'Position', get(0, 'Screensize'));
    drawnow
end