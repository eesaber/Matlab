function Pauli_decomp(R, G, B, varargin)
    parse_ = inputParser;
	validationFcn_1_ = @(x) validateattributes(x,{'numeric'},{'nonempty'});
    validationFcn_2_ = @(x) validateattributes(x,{'logical'},{});
    validationFcn_3_ = @(x) validateattributes(x,{'char'},{});
	addParameter(parse_,'Contour',[],validationFcn_1_);
    addParameter(parse_,'Saibu',false,validationFcn_2_);
    addParameter(parse_,'Filename','',validationFcn_3_);
	parse(parse_,varargin{:})

    %%
    up_ = 10; low_ = -50;
    
	Pauli = zeros([size(R), 3]);	
	% |S_vv - S_hh|^2 -> double bounce scattering 
	t_p = 10*log10(R);	
	t_p(t_p < low_) = low_;
	t_p(t_p > up_ ) = up_;
	Pauli(:,:,1) = (t_p-low_)/(up_-low_);	
	% |S_hv|^2 -> volume scattering
	t_p= 10*log10(G);
	t_p(t_p < low_) = low_;
	t_p(t_p > up_ ) = up_;
	Pauli(:,:,2) = (t_p-low_)/(up_-low_);
	% |S_vv + S_hh|^2 -> single scattering
	t_p = 10*log10(B);
	t_p(t_p < low_) = low_;
	t_p(t_p > up_ ) = up_;
	Pauli(:,:,3) = (t_p-low_)/(up_-low_);
    image(Pauli)
	
    %image((cat(3, 10*log10(sqrt(R)), 10*log10(sqrt(G)), 10*log10(sqrt(B)))-low_)/(up_-low_))

    %image(cat(3, 10*log10(sqrt(R)), 10*log10(sqrt(G)), 10*log10(sqrt(B))))
    %image((cat(3, R, G, B)+0)/(0.05))
    if numel(parse_.Results.Contour) ~= 0
        hold on 
        [rr1,rr2] = contour(parse_.Results.Contour,100:100:500,'LineColor','y','Linewidth',1,'ShowText','on');
        clabel(rr1,rr2,'Color','w','Fontsize',18)
        hold off
    end
    if numel(parse_.Results.Filename)
        plot_para('Maximize',true,'Filename', parse_.Results.Filename);
    end
    
    if parse_.Results.Saibu
        %% Plot the dominant channel
        dum = ones(size(R));
        figure
        imagesc(dum.*(R>G).*(R>B))
        colormap gray
        set(gca,'Ydir','normal')
        xlabel('Azimuth (pixel)', 'Fontsize', 40)
        ylabel('Range (pixel)', 'Fontsize', 40)
        if numel(parse_.Results.Filename)
            plot_para('Maximize',true,'Filename',[parse_.Results.Filename, '_r'])
            movefile([parse_.Results.Filename, '_r.jpg'],  'output/')
        end
        figure
        imagesc(dum.*(G>R).*(G>B))
        colormap gray
        set(gca,'Ydir','normal')
        xlabel('Azimuth (pixel)', 'Fontsize', 40)
        ylabel('Range (pixel)', 'Fontsize', 40)
        if numel(parse_.Results.Filename)
            plot_para('Maximize',true,'Filename',[parse_.Results.Filename, '_g'])
            movefile([parse_.Results.Filename, '_g.jpg'],  'output/')
        end
        figure
        imagesc(dum.*(B>R).*(B>G))
        colormap gray
        set(gca,'Ydir','normal')
        xlabel('Azimuth (pixel)', 'Fontsize', 40)
        ylabel('Range (pixel)', 'Fontsize', 40)
        if numel(parse_.Results.Filename)
            plot_para('Maximize',true,'Filename',[parse_.Results.Filename, '_b'])
            movefile([parse_.Results.Filename, '_b.jpg'],  'output/')
        end

    end
end