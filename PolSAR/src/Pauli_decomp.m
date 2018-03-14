%function Pauli_decomp(R, G, B, name, saibu, dem)
function Pauli_decomp(R, G, B, name, saibu)
     chk_pw()
    %%
    up_ = 10; low_ = -20;
	Pauli = zeros([size(R), 3]);	
	% |S_vv - S_hh|^2 -> double bounce scattering 
	t_p = 10*log10(sqrt(R));	
	t_p(t_p < low_) = low_;
	t_p(t_p > up_ ) = up_;
	Pauli(:,:,1) = (t_p-low_)/(up_-low_);	
	% |S_hv|^2 -> volume scattering
	t_p= 10*log10(sqrt(G));
	t_p(t_p < low_) = low_;
	t_p(t_p > up_ ) = up_;
	Pauli(:,:,2) = (t_p-low_)/(up_-low_);
	% |S_vv + S_hh|^2 -> single scattering
	t_p = 10*log10(sqrt(B));
	t_p(t_p < low_) = low_;
	t_p(t_p > up_ ) = up_;
	Pauli(:,:,3) = (t_p-low_)/(up_-low_);
	figure
    %hold on 
    img = image(Pauli);
    %[rr1,rr2] = contour(dem(2300:2700,7000:7900),[35,40,50,60,70],'LineColor','b','Linewidth',2,'ShowText','on');
    %clabel(rr1,rr2,'Color','b','Fontsize',18)
    %hold off
    set(gca,'Ydir','normal')
    xlabel('azimuth (pixel)', 'Fontsize', 40)
    ylabel('range (pixel)', 'Fontsize', 40)
    plot_para('Maximize',true,'Filename',name, 'Ratio', [4 3 1]);
    movefile([name, '.jpg'],  'output/')
    
    
    if saibu
        %% Plot the dominant channel
        dum = ones(size(R));
        figure
        imagesc(dum.*(R>G).*(R>B))
        colormap gray
        set(gca,'Ydir','normal')
        xlabel('azimuth (pixel)', 'Fontsize', 40)
        ylabel('range (pixel)', 'Fontsize', 40)
        plot_para('Maximize',true,'Filename',[name, '_r'])
        movefile([name, '_r.jpg'],  'output/')
        figure
        imagesc(dum.*(G>R).*(G>B))
        colormap gray
        set(gca,'Ydir','normal')
        xlabel('azimuth (pixel)', 'Fontsize', 40)
        ylabel('range (pixel)', 'Fontsize', 40)
        plot_para('Maximize',true,'Filename',[name, '_g'])
        movefile([name, '_g.jpg'],  'output/')
        figure
        imagesc(dum.*(B>R).*(B>G))
        colormap gray
        set(gca,'Ydir','normal')
        xlabel('azimuth (pixel)', 'Fontsize', 40)
        ylabel('range (pixel)', 'Fontsize', 40)
        plot_para('Maximize',true,'Filename',[name, '_b'])
        movefile([name, '_b.jpg'],  'output/')

    end
    
end