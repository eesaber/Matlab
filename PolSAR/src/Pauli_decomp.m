function Pauli_decomp(R, G, B, name)
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
    image(Pauli)
    set(gca,'Ydir','normal')
    xlabel('Azimuth (pixel)', 'Fontsize', 40)
    ylabel('Range (pixel)', 'Fontsize', 40)
    plot_para('Maximize',true,'Filename',name, 'Ratio', [4 3 1]);
    movefile([name, '.jpg'],  'output/')
end