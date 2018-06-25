function FourComp_decomp(hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv, filename, varargin)
% FOUCOMP_DECOMP implements the original 4-component decomposition.
%
% FOUCOMP_DECOMP(C), C is a matrix with size row x col x 6, and C(:,:,1) is
% |S_{hh}|^2, C(:,:,2) is |S_{hv}|^2, C(:,:,3) is |S_{vv}|^2, C(:,:,4) is
% S_{hh}S^*_{hv}, C(:,:,5) is S_{hh}S_{vv}^2 and C(:,:,6) is
% S_{hv}S^*_{vv}.
    parse_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'function_handle'},{});
    validationFcn_2_ = @(x) validateattributes(x,{'numeric'},{});
    addParameter(parse_,'Plotsetting',@(x) 0,validationFcn_1_);
    addParameter(parse_,'Powrange',[],validationFcn_2_);
    parse(parse_,varargin{:})
    
    % f_ is the scattering matrix coefficient. The subscript
    % f_s: surface, f_d: double-bounce, f_v: volume, f_c: helix 
    % f_t: total power
    f_t = hh_hh + vv_vv + 2*hv_hv;
    f_c = 2 *abs(imag(hh_hv - conj(hv_vv)));

    % Decide which volume scattering model is used.
    % With three interval (-infty,-2dB) [-2dB, 2dB] (2dB, infty)
    crit_fv = 10*log10(vv_vv./hh_hh);
    image_range = ~(ones(size(hh_hh)).*isnan(crit_fv)); % If the image is rotated,
    
    % thershold < -2dB
    f_v = (15/2*hv_hv-15/8*f_c).*(crit_fv < -2);
    % -2dB < thershold < 2dB
    f_v = f_v + (8*hv_hv-2*f_c).*(crit_fv >= -2).*(crit_fv<=2);
    % thershold > 2 dB
    f_v = f_v + (15/2*hv_hv-15/8*f_c).*(crit_fv > 2);
    f_v = f_v.*image_range;
    % If volume power is 0 or negative, set helix power to 0
    f_c(f_v<=0) = 0;
    % Re-compute those f_v <= 0 point.
    f_v = zeros(size(hh_hh));
    f_v = (15/2*hv_hv-15/8*f_c).*(crit_fv < -2);% thershold < -2dB
    f_v = f_v + (8*hv_hv-2*f_c).*(crit_fv >= -2).*(crit_fv<=2);% -2dB < thershold < 2dB
    f_v = f_v + (15/2*hv_hv-15/8*f_c).*(crit_fv > 2);% thershold > 2 dB
    f_v = f_v.*image_range;
    
    % Decide double scattering or single scattering domainate 
    S = (1/2*(hh_hh+vv_vv+2*real(hh_vv))-1/2*f_v).*((crit_fv < -2)+(crit_fv > 2)); % thershold < -2dB
    S = S + (1/2*(hh_hh+vv_vv+2*real(hh_vv))-4*hv_hv+f_c).*(crit_fv >= -2).*(crit_fv <= 2); % -2dB < thershold < 2dB
    S = S.*image_range;
    D = (1/2*(hh_hh+vv_vv-2*real(hh_vv))-7/4*hv_hv-1/16*f_c).*((crit_fv < -2)+(crit_fv > 2)); % thershold < -2dB or thershold > 2dB
    D = D + (1/2*(hh_hh+vv_vv-2*real(hh_vv))-2*hv_hv).*(crit_fv >= -2).*(crit_fv <= 2); % -2dB < thershold < 2dB
    D = D.*image_range;
    C = (1/2*(hh_hh-vv_vv-1j*2*imag(hh_vv))-1/6*f_v).*(crit_fv<-2);
    C = C + 1/2*(hh_hh-vv_vv-1j*2*imag(hh_vv)).*(crit_fv >= -2).*(crit_fv <= 2);
    C = C + (1/2*(hh_hh-vv_vv-1j*2*imag(hh_vv))+1/6*f_v).*(crit_fv>2);
    C = C.*image_range;
    clear crit_fv
    % check if the helix power + volume power is smaller than total power
	% check C_0, to single-bounce dominates and double-bounce dominates
    C_0 = hh_vv - hv_hv + 1/2*f_c; 
    f_s = (S - abs(C).^2./D).*(real(C_0)<0);
    f_d = (D + abs(C).^2./D).*(real(C_0)<0);
    f_s = f_s + (S + abs(C).^2./S).*(real(C_0)>=0);
    f_d = f_d + (D - abs(C).^2./S).*(real(C_0)>=0);
    f_s(f_t - f_c - f_v < 0) = 0;
    f_d(f_t - f_c - f_v < 0) = 0;
	% If both f_s and f_d is negative, set both of them to 0.
	temp = logical((f_s<0).*(f_d<0));
	f_s(temp) = 0;
	f_d(temp) = 0;
    f_v(f_t - f_c - f_v < 0) = f_t(f_t - f_c - f_v < 0) - f_c(f_t - f_c - f_v < 0);
	f_v(temp) = f_t(temp) - f_c(temp);
    clear S D C temp
    % Set the negative power to 0
    ind_d = logical((f_d<0).*(f_t > (f_c + f_v)));
    f_s(ind_d) = f_t(ind_d) - f_v(ind_d) - f_c(ind_d);
    f_d(ind_d) = 0;
    ind_d = logical((f_s<0).*(f_t > (f_c + f_v)));
    f_d(ind_d) = f_t(ind_d) - f_v(ind_d) - f_c(ind_d);
    f_s(ind_d) = 0;
    clear ind_d
    
    
    if(1)
        figure
            imagesc(10*log10(abs(f_s)))
            parse_.Results.Plotsetting(parse_.Results.Powrange, 1, 'Colorbar_unit', "(dB)")
			%title('single', 'Interpreter', 'latex')
            plot_para('Maximize',true,'Filename', [filename,'_s']);
        figure 
            imagesc(10*log10(abs(f_d)))
            parse_.Results.Plotsetting(parse_.Results.Powrange, 1, 'Colorbar_unit', "(dB)")
			%title('double', 'Interpreter', 'latex')
            plot_para('Maximize',true,'Filename',[filename, '_d']);
        figure
            imagesc(10*log10(abs(f_v)))
            parse_.Results.Plotsetting(parse_.Results.Powrange, 1, 'Colorbar_unit', "(dB)")
			%title('volume', 'Interpreter', 'latex')
            plot_para('Maximize',true, 'Filename', [filename, '_v']);
        figure
            imagesc(10*log10(abs(f_c)))
            parse_.Results.Plotsetting(parse_.Results.Powrange, 1, 'Colorbar_unit', "(dB)")
			%title('volume', 'Interpreter', 'latex')
            plot_para('Maximize',true, 'Filename', [filename, '_c']);
    end

    if(1)	% Plot the 4-component decomposition.
        img = single(zeros([size(hh_hh), 3]));
        up_ = 0; low_ = -40;
        img(:,:,1) = 10*log10(abs(f_d));
        img(:,:,2) = 10*log10(abs(f_v));
        img(:,:,3) = 10*log10(abs(f_s));
        %%clear P_d  f_v P_s f_c f_d f_s 
        img(img < low_) = low_;
        img(img > up_) = up_;
        img = (img-low_)/(up_-low_);
        figure
            image(img)
            %set(gca,'Ydir','normal','View',[90 90])
            set(gca,'Ydir','normal')
            xlabel('Azimuth (pixel)', 'Interpreter', 'latex')
            ylabel('Range (pixel)', 'Interpreter', 'latex')
            plot_para('Filename',filename,'Maximize',true);
        clear img 
    end

end