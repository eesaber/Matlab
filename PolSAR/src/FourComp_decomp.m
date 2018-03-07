function FourComp_decomp(hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv, filename)
% FOUCOMP_DECOMP implements the original 4-component decomposition.
%
% FOUCOMP_DECOMP(C), C is a matrix with size row x col x 6, and C(:,:,1) is
% |S_{hh}|^2, C(:,:,2) is |S_{hv}|^2, C(:,:,3) is |S_{vv}|^2, C(:,:,4) is
% S_{hh}S^*_{hv}, C(:,:,5) is S_{hh}S_{vv}^2 and C(:,:,6) is
% S_{hv}S^*_{vv}.

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
    % thershold < -2dB
    f_v = (15/2*hv_hv-15/8*f_c).*(crit_fv < -2);
    % -2dB < thershold < 2dB
    f_v = f_v + (8*hv_hv-2*f_c).*(crit_fv >= -2).*(crit_fv<=2);
    % thershold > 2 dB
    f_v = f_v + (15/2*hv_hv-15/8*f_c).*(crit_fv > 2);
    f_v = f_v.*image_range;
    
    %{
    f_v(logical((f_v<=0).*(crit_fv < -2))) = 15/2*hv_hv(logical((f_v<=0).*(crit_fv < -2))); % thershold < -2dB
    f_v(logical((f_v<=0).*(crit_fv >= -2).*(crit_fv<=2))) = 8*hv_hv(logical((f_v<=0).*(crit_fv >= -2).*(crit_fv<=2))); % -2dB < thershold < 2dB
    f_v(logical((f_v<=0).*(crit_fv > 2))) = 15/2*hv_hv(logical((f_v<=0).*(crit_fv > 2))); % thershold > 2 dB
    %}
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
    temp = f_t > (f_c + f_v);
    f_v(temp) = f_t(temp)-f_c(temp);
    
    C_0 = hh_vv - hv_hv + 1/2*f_c; 
    f_s = (S - abs(C).^2./D).*(real(C_0)<0);
    f_d = (D - abs(C).^2./D).*(real(C_0)<0);
    f_s = f_s + (S + abs(C).^2./S).*(real(C_0)>=0);
    f_d = f_d + (D - abs(C).^2./S).*(real(C_0)>=0);
    f_s = f_s.*temp.*image_range;
    f_d = f_d.*temp.*image_range;
    f_s(~temp) = 0;
    f_d(~temp) = 0;
    clear S D C
    % Set the negative power to 0    
    ch_ind = f_d<0;
    f_s(ch_ind) = f_t(ch_ind) - f_v(ch_ind) - f_c(ch_ind);
    f_d(ch_ind) = 0;
    ch_ind = f_s<0; 
    f_d(ch_ind) = f_t(ch_ind) - f_v(ch_ind) - f_c(ch_ind);
    f_s(ch_ind) = 0;
      
    if(1)
        figure
            imagesc(10*log10(f_s))
            xlabel('Azimuth')
            set(gca,'Ydir','normal'), colorbar, colormap jet;
            caxis([-30 20])
            plot_para('Maximize',true)
        figure 
            imagesc(10*log10(f_d))
            xlabel('Azimuth (pixel)', 'Interpreter', 'latex')
            ylabel('Range (pixel)', 'Interpreter', 'latex')
            set(gca,'Ydir','normal'), colorbar, colormap jet;
            caxis([-30 20])
            plot_para('Maximize',true)
        figure
            imagesc(10*log10(f_v))
            xlabel('Azimuth (pixel)', 'Interpreter', 'latex')
            ylabel('Range (pixel)', 'Interpreter', 'latex')
            set(gca,'Ydir','normal'), colorbar, colormap jet;
            caxis([-30 10])
            plot_para('Maximize',true)
    end

    if(1)	% Plot the 4-component decomposition.
        img = single(zeros([size(hh_hh), 3]));
        up_ = 20; low_ = -30;
        img(:,:,1) = 10*log10(f_d);
        img(:,:,2) = 10*log10(f_v);
        img(:,:,3) = 10*log10(f_s);
        %%clear P_d  f_v P_s f_c f_d f_s 
        img(img < low_) = low_;
        img(img > up_) = up_;
        img = (img-low_)/(up_-low_);
        figure(10)
            image(img)
            set(gca,'Ydir','normal')
            xlabel('Azimuth (pixel)', 'Interpreter', 'latex')
            ylabel('Range (pixel)', 'Interpreter', 'latex')
            plot_para('Filename',filename,'Maximize',true,'Ratio', [4 3 1]);
        movefile([filename, '.jpg'], 'output/')
        clear img 
    end

end