function fourComponentDecomposition(obj, varargin)
    % FOURCOMPONENTDECOMPOSITION implements the original 4-component decomposition.
    %
    % syntax:
    % FOURCOMPONENTDECOMPOSITION(varargin)
    % Specifies function properties using one or more Name,Value pair arguments. 
    % Name-value pair settings apply to all the lines plotted. 
    %
    % Inputs:
    % Name-Value Pair Arguments:
    %    ('vv', s) - The pixel above SNR thershold is marked. If s is ture, 
    %                function will show the figure.
    %    ('hh', s) - The pixel above SNR thershold is marked. If s is ture, 
    %                function will show the figure.
    %    ('hv', s) - The pixel above SNR thershold is marked. If s is ture, 
    %                function will show the figure.
    %
    % Example: 
    %    1.
    %
    % Other m-files required: plot_para
    % Subfunctions: none
    % MAT-files required: none
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------

    parse_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'function_handle'},{});
    validationFcn_2_ = @(x) validateattributes(x,{'numeric'},{});
    validationFcn_3_ = @(x) validateattributes(x,{'char'},{});
    addParameter(parse_,'Plotsetting',@(x) 0,validationFcn_1_);
    addParameter(parse_,'Powrange',[],validationFcn_2_);
    addParameter(parse_,'Filename','',validationFcn_3_);
    parse(parse_,varargin{:})
    
    % f_ is the scattering matrix coefficient. The subscript
    % f_s: surface, f_d: double-bounce, f_v: volume, f_c: helix 
    % f_t: total power
    f_t = obj.obj.hh_hh + obj.vv_vv + 2*obj.hv_hv;
    f_c = 2 *abs(imag(obj.hh_hv - conj(obj.hv_vv)));

    % Decide which volume scattering model is used.
    % With three interval (-infty,-2dB) [-2dB, 2dB] (2dB, infty)
    crit_fv = 10*log10(obj.vv_vv./obj.hh_hh);
    image_range = ~(ones(size(obj.hh_hh)).*isnan(crit_fv)); % If the image is rotated,
    
    % thershold < -2dB
    f_v = (15/2*obj.hv_hv-15/8*f_c).*(crit_fv < -2);
    % -2dB < thershold < 2dB
    f_v = f_v + (8*obj.hv_hv-2*f_c).*(crit_fv >= -2).*(crit_fv<=2);
    % thershold > 2 dB
    f_v = f_v + (15/2*obj.hv_hv-15/8*f_c).*(crit_fv > 2);
    f_v = f_v.*image_range;
    % If volume power is 0 or negative, set helix power to 0
    f_c(f_v<=0) = 0;
    % Re-compute those f_v <= 0 point.
    f_v = zeros(size(obj.hh_hh));
    f_v = (15/2*obj.hv_hv-15/8*f_c).*(crit_fv < -2);% thershold < -2dB
    f_v = f_v + (8*obj.hv_hv-2*f_c).*(crit_fv >= -2).*(crit_fv<=2);% -2dB < thershold < 2dB
    f_v = f_v + (15/2*obj.hv_hv-15/8*f_c).*(crit_fv > 2);% thershold > 2 dB
    f_v = f_v.*image_range;
    
    % Decide double scattering or single scattering domainate 
    S = (1/2*(obj.hh_hh+obj.vv_vv+2*real(obj.hh_vv))-1/2*f_v).*((crit_fv < -2)+(crit_fv > 2)); % thershold < -2dB
    S = S + (1/2*(obj.hh_hh+obj.vv_vv+2*real(obj.hh_vv))-4*obj.hv_hv+f_c).*(crit_fv >= -2).*(crit_fv <= 2); % -2dB < thershold < 2dB
    S = S.*image_range;
    D = (1/2*(obj.hh_hh+obj.vv_vv-2*real(obj.hh_vv))-7/4*obj.hv_hv-1/16*f_c).*((crit_fv < -2)+(crit_fv > 2)); % thershold < -2dB or thershold > 2dB
    D = D + (1/2*(obj.hh_hh+obj.vv_vv-2*real(obj.hh_vv))-2*obj.hv_hv).*(crit_fv >= -2).*(crit_fv <= 2); % -2dB < thershold < 2dB
    D = D.*image_range;
    C = (1/2*(obj.hh_hh-obj.vv_vv-1j*2*imag(obj.hh_vv))-1/6*f_v).*(crit_fv<-2);
    C = C + 1/2*(obj.hh_hh-obj.vv_vv-1j*2*imag(obj.hh_vv)).*(crit_fv >= -2).*(crit_fv <= 2);
    C = C + (1/2*(obj.hh_hh-obj.vv_vv-1j*2*imag(obj.hh_vv))+1/6*f_v).*(crit_fv>2);
    C = C.*image_range;
    clear crit_fv
    % check if the helix power + volume power is smaller than total power
	% check C_0, to single-bounce dominates and double-bounce dominates
    C_0 = obj.hh_vv - obj.hv_hv + 1/2*f_c; 
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
            parse_.Results.Plotsetting(parse_.Results.Powrange, 'Colorbar_unit', "(dB)")
			%title('single', 'Interpreter', 'latex')
            plot_para('Maximize',true,'Filename', [parse_.Results.Filename,'_s']);
        figure 
            imagesc(10*log10(abs(f_d)))
            parse_.Results.Plotsetting(parse_.Results.Powrange, 'Colorbar_unit', "(dB)")
			%title('double', 'Interpreter', 'latex')
            plot_para('Maximize',true,'Filename',[parse_.Results.Filename, '_d']);
        figure
            imagesc(10*log10(abs(f_v)))
            parse_.Results.Plotsetting(parse_.Results.Powrange, 'Colorbar_unit', "(dB)")
			%title('volume', 'Interpreter', 'latex')
            plot_para('Maximize',true, 'Filename', [parse_.Results.Filename, '_v']);
        figure
            imagesc(10*log10(abs(f_c)))
            parse_.Results.Plotsetting(parse_.Results.Powrange, 'Colorbar_unit', "(dB)")
			%title('volume', 'Interpreter', 'latex')
            plot_para('Maximize',true, 'Filename', [parse_.Results.Filename, '_c']);
    end

    if(1)	% Plot the 4-component decomposition.
        img = single(zeros([size(obj.hh_hh), 3]));
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
            plot_para('Filename',parse_.Results.Filename,'Maximize',true);
        clear img 
    end
end