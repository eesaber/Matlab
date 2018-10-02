function varargout = logCumulant(obj)
    % LOGCUMULANT calculate the first, second and third order log
    % cumulant of a 2x2 matrix.
    %
    % Syntax:
	%	[mu_1, kai_2, kai_3] = LOGCUMULANT()
	%
    % Description:
    %   LOGCUMULANT()
    %   * Because S_hv is very noisy (below SNR), the 2x2 covariance 
    %     matrix is used. A 3x3 covaraince can be used by setting row to 3.
    %
    %   * If IS_BIGFILE is true, the elements of covaraince matrix will
    %     be set to zero.
    %
    % Outputs:
    %    mu_1  - First order log-cumulants.
    %    kai_2 - Second order log-cumulants.
    %    kai_3 - Third order log-cumulants.
    %
    % Other m-files required: none
    % Subfunctions: none
    % MAT-files required: none
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------

    minArgs=1;  
    maxArgs=6;
    nargoutchk(minArgs,maxArgs)
    row = 2;
    size_N = obj.IMAGE_SIZE(1)*obj.IMAGE_SIZE(2); 
    A = cat(1,cat(2, reshape(obj.hh_hh,[1,1,size_N]), reshape(obj.hh_vv,[1,1,size_N])), ...
        cat(2,reshape(conj(obj.hh_vv),[1,1,size_N]), reshape(obj.vv_vv,[1,1,size_N])));
    if obj.IS_BIGFILE
        setCov2Zero(obj)
    end
    if row == 3
        m_det = @(x) x(1,1,:).*x(2,2,:).*x(3,3,:)+x(1,2,:).*x(2,3,:).*x(3,1,:)+x(1,3,:).*x(2,1,:).*x(3,2,:)......
                -(x(1,3,:).*x(2,2,:).*x(3,1,:)+x(1,1,:).*x(2,3,:).*x(3,2,:)+x(1,2,:).*x(2,1,:).*x(3,3,:));
    else
        m_det = @(x) x(1,1,:).*x(2,2,:) - x(1,2,:).*x(2,1,:);
    end
    filtsize = 3;
    mask = ones(3,3)/(3*3);
    % μ_1
    mu_1 = real(squeeze(log(m_det(A))));
    mu_1 = reshape(mu_1, obj.IMAGE_SIZE);
    mu_1 = imboxfilt(mu_1,filtsize,'Padding','replicate');
    % μ_2
    mu_2 = real(squeeze(log(m_det(A)).^2));
    mu_2 = reshape(mu_2, obj.IMAGE_SIZE);
    mu_2 = imboxfilt(mu_2,filtsize,'Padding','replicate');
    % μ_3
    mu_3 = real(squeeze(log(m_det(A)).^3));
    mu_3 = reshape(mu_3, obj.IMAGE_SIZE);
    mu_3 = imboxfilt(mu_3,filtsize,'Padding','replicate');
    % μ_4
    mu_4 = real(squeeze(log(m_det(A)).^4));
    mu_4 = reshape(mu_4, obj.IMAGE_SIZE);
    mu_4 = imboxfilt(mu_4,filtsize,'Padding','replicate');
    % μ_5
    mu_5 = real(squeeze(log(m_det(A)).^5));
    mu_5 = reshape(mu_5, obj.IMAGE_SIZE);
    mu_5 = imboxfilt(mu_5,filtsize,'Padding','replicate');
    % μ_6
    mu_6 = real(squeeze(log(m_det(A)).^6));
    mu_6 = reshape(mu_6, obj.IMAGE_SIZE);
    mu_6 = imboxfilt(mu_6,filtsize,'Padding','replicate');
    % κ_2 to κ_6
    kai_2 = mu_2 - mu_1.^2;
    kai_3 = mu_3 - 3*mu_1.*mu_2 + 2*mu_1.^3;
    kai_4 = mu_4 - 4*mu_3.*mu_1 - 3*mu_2.^2 + 12*mu_2.*mu_1.^2 - 6*mu_1.*4;
    kai_5 = mu_5 - 5*mu_1.*mu_4 - 10*mu_2.*mu_3 + 20*mu_1.^2.*mu_3 ...
            + 30*mu_1.*mu_2.^2 - 60*mu_1.^3.*mu_2 + 24*mu_1.^5;
    kai_6 = mu_6 - 6*mu_1.*mu_5 - 15*mu_2.*mu_4 + 30*mu_1.^2.*mu_4 ...
        - 10*mu_3.^2 + 120*mu_1.*mu_2.*mu_3 - 120*mu_1.^3.*mu_3 ...
        + 30*mu_2.^3 - 270*mu_1.^2.*mu_2.^2 + 360*mu_1.^4.*mu_2 - 120*mu_1.^6;
    showIm(obj, mu_1, kai_2, kai_3) % Display RGB image
    % showIm(obj, mu_1, kai_2, kai_4) % Display RGB image
    % showIm(obj, mu_1, kai_3, kai_4) % Display RGB image
    % showIm(obj, kai_2, kai_3, kai_4) % Display RGB image
    
    % showIm(obj, mu_1, kai_2, kai_3) % Display RGB image
    if nargout>=1, varargout{1} = mu_1;
        if nargout>=2, varargout{2} = kai_2;
            if nargout>=3, varargout{3} = kai_3;
                if nargout>=4, varargout{4} = kai_4;
                    if nargout>=5, varargout{5} = kai_5;
                        if nargout>=6, varargout{6} = kai_6;
                        end
                    end
                end
            end
        end
    end
end
function showIm(obj, kai_1, kai_2, kai_3)
    figure
    if 1
        image(cat(3, reshape(intensityMapping(kai_1,'Bit','uint16','Method','Sigma','std',2.5), size(kai_1)),...
            reshape(intensityMapping(kai_2,'Bit','uint16','Method','Sigma','std',2.5), size(kai_2)), ...
            reshape(intensityMapping(kai_3,'Bit','uint16','Method','Sigma','std',2.5), size(kai_3))))
    elseif 0
        image(cat(3, histeq(rescale(kai_1)), histeq(rescale(kai_2)), histeq(rescale(kai_3))))
    else
        image(cat(3, reshape(intensityMapping(kai_1,'Bit','uint16','Method','truncate','Range',[-14 -4]), size(kai_1)),...
        reshape(intensityMapping(kai_2,'Bit','uint16','Method','truncate','Range',[0 3]), size(kai_2)), ...
        reshape(intensityMapping(kai_3,'Bit','uint16','Method','truncate','Range',[-1 1]), size(kai_3))))
    end
    set(gca,'Ydir','normal')
    plot_para('Maximize',true,'Filename',[obj.OUTPUT_PATH, '/logcumulant'])
end