function [mu_1, kai_2, kai_3] = logCumulant(obj)
    row = 2;
    size_N = obj.IMAGE_SIZE(1)*obj.IMAGE_SIZE(2); 
    A = cat(1,cat(2, reshape(obj.hh_hh,[1,1,size_N]), reshape(obj.hh_vv,[1,1,size_N])), ...
        cat(2,reshape(conj(obj.hh_vv),[1,1,size_N]), reshape(obj.vv_vv,[1,1,size_N])));
    if row == 3
        m_det = @(x) x(1,1,:).*x(2,2,:).*x(3,3,:)+x(1,2,:).*x(2,3,:).*x(3,1,:)+x(1,3,:).*x(2,1,:).*x(3,2,:)......
                -(x(1,3,:).*x(2,2,:).*x(3,1,:)+x(1,1,:).*x(2,3,:).*x(3,2,:)+x(1,2,:).*x(2,1,:).*x(3,3,:));
    else
        m_det = @(x) x(1,1,:).*x(2,2,:) - x(1,2,:).*x(2,1,:);
    end
    mask = ones(2,4)/(2*4);
    % μ_1
    mu_1 = real(squeeze(log(m_det(A))));
    mu_1 = reshape(mu_1, obj.IMAGE_SIZE);
    mu_1 = conv2(mu_1, mask, 'same');
    % μ_2
    mu_2 = real(squeeze(log(m_det(A)).^2));
    mu_2 = reshape(mu_2, obj.IMAGE_SIZE);
    mu_2 = conv2(mu_2, mask, 'same');
    % μ_3
    mu_3 = real(squeeze(log(m_det(A)).^3));
    mu_3 = reshape(mu_3, obj.IMAGE_SIZE);
    mu_3 = conv2(mu_3, mask, 'same');
    % κ_2 and κ_3
    kai_2 = mu_2 - mu_1.^2;
    kai_3 = mu_3 - 3*mu_1.*mu_2 + 2*mu_1.^3;
end