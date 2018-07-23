function [mu_1, kai_2, kai_3] = Logcumulant(A)
    global im_size
    [row,~,pixel_size] = size(A);
    if row == 3
        m_det = @(x) x(1,1,:).*x(2,2,:).*x(3,3,:)+x(1,2,:).*x(2,3,:).*x(3,1,:)+x(1,3,:).*x(2,1,:).*x(3,2,:)......
                -(x(1,3,:).*x(2,2,:).*x(3,1,:)+x(1,1,:).*x(2,3,:).*x(3,2,:)+x(1,2,:).*x(2,1,:).*x(3,3,:));
    else
        m_det = @(x) x(1,1,:).*x(2,2,:) - x(1,2,:).*x(2,1,:);
    end
    mask = ones(9,9)/9^2;
    % μ_1
    mu_1 = real(squeeze(log(m_det(A))));
    mu_1 = reshape(mu_1, im_size);
    mu_1 = conv2(mu_1, mask, 'same');
    % μ_2
    mu_2 = real(squeeze(log(m_det(A)).^2));
    mu_2 = reshape(mu_2, im_size);
    mu_2 = conv2(mu_2, mask, 'same');
    % μ_3
    mu_3 = real(squeeze(log(m_det(A)).^3));
    mu_3 = reshape(mu_3, im_size);
    mu_3 = conv2(mu_3, mask, 'same');
    % κ_2 and κ_3
    kai_2 = mu_2 - mu_1.^2;
    kai_3 = mu_3 - 3*mu_1.*mu_2 + 2*mu_1;
end