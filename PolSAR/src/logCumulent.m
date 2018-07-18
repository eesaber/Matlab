function [mu_1, kai_2] = logCumulent(A)
    global im_size
    [row,~,pixel_size] = size(A);
    if row == 3
        m_det = @(x) x(1,1,:).*x(2,2,:).*x(3,3,:)+x(1,2,:).*x(2,3,:).*x(3,1,:)+x(1,3,:).*x(2,1,:).*x(3,2,:)......
                -(x(1,3,:).*x(2,2,:).*x(3,1,:)+x(1,1,:).*x(2,3,:).*x(3,2,:)+x(1,2,:).*x(2,1,:).*x(3,3,:));
    else
        m_det = @(x) x(1,1,:).*x(2,2,:) - x(1,2,:).*x(2,1,:);
    end
    mask = ones(9,9)/9^2;
    mu_1 = real(squeeze(log(m_det(A))));
    mu_1 = reshape(mu_1, im_size);
    mu_1 = conv2(mu_1, mask, 'same');
    
    mu_2 = real(squeeze(log(m_det(A)).^2));
    mu_2 = reshape(mu_2, im_size);
    mu_2 = conv2(mu_2, mask, 'same');
    
    kai_2 = mu_2 - mu_1.^2;
    
    %{
    mu_1 = zeros(1,pixel_size);
    mu_2 = zeros(1,pixel_size);
    kai_2 = zeros(1,pixel_size);
    parfor n = 1 : pixel_size
        mu_1(n) = det(log(A(:,:,n)));
        mu_2(n) = mu_1(n)^2;
        kai_2(n) = mu_2(n) - mu_1(n)^2;
    end
    %}
    
    kai_2 = reshape(kai_2, im_size);
end