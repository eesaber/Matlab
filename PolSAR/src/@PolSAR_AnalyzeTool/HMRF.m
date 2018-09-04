function [X] = HMRF(obj, kai_1, kai_2, kai_3, labels, nColors)
    %% HMRF implement hidden Markov random field.

    EM_iter=10; % max num of iterations
    MAP_iter=10; % max num of iterations
    % fprintf('Performing k-means segmentation\n');
    % [X, mu, sigma]=image_kmeans(Y,k);
    % imwrite(uint8(X*120),'initial labels.png');
    im =cat(3, reshape(Stat(kai_1,'Bit','uint16'), size(kai_1)),...
        reshape(Stat(kai_2,'Bit','uint16'), size(kai_2)), ...
        reshape(Stat(kai_3,'Bit','uint16'), size(kai_3)));
    set(gca,'Ydir','normal')
    im = rgb2gray(im);
    Z = edge(im, 'Canny', 0.7);
    [X, mu, sigma]=HMRF_EM(labels, im , Z, mu, sigma, nColors, EM_iter, MAP_iter);
    
end