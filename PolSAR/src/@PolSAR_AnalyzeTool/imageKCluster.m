function [labels] = imageKCluster(obj, im, nColors)
    % IMAGEKCLUSTER Implement the classification by using k_means 
    %
    % Syntax: IMAGEKCLUSTER(kai_1,kai_2,kai_3)
    %
    % Inputs:
	%    im      - image for k-means cluster.
    %    nColors - numbers of center of k-means cluster.
    %
    % Other m-files required: 
    % Subfunctions:
    % MAT-files required: none
    %
    % Reference 
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------

    lab_I = rgb2lab(im);
    ab = lab_I(:,:,2:3);
    nrows = size(ab,1);
    ncols = size(ab,2);
    ab = reshape(ab, nrows*ncols, 2);
    %ab = reshape(lab_I, size(ab,1)* size(ab,2),3);
    nColors = 3;
    % repeat the clustering 3 times to avoid local minima
    [cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','cityblock', ...
                                          'Replicates',5,'MaxIter',100);
    labels = uint8(reshape(cluster_idx,nrows,ncols));
end