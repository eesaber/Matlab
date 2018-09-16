function [label, u] = myFCM (obj, x_input, clusterNum, max_iter)
    % MYFCM implement the fuzzy c-means clustering
    %
    % Syntax:
	%	* [label] = MYFCM(x, clusterNum, max_iter)
	%
    % Inputs:
    %	* x           : Traing data
    %	* clusterNum  : number of center of cluster.
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------
    
    % Data 
    c = ones(clusterNum, 1); % initialize cluster centroid
    u = rand(numel(x), clusterNum); % initialzie membership matrix
    DIM = shape(x_input,3);
    N =  numel(x_input)/dim;
    u = u./repmat(sum(u,2), [1, dim]);
    x = reshape(x_input, [N, DIM]);
    
    % Algorithm parameters
    epsilon = 1e-6;
    obf = 0;
    obf_prev = 0;
    %my_dis
    
    for it = 1 : max_iter    
        c = sum(u.*x,1) ./sum(u,1);
        dist = pdist2(x,c,'euclidean');
        temp = dist.^(1/(1-clusterNum));
        u = temp ./ repmat(sum(temp,2),[1, DIM]);
        obf = sum(sum(u.*d));
        if abs(obf-obf_prev) < epsilon
            break;
        end
        obf_prev = obf;
    end
    [~, label] = max(u,2);
    label = uint8(reshape(label, shape(x_input(:,:,1))));
    
end
function [w] = spatialDist(obj, n_1, n_2, row)
    dist = norm(ind2ij(n_1, row) - ind2ij(n_2, row),2);
    w = 1/sqrt(2*pi*sigma^2)*exp(-dist/(2*sigam^2));
end
function [i, j] = ind2ij(ind, m)
    i=mod(ind-1,m)+1;
    j=floor((ind-1)/m)+1;
end