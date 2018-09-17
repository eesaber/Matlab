function [label, c] = myFCM (obj, x, clusterNum, max_iter)
    % MYFCM implement the fuzzy c-means clustering
    %
    % Syntax:
	%	* [label] = MYFCM(x, clusterNum, max_iter)
	%
    % Inputs:
    %	* x           : Traing data
    %	* clusterNum  : number of center of cluster.
	%   * m           : m is fuzzy partition matrix exponent for controlling
	%                   the degree of fuzzy overlap, with m > 1. Fuzzy overlap
	%					refers to how fuzzy the boundaries between clusters 
	%					are, that is the number of data points that have 
	%				    significant membership in more than one cluster.
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------
    
    % Data 
    DIM = size(x,2);
	N = size(x,1);
	c = zeros(clusterNum, DIM); % cluster centroid
	u = rand(N, clusterNum); % initialzie membership matrix
    u = u./repmat(sum(u,2), [1, DIM]);	
    % Algorithm parameters
    epsilon = 1e-5;
    
    label = uint8(max(u,[],2));
	%label = reshape(label, obj.IMAGE_SIZE);
end
%% Generalized FCM
function [u, c] = GFCM(obj, DIM, N, c, u, epsilon)
	window_size = 5;
    obf = 0; obf_prev = 0;
	weight = zeros(N,1);
	for it_n : 1 : N
		[w, row_n] = spatialDist(obj, it_n, window_size);
	end
	for it = 1 : max_iter
		for it_c = 1 : clusterNum	
			c(it_c,:) = sum(x.*repmat(u(:,it_c).^m,[1, DIM]),1) ...
				./sum(u(:,it_c).^m);
		end
        dist = pdist2(x,c,'euclidean');
        temp = dist.^(1/(1-m));
        u = temp ./ repmat(sum(temp,2),[1, DIM]);
        obf = sum(sum(u.^m.*dist));
		%fprintf('Iteration count = %i, obj. fcn = %f \n', it, obf);
        if abs(obf-obf_prev) < epsilon
            break;
        end
        obf_prev = obf;
    end
end
%% Hierarchical FCM
function [u, c] = HFCM(obj, DIM, N, c, u, epsilon)
end
%% Generalized Hierarchical FCM
function [u, c] = GHFCM(obj, DIM, N, c, u, epsilon)
end
%% Conventional FCM
function [u, c] = FCM(obj, DIM, N, c, u, epsilon)
    obf = 0;
    obf_prev = 0;
	for it = 1 : max_iter    
		for it_c = 1 : clusterNum
			c(it_c,:) = sum(x.*repmat(u(:,it_c).^m,[1, DIM]),1) ...
				./sum(u(:,it_c).^m);
		end
        dist = pdist2(x,c,'euclidean');
        temp = dist.^(1/(1-m));
        u = temp ./ repmat(sum(temp,2),[1, DIM]);
        obf = sum(sum(u.^m.*dist));
		%fprintf('Iteration count = %i, obj. fcn = %f \n', it, obf);
        if abs(obf-obf_prev) < epsilon
            break;
        end
        obf_prev = obf;
    end
end
%% ------------- Utility function -------------- %%
function [weight, row_n] = spatialDist(obj, n_1, window_size)
	% Inputs:
    %	* n_1 : scalar
	% Outputs:
    %	* weight : 
    %------------- <<<<< >>>>>--------------
	
	sigma = (window_size-1)/4; % Spatial relation
	box = -(window_size-1/2):(window_size-1/2);
	[X,Y] = meshgrid(box,box);
	% obtain the index of the center of the windows.
	idx = idx2ij(obj, n_1);
	% obtain the index inside windows
	idx = [idx(1)+X(:), idx(2)+Y(:)];
	% exclude the elements which is not in the matrix
	idx = idx(~idx(idx(:,1)<0 || idx(:,1)>obj.IMAGE_SIZE(1) || ...
			idx(:,2)<0 || idx(:,2)>obj.IMAGE_SIZE(2))); 
	% obtain the distance between the center and the elements in windows
	dist = pdist2(idx2ij(obj, n_1), idx);
	% obtain the weight of each element
    weight = 1/sqrt(2*pi*sigma^2)*exp(-dist/(2*sigam^2));
	% obtain the row number of each elements
	row_n = ij2idx(obj, idx);
end
function [i, j] = idx2ij(obj, ind)
    i = mod(ind-1, obj.IMAGE_SIZE(1))+1;
    j = floor((ind-1)/obj.IMAGE_SIZE(1))+1;
end
function [n] = ij2idx(obj, i, j)
	n = i*obj.IMAGE_SIZE(1)+j*obj.IMAGE_SIZE(2);
end