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
	N = size(x,1);
    DIM = size(x,2);
	c = zeros(clusterNum, DIM); % cluster centroid
	u = rand(N, clusterNum); % initialzie membership matrix
    u = u./sum(u,2);
    m = 2;
    % Algorithm parameters
    epsilon = 1e-5;
    
    [u, c] = GFCM(obj, x, m, c, u, N, clusterNum, max_iter, epsilon);
    [~, label] = max(u,[],2);
	%label = reshape(label, obj.IMAGE_SIZE);
end

%% Generalized FCM
function [u, c] = GFCM(obj, x, m, c, u, N, clusterNum, max_iter, epsilon)
	window_size = 5;
    obf = 0; obf_prev = 0;
    weight = zeros(N+1,window_size^2);
    row = uint32(zeros(N+1,window_size^2));
    x = [x; zeros(1,size(x,2))];
    dim = size(x,2);
    parfor it_n = 1 : N
        [w, r] = spatialDist(obj, it_n, window_size);
        weight(it_n,:) = w;
        row(it_n,:) = r;
    end    
    t_w = weight(1:end-1,:).';
    t_r = row(1:end-1,:).';
    for it = 1 : max_iter
        %% update center

        for it_c = 1 : clusterNum
            % Straight Forward
            %{
            dividend = 0; divisor = 0;
            parfor it_n = 1 : N
                [w, r] = spatialDist(obj, it_n, window_size);
                dividend= dividend + (u(it_n, it_c).^m) *sum(x(r,:).* ...
                                           repmat(w.',[1,dim]), 1);
                
                divisor = divisor + sum(sum(w)*u(it_n, it_c).^m);
            end
            %}
            % Matrix version
            
            div = repmat(t_w(:),[1, dim]).*x(t_r(:),:);
            div = reshape(sum(reshape(div, window_size^2, [])), [], dim);
            div = sum(repmat(u(:,it_c).^m,[1, dim]).*div, 1);
			c(it_c,:) = div/sum((u(:,it_c).^m).*sum(weight(1:end-1,:),2));
        end

        %% update membership (phase 1)
        dist = pdist2(x,c,'euclidean'); % dist(i,j) is the distance
                                        % between observation i in X 
                                        % and observation j in Y.
        % Straight forward
        %{
        for it_n = 1 : N
            [w, r] = spatialDist(obj, it_n, window_size);
            for it_c = 1 : clusterNum
                % denominator
                temp = 0;
                for it_cc = 1 : clusterNum
                    temp = temp + sum((w.*dist(r, it_cc)).^(1/(1-m)));
                end
                u(it_n, it_c) = sum(dist(r, it_c).*w).^(1/(1-m))/temp;
            end
        end
        %}
        % Matrix version
        for it_c = 1 : clusterNum
           div = t_w(:).*dist(t_r(:),it_c);
           u(:,it_c) = reshape(sum(reshape(div, window_size^2, [])),...
                                    [], dim).^(1/(1-m));
        end
        u = u./(repmat(sum(u,2), [1, clusterNum]));

        %% update membership (phase 2)
        % matrix version
        for it_c = 1 : clusterNum
            temp = t_w(:).u(t_r(:),it_c);
            u(:,it_c) = reshape(sum(reshape(temp, window_size^2, [])),...
                                    [], dim);
        end
        u = u./(repmat(sum(u,2), [1, clusterNum]));
        
        % objective function
        dist4obj = zeros(N,clusterNum);
        for it_c = 1 : clusterNum
            temp = t_w(:).*dist(t_r(:),it_c);
            dist4obj(:,it_c) = reshape(sum(reshape(temp, window_size^2,...
                                [])), [], dim);
        end
        obf = sum(sum(((u.^m).*dist4obj)));
        if abs(obf-obf_prev) < epsilon
            break;
        end
        obf_prev = obf;
    end
 
end

%% Hierarchical FCM
function [u, c] = HFCM(obj, x, m, c, u, N, clusterNum, max_iter, epsilon)
end

%% Generalized Hierarchical FCM
function [u, c] = GHFCM(obj, x, m, c, u, N, clusterNum, max_iter, epsilon)
end

%% Conventional FCM
function [u, c] = CFCM(obj, x, m, c, u, N, clusterNum, max_iter, epsilon)
    obf = 0;
    obf_prev = 0;
    halt = 0;
    for it = 1 : max_iter
		for it_c = 1 : clusterNum
			c(it_c,:) = sum(x.*(u(:,it_c).^m), 1) ...
				./sum(u(:,it_c).^m);
		end
        dist = pdist2(x,c,'euclidean');
        temp = dist.^(1/(1-m));
        u = temp ./ repmat(sum(temp,2), [1, size(x, 2)]);
        u = u./sum(u,2);
        obf = sum(sum((u.^m).*dist));
		%fprintf('@Iteration count = %i, obj. fcn = %f \n', it, obf);
        if abs(obf-obf_prev) < epsilon
            halt = halt+1;
            if halt> 3
                break;
            end
        end
        obf_prev = obf;
    end
end

%% ------------- Utility function -------------- %%
function [weight, row_n] = spatialDist(obj, n_1, window_size)
	% Inputs:
    %	* n_1 : nth data
	% Outputs:
    %	* weight : the distance weight array with size 1xN  
    %   * row_n  : A array with size 1xN which containing
    %              the row number of the pixel in the window.
    %------------- <<<<< >>>>>--------------
	
	sigma = (window_size-1)/4; % Spatial relation
	box = -(window_size-1)/2:1:(window_size-1)/2;
	[X,Y] = meshgrid(box,box);
	% obtain the index of the center of the windows.
	[row_, col_]= idx2nm(obj, n_1); % idx = (row, col)
	% obtain the index inside windows
	idx = [col_+X(:), row_+Y(:)]; 
    
	% exclude the elements which is not in the matrix
	idx = idx(( (idx(:,1)<=0)+(idx(:,1)>obj.IMAGE_SIZE(1))+...
			    (idx(:,2)<=0)+(idx(:,2)>obj.IMAGE_SIZE(2)) )==0, :); 
	% obtain the distance between the center and the elements in windows
	dist = pdist2([row_, col_], idx);
	% obtain the weight of each element
    weight = 1/sqrt(2*pi*sigma^2)*exp(-dist/(2*sigma^2));
    padding = zeros(1,window_size^2-numel(weight));
    weight = [weight, padding];
	% obtain the row number of each elements
	row_n = uint32([nm2idx(obj, idx(:,2), idx(:,1)).', ...
        obj.IMAGE_SIZE(1)*obj.IMAGE_SIZE(2)+1+padding]);
end
function [n, m] = idx2nm(obj, ind)
    % Inputs:
    %  * ind : nth data
    n = mod(ind-1, obj.IMAGE_SIZE(1))+1; % row 
    m = floor((ind-1)/obj.IMAGE_SIZE(1))+1; % column
end
function [row] = nm2idx(obj, n, m)
    % Inputs:
    %  * n : row 
    %  * m : column
	row = (m-1)*obj.IMAGE_SIZE(1)+n;
end