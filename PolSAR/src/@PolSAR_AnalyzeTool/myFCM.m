function [label, c] = myFCM (obj, x, clusterNum, max_iter, subClusterNum, algo)
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
    
    % input arugument
    switch nargin
        case 4
            subClusterNum  = 1;
            algo = 'CFCM';
        case 5
            algo = 'CFCM';
    end

    % Algorithm parameters
	N = size(x,1);
    rng(5)
	u = rand(N, clusterNum); % initialzie membership matrix
    u = u./sum(u,2);
    drg_m = 2;
    drg_n = 2;
    epsilon = 1e-10;
    
    switch algo
        case 'GFCM'
            [u, c] = GFCM(obj, x, u, drg_m, clusterNum, max_iter, epsilon);
        case 'HFCM'
            [u, c] = HFCM(obj, x, u, drg_m, drg_n, clusterNum, subClusterNum, max_iter, epsilon);
        case 'GHFCM'
            [u, c] = GHFCM(obj, x, u, drg_m, drg_n, clusterNum, subClusterNum, max_iter, epsilon);
        otherwise
            [u, c] = CFCM(obj, x, u, drg_m, clusterNum, max_iter, epsilon);
    end
    [~, label] = max(u,[],2);
    %label = reshape(label, obj.IMAGE_SIZE);
    
end

%% Generalized FCM
function [u, c] = GFCM(obj, x, u, drg_m, clusterNum, max_iter, epsilon)
    N_DATA = size(x,1);
    DIM = size(x,2);
	window_size = 5;
    obf = 0; obf_prev = 0;
    c = zeros(clusterNum, DIM); % cluster centroid
    weight = zeros(N_DATA+1,window_size^2);
    row = uint32(zeros(N_DATA+1,window_size^2));
    x = [x; zeros(1,size(x,2))];

    parfor it_n = 1 : N_DATA
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
            parfor it_n = 1 : N_DATA
                [w, r] = spatialDist(obj, it_n, window_size);
                dividend= dividend + (u(it_n, it_c).^drg_m) *sum(x(r,:).* ...
                                           repmat(w.',[1,DIM]), 1);
                
                divisor = divisor + sum(sum(w)*u(it_n, it_c).^drg_m);
            end
            %}
            % Matrix version
            
            div = repmat(t_w(:),[1, DIM]).*x(t_r(:),:);
            div = reshape(sum(reshape(div, window_size^2, [])), [], DIM);
            div = sum(repmat(u(:,it_c).^drg_m,[1, DIM]).*div, 1);
			c(it_c,:) = div/sum((u(:,it_c).^drg_m).*sum(weight(1:end-1,:),2));
        end

        %% update membership (phase 1)
        dist = pdist2(x,c,'euclidean').^2;; % dist(i,j) is the distance
                                        % between observation i in X 
                                        % and observation j in Y.
        % Straight forward
        %{
        for it_n = 1 : N_DATA
            [w, r] = spatialDist(obj, it_n, window_size);
            for it_c = 1 : clusterNum
                % denominator
                temp = 0;
                for it_cc = 1 : clusterNum
                    temp = temp + sum((w.*dist(r, it_cc)).^(1/(1-drg_m)));
                end
                u(it_n, it_c) = sum(dist(r, it_c).*w).^(1/(1-drg_m))/temp;
            end
        end
        %}
        % Matrix version
        for it_c = 1 : clusterNum
           div = t_w(:).*dist(t_r(:),it_c);
           u(:,it_c) = reshape(sum(reshape(div, window_size^2, [])),...
                                    [], DIM).^(1/(1-drg_m));
        end
        u = u./(repmat(sum(u,2), [1, clusterNum]));

        %% update membership (phase 2)
        % matrix version
        for it_c = 1 : clusterNum
            temp = t_w(:).u(t_r(:),it_c);
            u(:,it_c) = reshape(sum(reshape(temp, window_size^2, [])),...
                                    [], DIM);
        end
        u = u./(repmat(sum(u,2), [1, clusterNum]));
        
        % objective function
        dist4obj = zeros(N_DATA,clusterNum);
        for it_c = 1 : clusterNum
            temp = t_w(:).*dist(t_r(:),it_c);
            dist4obj(:,it_c) = reshape(sum(reshape(temp, window_size^2,...
                                [])), [], DIM);
        end
        obf = sum(sum(((u.^drg_m).*dist4obj)));
        if abs(obf-obf_prev) < epsilon
            break;
        end
        obf_prev = obf;
    end
 
end

%% Hierarchical FCM
function [u, c] = HFCM(obj, x, u, drg_m, drg_n, clusterNum, subClusterNum, max_iter, epsilon)
    DIM = size(x,2);
    N_DATA = size(x,1);
    v = rand(N_DATA, subClusterNum, clusterNum); % v_{i,k,j}
    v = v./repmat(sum(v,2),[1,subClusterNum,1]);
    c = zeros(clusterNum,DIM,subClusterNum);
    obf = 0; obf_prev = 0; halt = 0;
    for it = 1 : max_iter
        % update centroid
        for it_c = 1 : clusterNum
            for it_cc = 1 : subClusterNum
                c(it_c,:,it_cc) = sum(repmat(u(:,it_c).^drg_m.*...
                    v(:,it_cc, it_c).^drg_n, [1, DIM]).*x, 1);
                c(it_c,:,it_cc) = c(it_c,:,it_cc)/...
                    sum(u(:,it_c).^drg_m.*v(:,it_cc, it_c).^drg_n);
            end
        end
        % calculate distance 
        dist = zeros(N_DATA, subClusterNum, clusterNum);
        for it_c = 1 : clusterNum
            for it_cc = 1 : subClusterNum
                dist(:,it_cc,it_c) = pdist2(x,c(it_c,:,it_cc),'squaredeuclidean');
            end
        end
        % update membership
        for it_c = 1 : clusterNum
            u(:,it_c) = squeeze(sum(v(:,:,it_c).^drg_n.*...
                        dist(:,:,it_c),2)).^(1/(1-drg_m));
        end
        u = u./repmat(sum(u,2),[1,clusterNum]);       
        % update sub-membership
        for it_c = 1 : clusterNum
            v(:,:,it_c) = (repmat(u(:,it_c).^drg_m,[1, subClusterNum])...
                    .*dist(:,:,it_c)).^(1/(1-drg_n));
        end
        v = v./repmat(sum(v,2),[1, subClusterNum, 1]);
        
        % objective function
        obf = sum(sum(u.^drg_m.*squeeze(sum(v.^drg_n.*dist,2))));
        %fprintf('HFCM@Iteration count = %i, obj. fcn = %f \n', it, obf);
        if (obf-obf_prev)<epsilon
            halt = halt+1;
            if halt> 100
                break;
            end
        end
        obf_prev = obf;
    end
    fprintf('HFCM@Iteration count = %i, obj. fcn = %f \n', it, obf);
end

%% Generalized Hierarchical FCM
function [u, c] = GHFCM(obj, x, u, drg_m, drg_n, clusterNum, subClusterNum, max_iter, epsilon)
    N_DATA = size(x,1);
    DIM = size(x,2);
	window_size = 5;
    obf = 0; obf_prev = 0;
    v = rand(N_DATA, subClusterNum, clusterNum);
    v = v./repmat(sum(v,2),[1,subClusterNum,1]); % submembership 
    c = zeros(clusterNum,DIM,subClusterNum); % cluster centroid
    weight = zeros(N_DATA+1,window_size^2);
    row = uint32(zeros(N_DATA+1,window_size^2));
    x = [x; zeros(1, DIM)];

    parfor it_n = 1 : N_DATA
        [w, r] = spatialDist(obj, it_n, window_size);
        weight(it_n,:) = w;
        row(it_n,:) = r;
    end    
    weight = weight(1:end-1,:).';
    row = row(1:end-1,:).';
    
    % Iteration 
    for it = 1 : max_iter
        % update centroid
        for it_c = 1 : clusterNum
            for it_cc = 1 : subClusterNum
                c(it_c,:,it_cc) = reshape(sum(reshape(x(row(:),:), ... 
                        window_size^2,[])), [], DIM);
                c(it_c,:,it_cc) = c(it_c,:,it_cc).*sum(u(:,it_c).^drg_m ...
                        .*v(:,it_cc,it_c).^drg_n);
                c(it_c,:,it_cc) = c(it_c,:,it_cc)/...
                    sum(u(:,it_c).^drg_m.*v(:,it_cc, it_c).^drg_n);
            end
        end
        % calculate distance 
        dist = zeros(N_DATA+1, subClusterNum, clusterNum);
        for it_c = 1 : clusterNum
            for it_cc = 1 : subClusterNum
                dist(1:end-1,it_cc,it_c) = pdist2(x(1:end-1,:), ...
                        c(it_c,:,it_cc),'squaredeuclidean');
            end
        end
        % update membership (phase 1)
        for it_c = 1 : clusterNum
            u(:,it_c) = (v(:,:,it_c).^n.*sum(weight(:).*dist(:,:it_c) )).^(1/(1-drg_m));
        end
        u = bsxfun(@rdivide, u, sum(u,2));
        % update membership (phase 2)
        temp = u;
        for it_c = 1 : clusterNum
            temp(:,it_c) = 
        end
        u = temp;
        u = bsxfun(@rdivide, u, sum(u,2));
        % update submembership (phase 1)

        % update submembership (phase 2)

    end
end

%% Conventional FCM
function [u, c] = CFCM(obj, x, u, drg_m, clusterNum, max_iter, epsilon)
    obf = 0; obf_prev = 0;
    halt = 0;
    DIM = size(x,2);
    c = zeros(clusterNum, DIM); % cluster centroid

    for it = 1 : max_iter
		for it_c = 1 : clusterNum
			c(it_c,:) = sum(x.*(u(:,it_c).^drg_m), 1) ...
				./sum(u(:,it_c).^drg_m);
		end
        dist = pdist2(x,c,'squaredeuclidean');
        temp = dist.^(1/(1-drg_m));
        u = temp ./ repmat(sum(temp,2), [1, size(x, 2)]);
        u = u./sum(u,2);
        obf = sum(sum((u.^drg_m).*dist));
		fprintf('CFCM@Iteration count = %i, obj. fcn = %f \n', it, obf);
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