function varargout = myFCM (obj, x, clusterNum, max_iter, subClusterNum, algo)
    % MYFCM implement the fuzzy c-means clustering
    %
    % Syntax:
	%  * [label] = MYFCM(x, clusterNum, max_iter)
	%  * [label,c] = MYFCM(x, clusterNum, max_iter)
	%  * [label,c,u] = MYFCM(x, clusterNum, max_iter)
	%
    % Inputs:
    %  * x           : Data, specified as a numeric matrix. The rows of x 
    %                  correspond to observations, and the columns correspond
    %                  to variables.
    %  * clusterNum  : Number of clusters in the data, specified as a positive
    %                  integer.
    %  * subClusterNum: Number of sub-clusters of a cluster, specified as a positive integer.
    %  * max_iter    : Maximum number of iterations, specified as a positive integer.
    %
    % Outputs:
    %  * label  : Final label of each observation, returened as a matrix with N
    %             rows where N is the number of observation.
    %  * c      : Cluster centers, returned as a matrix with clusterNum rows containing 
    %             the coordinates of each cluster center.
    %  * u      : Memebership matrix, or fuzzy partition matrix, returned as a matrix 
    %             with N rows and clusterNum columns. 
    %
    % References:
    % [1] Image segmentation by generalized hierarchical fuzzy C-means algorithm
    %
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
    minArgs=1;  
    maxArgs=3;
    nargoutchk(minArgs,maxArgs)

    % Algorithm parameters
	N = size(x,1);
    rng(5)
	u = rand(N, clusterNum); % membership matrix
    u = u./sum(u,2);
    drg_m = 2; % drg_m is fuzzy partition matrix exponent for controlling
               % the degree of fuzzy overlap, with m > 1. Fuzzy overlap
	           % refers to how fuzzy the boundaries between clusters 
	           % are, that is the number of data points that have 
	           % significant membership in more than one cluster.
    drg_n = 2; 
    epsilon = 1e-5;

    % Execute fuzzy c-means clustering
    switch algo
        case 'GFCM'
            [u, c] = GFCM(obj, x, u, drg_m, clusterNum, max_iter, epsilon);
        case 'HFCM'
            [u, c] = HFCM(obj, x, u, drg_m, drg_n, clusterNum, subClusterNum, max_iter, epsilon);
        case 'GHFCM'
            [u, c] = GHFCM(obj, x, u, drg_m, drg_n, clusterNum, subClusterNum, max_iter, epsilon);
        case 'CFCM'
            [u, c] = CFCM(obj, x, u, drg_m, clusterNum, max_iter, epsilon);
        otherwise
            error('Incorrect algorithm name.')
    end

    % Output results
    [~, label] = max(u,[],2);
    if nargout>=1, varargout{1} = label;
        if nargout>=2, varargout{2} = c;
            if nargout>=3, varargout{3} = u;
            end
        end
    end
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
	window_size = 3;
    obf_prev = 0;
    halt = 3;
    v = rand(N_DATA, subClusterNum, clusterNum);
    v = v./repmat(sum(v,2),[1,subClusterNum,1]); % submembership 
    c = zeros(clusterNum,DIM,subClusterNum); % cluster centroid
    weight = zeros(window_size^2,N_DATA);
    row = uint32(zeros(window_size^2,N_DATA));
    x = [x; zeros(1, DIM)];

    parfor it_n = 1 : N_DATA
        [w, r] = spatialDist(obj, it_n, window_size);
        weight(:,it_n) = w;
        row(:,it_n) = r;
    end
    
    % Iteration 
    for it = 1 : max_iter
        % update centroid
        for it_c = 1 : clusterNum
            for it_cc = 1 : subClusterNum
                temp = squeeze(sum(reshape(x(row(:),:),window_size^2,[],DIM),1));
                c(it_c,:,it_cc) = sum(u(:,it_c).^drg_m ...
                        .*v(:,it_cc,it_c).^drg_n.*temp);
                c(it_c,:,it_cc) = c(it_c,:,it_cc)/...
                    sum(u(:,it_c).^drg_m.*v(:,it_cc, it_c).^drg_n);
            end
        end
        clear temp
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
            u(:,it_c) = sum((v(:,:,it_c).^drg_n.* ...
                reshape(sum(repmat(weight,[1,subClusterNum]).* ...
                    reshape(dist(row(:),:,it_c), window_size^2,[]),1),...
                    [],subClusterNum)),2).^(1/(1-drg_m));
        end
        u = bsxfun(@rdivide, u, sum(u,2));
        % update membership (phase 2)
        temp = [u; zeros(1,clusterNum)];
        for it_c = 1 : clusterNum
            u(:,it_c) = sum(weight.* ...
                reshape(temp(row(:),it_c), window_size^2, []),1);
        end
        u = bsxfun(@rdivide, u, sum(u,2));
        clear temp
        % update submembership (phase 1)
        for it_c = 1 : clusterNum
            for it_cc = 1 : subClusterNum
                v(:,it_cc,it_c) = (u(:,it_c).^drg_m.* ...
                    sum(reshape(dist(row(:),it_cc,it_c),window_size^2,[]).* ...
                    weight, 1).').^(1/(1-drg_n));
            end
        end
        v = v./repmat(sum(v,2),[1, subClusterNum, 1]);
        % update submembership (phase 2)
        temp = cat(1,v,zeros(1,subClusterNum,clusterNum));
        for it_c = 1 : clusterNum
            for it_cc = 1 : subClusterNum
                v(:, it_cc, it_c) = sum(weight.* ...
                    reshape(temp(row(:),it_cc, it_c), window_size^2, []),1);
            end
        end
        clear temp
        v = v./repmat(sum(v,2),[1, subClusterNum, 1]);
        % objective function 
        obf = 0;
        for it_c = 1 : clusterNum            
            obf = obf + sum(u(:,it_c).^drg_m.*sum(v(:,:,it_c).^drg_n.* ...
                reshape(sum(repmat(weight,[1,subClusterNum]).* ... 
                reshape(dist(row(:),:,it_c), window_size^2, []),1), ...
                [], subClusterNum),2));
        end
        %fprintf('GHFCM@Iteration count = %i, obj. fcn = %f \n', it, obf);
        if abs(obf-obf_prev) < epsilon
            halt = halt+1;
            if halt> 3
                break;
            end
        end
        obf_prev = obf;
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
        u = dist.^(1/(1-drg_m));
        u = u./repmat(sum(u,2),[1, clusterNum]);
        obf = sum(sum((u.^drg_m).*dist));
		%fprintf('CFCM@Iteration count = %i, obj. fcn = %f \n', it, obf);
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
	[row_, col_]= idx2nm(obj, n_1); 
	% obtain the index inside windows
	idx = [row_+Y(:), col_+X(:)]; % idx = (row, col)
    
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