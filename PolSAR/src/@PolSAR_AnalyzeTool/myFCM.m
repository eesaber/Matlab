function varargout = myFCM(obj, x, varargin)
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
    
    %% Parse input/output arguments
    % parse input arguments
    p_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'int32'},{'scalar'});
    validationFcn_2_ = @(x) validateattributes(x,{'numeric'},{'scalar'});
    validationFcn_3_ = @(x) validateattributes(x,{'char'},{'nonempty'});
    validationFcn_4_ = @(x) validateattributes(x,{'char','function_handle'},{'nonempty'});
    validationFcn_5_ = @(x) validateattributes(x,{'numeric'},{});
    addParameter(p_,'clusterNum', 4,validationFcn_1_);
    addParameter(p_,'subClusterNum', 0,validationFcn_1_);
    addParameter(p_,'max_iter', 100,validationFcn_1_);
    addParameter(p_,'drg_m', 2,validationFcn_2_);
    addParameter(p_,'drg_n', 2,validationFcn_2_);
    addParameter(p_,'algo', 'CFCM',validationFcn_3_);
    addParameter(p_,'distance', 'squaredeuclidean',validationFcn_4_);
    addParameter(p_,'covariance_matirx', [],validationFcn_5_);
    parse(p_,varargin{:})
    p_ = p_.Results;
    if strcmp(p_.distance, 'wishart')
        p_.distance = @wishart;
    end
    % parse output arguments
    minArgs=1;
    maxArgs=3;
    nargoutchk(minArgs,maxArgs)
    
    %% Algorithm parameters
    % drg_m is fuzzy partition matrix exponent for controlling
    % the degree of fuzzy overlap, with m > 1. Fuzzy overlap
    % refers to how fuzzy the boundaries between clusters 
    % are, that is the number of data points that have 
    % significant membership in more than one cluster.
    N = size(x,1);
    %rng(3)
    rng(2)
    u = rand(N, p_.clusterNum); % membership matrix
    u = u./sum(u,2);
    epsilon = 1e-10;

    %% Execute fuzzy c-means clustering
    switch p_.algo
        case 'GFCM'
            [u, c] = GFCM(obj, x, u, p_.drg_m, p_.clusterNum, p_.max_iter, epsilon, p_.distance, p_.covariance_matirx);
        case 'HFCM'
            [u, c] = HFCM(obj, x, u, p_.drg_m, p_.drg_n, p_.clusterNum, p_.subClusterNum, p_.max_iter, epsilon, p_.distance);
        case 'GHFCM'
            [u, c] = GHFCM(obj, x, u, p_.drg_m, p_.drg_n, p_.clusterNum, p_.subClusterNum, p_.max_iter, epsilon, p_.distance);
        case 'CFCM'
            [u, c] = CFCM(obj, x, u, p_.drg_m, p_.clusterNum, p_.max_iter, epsilon, p_.distance, p_.covariance_matirx);
        otherwise
            error('Incorrect algorithm name.')
    end
    disp('Done')
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
function [u, c] = GFCM(obj, x, u, drg_m, clusterNum, max_iter, epsilon, distance_metric, covariance_matirx)
    N_DATA = size(x,1);
    DIM = size(x,2);
    window_size = 3;
    obf = 0; obf_prev = 0; halt = 0;
    c = zeros(clusterNum, DIM); % cluster centroid
    c_cov = zeros(3,3,clusterNum); % cluster centroid in covariance 
    weight = zeros(window_size^2,N_DATA);
    row = uint32(zeros(window_size^2,N_DATA));
    x = [x; zeros(1, DIM)];
    rebuildCov = @(x) [x(1), sqrt(2)*(x(4)+1j*x(5)), x(6)+1j*x(7);...
        sqrt(2)*(x(4)-1j*x(5)), 2*x(2), sqrt(2)*(x(8)+1j*x(9));...
        x(6)-1j*x(7), sqrt(2)*(x(8)-1j*x(9)), x(3)];
    parfor it_n = 1 : N_DATA
        [w, r] = spatialDist(obj, it_n, window_size);
        weight(:,it_n) = w;
        row(:,it_n) = r;
    end    
    
    for it = 1 : max_iter
        % update centroid
        
        for it_c = 1 : clusterNum
            temp = sum(squeeze(sum(reshape(x(row(:),:).* ...
                reshape(weight,[],1),window_size^2,[],DIM),1)).* ...
                u(:,it_c).^drg_m ,1);
            c(it_c,:) = temp./sum(u(:,it_c).^drg_m.* ...
                sum(weight,1).',1);
            if isa(distance_metric, 'function_handle')
                c_cov(:,:,it_c) = rebuildCov(c(it_c,:));
            end
        end
        if isa(distance_metric, 'function_handle')  
            dist = [wishart(covariance_matirx, c_cov); zeros(1,clusterNum)];
        else
            dist = [pdist2(x,c, distance_metric); zeros(1,clusterNum)];
        end
        % dist(i,j) is the distance between observation i in X 
        % and observation j in Y.
        % update membership (phase 1)
        for it_c = 1 : clusterNum
           u(:,it_c) = sum(weight.* ...
               reshape(dist(row(:),it_c), window_size^2, []),1).^(1/(1-drg_m));
        end
        u = u./(repmat(sum(u,2), [1, clusterNum]));
        % update membership (phase 2)
        temp = [u; zeros(1,clusterNum)];
        for it_c = 1 : clusterNum
            u(:,it_c) = sum(reshape(temp(row(:),it_c),window_size^2,[]).* ...
                weight,1);
        end
        u = u./(repmat(sum(u,2), [1, clusterNum]));
        clear temp        
        % objective function
        obf = 0;
        for it_c = 1 : clusterNum
            obf = obf + sum(u(:,it_c).^drg_m.*sum(weight.* ...
                reshape(dist(row(:),it_c),window_size^2,[]),1).');
        end
        if abs(obf-obf_prev) < epsilon
            halt = halt + 1;
            if halt>3, break; end
        end
        obf_prev = obf;
    end
    fprintf('GFCM@Iteration count = %i, objective function = %f \n', it, obf);
end

%% Hierarchical FCM
function [u, c] = HFCM(obj, x, u, drg_m, drg_n, clusterNum, subClusterNum, max_iter, epsilon, distance_metric)
    DIM = size(x,2);
    N_DATA = size(x,1);
    v = rand(N_DATA, subClusterNum, clusterNum); % v_{i,k,j}
    v = v./repmat(sum(v,2),[1,subClusterNum,1]);
    c = zeros(subClusterNum, DIM, clusterNum);
    obf = 0; obf_prev = 0; halt = 0;
    for it = 1 : max_iter
        % update centroid
        for it_c = 1 : clusterNum
            for it_cc = 1 : subClusterNum
                c(it_cc,:,it_c) = sum(repmat(u(:,it_c).^drg_m.*...
                    v(:,it_cc, it_c).^drg_n, [1, DIM]).*x, 1);
                c(it_cc,:,it_c) = c(it_cc,:,it_c)/...
                    sum(u(:,it_c).^drg_m.*v(:,it_cc, it_c).^drg_n);
            end
        end
        % calculate distance 
        dist = zeros(N_DATA, subClusterNum, clusterNum);
        for it_c = 1 : clusterNum
            for it_cc = 1 : subClusterNum
                dist(:,it_cc,it_c) = pdist2(x,c(it_cc,:,it_c), distance_metric);
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
            if halt> 3, break; end
        end
        obf_prev = obf;
    end
    fprintf('HFCM@Iteration count = %i, objective function = %f \n', it, obf);
end

%% Generalized Hierarchical FCM
function [u, c] = GHFCM(obj, x, u, drg_m, drg_n, clusterNum, subClusterNum, max_iter, epsilon, distance_metric)
    N_DATA = size(x,1);
    DIM = size(x,2);
    window_size = 3;
    obf = 0; obf_prev = 0; halt = 0;
    v = rand(N_DATA, subClusterNum, clusterNum);
    v = v./repmat(sum(v,2),[1,subClusterNum,1]); % submembership 
    c = zeros(subClusterNum,DIM,clusterNum); % cluster centroid
    weight = zeros(window_size^2,N_DATA);
    row = uint32(zeros(window_size^2,N_DATA));
    x = [x; zeros(1, DIM)];

    parfor it_n = 1 : N_DATA
        [w, r] = spatialDist(obj, it_n, window_size);
        weight(:,it_n) = w;
        row(:,it_n) = r;
    end
    figure
    h = animatedline;
    xlabel('iter time')
    ylabel('objective function')
    % Iteration 
    for it = 1 : max_iter
        % update centroid
        for it_c = 1 : clusterNum
            for it_cc = 1 : subClusterNum
                % original                 
                c(it_cc,:,it_c) = sum(repmat(u(:,it_c).^drg_m.* ...
                    v(:,it_cc,it_c).^drg_n, [1, DIM]).*...
                    squeeze(sum(reshape(x(row(:),:),window_size^2,[],DIM),1))...
                    ./repmat(sum(weight~=0,1).',[1,DIM]),1);
                c(it_cc,:,it_c) = c(it_cc,:,it_c)/...
                    sum(u(:,it_c).^drg_m.*v(:,it_cc, it_c).^drg_n);
                %{
                % modified (fail)
                temp = squeeze(sum(reshape(x(row(:),:).* ...
                    reshape(weight,[],1),window_size^2,[],DIM),1));
                temp = temp./repmat(sum(weight,1).',[1, DIM]);
                temp = sum(temp.*repmat(u(:,it_c).^drg_m.* ...
                    v(:,it_cc,it_c).^drg_n,[1 DIM]),1);
                c(it_cc,:,it_c) = temp./sum(u(:,it_c).^drg_m.* ...
                    v(:,it_cc,it_c));
                %}
            end
        end
        clear temp
        % calculate distance 
        dis = zeros(N_DATA+1, subClusterNum, clusterNum);
        for it_c = 1 : clusterNum
            for it_cc = 1 : subClusterNum
                dis(1:end-1,it_cc,it_c) = pdist2(x(1:end-1,:), ...
                    c(it_cc,:,it_c), distance_metric);
            end
        end
        
        % update membership (phase 1)
        for it_c = 1 : clusterNum
            %{
            mem = 0;
            for it_cc = 1 : subClusterNum 
                mem = mem + v(:,it_cc,it_c).^drg_n.*sum(weight.* ...
                    reshape(dis(row(:),it_cc,it_c),window_size^2,[]),1).';
            end
            %}
            u(:,it_c) = sum((v(:,:,it_c).^drg_n.* ...
                reshape(sum(repmat(weight,[1,subClusterNum]).* ...
                    reshape(dis(row(:),:,it_c), window_size^2,[]),1),...
                    [],subClusterNum)),2).^(1/(1-drg_m));
            
            %u(:,it_c) = mem.^(1/(1-drg_m));
        end
        u = u./repmat(sum(u,2),[1,clusterNum]);
        % update membership (phase 2)
        temp = [u; zeros(1,clusterNum)];
        for it_c = 1 : clusterNum
            u(:,it_c) = sum(weight.* ...
                reshape(temp(row(:),it_c), window_size^2, []),1);
        end
        u = u./repmat(sum(u,2),[1,clusterNum]);
        clear temp
        % update submembership (phase 1)
        for it_c = 1 : clusterNum
            for it_cc = 1 : subClusterNum
                v(:,it_cc,it_c) = (u(:,it_c).^drg_m.* ...
                    sum(reshape(dis(row(:),it_cc,it_c),window_size^2,[]).* ...
                    weight, 1).').^(1/(1-drg_n));
                % v(:,:,it_c) = (repmat(u(:,it_c).^drg_m,[1, subClusterNum])...
                    % .*dis(1:end-1,:,it_c)).^(1/(1-drg_n));    
            end
        end
        v = v./repmat(sum(v,2),[1, subClusterNum, 1]);
        % update submembership (phase 2)
        temp = cat(1,v,zeros(1,subClusterNum,clusterNum));
        for it_c = 1 : clusterNum
            for it_cc = 1 : subClusterNum
                v(:,it_cc,it_c) = sum(weight.* ...
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
                reshape(dis(row(:),:,it_c), window_size^2, []),1), ...
                [], subClusterNum),2));
        end
        addpoints(h,it,obf);
        drawnow
        if abs(obf-obf_prev) < epsilon
            halt = halt+1;
            if halt> 3, break; end
        end
        obf_prev = obf;
        if mod(it,100) == 0, disp(it); end
    end
    fprintf('GHFCM@Iteration count = %i, objective function = %f \n', it, obf);
end

%% Conventional FCM
function [u, c] = CFCM(obj, x, u, drg_m, clusterNum, max_iter, epsilon, distance_metric, covariance_matirx)
    obf = 0; obf_prev = 0;
    halt = 0;
    DIM = size(x,2);
    c = zeros(clusterNum, DIM); % cluster centroid
    c_cov = zeros(3,3,clusterNum); % cluster centroid in covariance 
    rebuildCov = @(x) [x(1), sqrt(2)*(x(4)+1j*x(5)), x(6)+1j*x(7);...
        sqrt(2)*(x(4)-1j*x(5)), 2*x(2), sqrt(2)*(x(8)+1j*x(9));...
        x(6)-1j*x(7), sqrt(2)*(x(8)-1j*x(9)), x(3)];
        
    for it = 1 : max_iter
        for it_c = 1 : clusterNum
            c(it_c,:) = sum(x.*(u(:,it_c).^drg_m), 1) ...
                ./sum(u(:,it_c).^drg_m);
            if isa(distance_metric, 'function_handle')
                c_cov(:,:,it_c) = rebuildCov(c(it_c,:));
            end
        end
        if isa(distance_metric, 'function_handle')
            dist = wishart(covariance_matirx, c_cov);
            dist = dist + 2*min(min(dist));
        else
            dist = pdist2(x,c, distance_metric);
        end
        u = dist.^(1/(1-drg_m));
        u = u./repmat(sum(u,2),[1, clusterNum]);
        obf = sum(sum((u.^drg_m).*dist));
        %fprintf('CFCM@Iteration count = %i, obj. fcn = %f \n', it, obf);
        if abs(obf-obf_prev) < epsilon
            %halt = halt+1;
            if halt> 3
                break;
            end
        end
        obf_prev = obf;
    end
    fprintf('CFCM@Iteration count = %i, objective function = %f \n', it, obf);
end

%% ------------- Utility function -------------- %%
function [weight, row_n] = spatialDist(obj, n_1, window_size)
    % SPATIALDIST 
    %
    % Inputs:
    %	* n_1 : linear index of data
    %	* window_size : the windows size.
    %
    % Outputs:
    %	* weight : the distance weight array with size 1xN  
    %   * row_n  : A array with size 1xN which containing
    %              the row number of the pixel in the window.
    %------------- <<<<< >>>>>--------------
    
    sigma = (window_size-1)/2; % Spatial relation
    box = -(window_size-1)/2:1:(window_size-1)/2;
    [X,Y] = meshgrid(box,box);
    % obtain the index of the center of the windows.
    [row_, col_]= ind2sub(obj.IMAGE_SIZE, n_1);
    % obtain the index inside windows
    idx = [row_+Y(:), col_+X(:)]; % idx = (row, col)
    % exclude the elements which is not in the matrix
    idx = idx(( (idx(:,1)<=0)+(idx(:,1)>obj.IMAGE_SIZE(1))+...
                (idx(:,2)<=0)+(idx(:,2)>obj.IMAGE_SIZE(2)) )==0, :);   
    % obtain the distance between the center and the elements in windows
    dist = pdist2([row_, col_], idx);
    padding = zeros(1,window_size^2-numel(dist));
    % obtain the weight of each element
    weight = [1/sqrt(2*pi*sigma^2)*exp(-dist/(2*sigma^2)), ...
            padding];
    % obtain the row number of each elements
    row_n = uint32([sub2ind(obj.IMAGE_SIZE, idx(:,1), idx(:,2)).' ...
         obj.IMAGE_SIZE(1)*obj.IMAGE_SIZE(2)+1+padding]);
end
function D2 = wishart(x, c)
    % WISHART implement Wishart distance
    %
    % Inputs:
    %	* x is a (3,3,n) ndarray containing the observations
    %	* c is a (3,3,m) ndarray containing the centroids
    %	* D2 is a (n,m) matrix where D2(i,j) represents the distance between 
    %       ith observation and mth center.
    %
    % Definition:
    %   * D(a, b) = ln(|b|) + Tr(b^{-1}*a)
    %     where a is observation and b is cluster centroid.
    %
    % Reference:
    % [1] Unsupervised classification using polarimetric decomposition
    %     and the complex Wishart classifier.
    %
    %------------- <<<<< >>>>>--------------
    
    D2 = single(zeros(size(x,3), size(c,3)));
    for it_c = 1 : size(c,3)
        ln_det_c = real(log(det(c(:,:,it_c))));
        inv_c = inv(c(:,:,it_c));
        temp = repmat(inv_c,1,1,size(D2,1)).*x;
        temp = temp(1,1,:)+temp(2,2,:)+temp(3,3,:);
        D2(:,it_c) = ln_det_c+squeeze(temp);
    end
    D2 = D2+2*min(min(D2));
end