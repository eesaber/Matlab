function varargout = myKmeans(obj, x , varargin)
    % MYKMEANS encapsulation k-means with custom distance metric
    %
    % Syntax:
    %  * [label] = MYFCM(x)
    %  * [label] = MYFCM(x, name, value)
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
    % Reference
    %  [1] K-means++: The Advantages of Careful Seeding.
    % 
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------

    %% Parse input/output arguments
    % parse input arguments
    p_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'int32'},{'scalar'});
    validationFcn_2_ = @(x) validateattributes(x,{'char','function_handle'},{'nonempty'});
    validationFcn_3_ = @(x) validateattributes(x,{'numeric'},{});
    addParameter(p_,'maxClusterNum', 10,validationFcn_1_);
    addParameter(p_,'max_iter', 100,validationFcn_1_);
    addParameter(p_,'replicates', 3,validationFcn_1_);
    addParameter(p_,'distance', 'squaredeuclidean',validationFcn_2_);
    addParameter(p_,'covariance', [],validationFcn_3_);
    parse(p_,varargin{:})
    p_ = p_.Results;
    minArgs=0;
    maxArgs=2;
    nargoutchk(minArgs,maxArgs)
    rng(10)
    %% k-means
    tic 
    if strcmp(p_.distance, 'wishart')
        [label, centroid] = kmeans4Wishart(obj, x, p_);
    else
        [label, centroid] = kmeans(im, p_.maxClusterNum,...
            'distance', p_.distance,...
            'Replicates', p_replicates,...
            'MaxIter', p_.max_iter,...
            'Options',statset('Display','final'));
    end
    toc
    %% Output arguments
    if nargout>=1, varargout{1} = uint8(label);
        if nargout>=2, varargout{2} = centroid;
        end
    end
end

function [label_ans, centroid] = kmeans4Wishart(obj, x, p_)
    label_ans = zeros(size(x,1),1);
    dist_prev = inf;
    rebuildCov = @(x) [x(1), sqrt(2)*(x(4)+1j*x(5)), x(6)+1j*x(7);...
        sqrt(2)*(x(4)-1j*x(5)), 2*x(2), sqrt(2)*(x(8)+1j*x(9));...
        x(6)-1j*x(7), sqrt(2)*(x(8)-1j*x(9)), x(3)];
    for it = 1 : p_.replicates 
        counter = 0; 
        %{
        centroid = single(zeros(3, 3, p_.maxClusterNum));
        [label, c] = kmeans(x, p_.maxClusterNum,...
                        'distance', 'sqeuclidean',...
                        'Replicates', 1,...
                        'MaxIter', 3,...
                        'Options',statset('Display','final'));
        label = single(label);
        for it_c = 1 : p_.maxClusterNum
           centroid(:,:,it_c) = rebuildCov(single(c(it_c,:)));
        end
        %}
        centroid = init(obj, p_.covariance, p_.maxClusterNum);
        num_pixel = size(x,3);
        label = single(ones(num_pixel));
        last = single(zeros(num_pixel));
        %% iteration 
        %textprogressbar('iteration: ',0);
        while any(label ~= last) && counter<p_.max_iter
            %textprogressbar(counter,p_.max_iter);
            %pause(0.1);
            last = label;
            [label, centroid, dist_sum] = find_label(obj, p_.covariance, centroid);
            counter = counter+1;
            fprintf('%i,', counter)
        end
        %textprogressbar('done',0);
        disp(' ')
        if dist_prev > dist_sum
            label_ans =  label;
        end
    end
end

function [new_label, centroid, dist_sum] = find_label(obj, x, c)
    % calculate revised Wishart distance (KL divergence)
    D2 = obj.wishartDist(x, c);
    % assign label
    [~, label]= min(D2,[],2);
    [avaliable_idx,~,new_label(:)] = unique(label); % remove empty cluster
    num_centroid = numel(avaliable_idx); % get number of cluster
    dist_sum = sum(sum(D2)); % sum of total distance
    centroid = single(zeros(3,3,num_centroid));
    for it_c = 1 : num_centroid
        centroid(:,:,it_c) = sum(x(:,:,new_label==it_c),3)./sum(new_label==it_c);
    end
end
function [centroid] = init(obj, x, num_c)
    disp('initializing...')
    tic
    centroid = single(zeros(3,3,num_c));
    num_c = single(num_c);
    
    if isempty(obj.H)
        error('Please calcuilate polarimetry entropy first.')
    end
    H = obj.H(:);
    %{
    %alpha = obj.alpha_bar(:);
    %{
    h = histogram2(x_c.H,x_c.alpha_bar,'YBinLimits',[0 90],'XBinLimits',[0 1],'NumBins',[450 1350]);
    [~, idx] = max(h.Values(:));
    [idx_row, idx_col] = ind2sub(size(h.Values), idx);
    fprintf('H in [%f,%f]\n', idx_row*(1-0)/size(h.Values,1), (idx_row+1)*(1-0)/size(h.Values,1))
    fprintf('alpha in [%f,%f]\n', idx_col*(90-0)/size(h.Values,2), (idx_col+1)*(90-0)/size(h.Values,2))
    range = (x_c.H(:)>=idx_row*(1-0)/size(h.Values,1)).*(x_c.H(:)<=(idx_row+1)*(1-0)/size(h.Values,1));
    range = range.*(x_c.alpha_bar(:)>=idx_col*(90-0)/size(h.Values,2)).*(x_c.alpha_bar(:)<= (idx_col+1)*(90-0)/size(h.Values,2));
    range = logical(range);
    %}
    %}
    global kpp
    if ~kpp 
        section = 1/num_c:1/num_c:1;
        for it = 1 : numel(section)
            tmp = x(:,:,H<0.5);
            centroid(:,:,it) = tmp(:,:,randi(size(tmp,3),1));
        end
    else 
        %% k-means++ algorithm [1]
        disp('k-means++ initialization used.')
        centroid(:,:,1) = x(:,:,randi(size(x,3))); % take one center uniformly.
        dist_ = single(zeros(size(x,3), num_c-1)); 
        % take the new center with probability D(x)
        for it = 2 : size(centroid,3)
            ln_det_c = log10(real(det(centroid(:,:,it-1)))); % b
            inv_c = inv(centroid(:,:,it-1));
            for it_x = 1 : size(x,3)
                ln_det_x = log10(real(det(x(:,:,it_x)))); % a 
                dist_(it_x,it-1) = ln_det_c - ln_det_x + real(trace(inv_c*x(:,:,it_x)));
            end
            prob = min(dist_(:,1:it-1),[],2).^2;
            centroid(:,:,it) = x(:,:,randsample(size(x,3),1,true,prob));
        end
        
    end
    toc
        disp('initializing finished...')
        disp(centroid)
end

%{
function [label, mu, energy] = kmeans(X, m)
% Perform kmeans clustering.
% Input:
%   X: d x n data matrix
%   m: initialization parameter
% Output:
%   label: 1 x n sample labels
%   mu: d x k center of clusters
%   energy: optimization target value
% Written by Mo Chen (sth4nth@gmail.com).
label = init(X, m);
n = numel(label);
idx = 1:n;
last = zeros(1,n);
while any(label ~= last)
    [~,~,last(:)] = unique(label);                  % remove empty clusters
    mu = X*normalize(sparse(idx,last,1),1);         % compute cluster centers 
    [val,label] = min(dot(mu,mu,1)'/2-mu'*X,[],1);  % assign sample labels
end
energy = dot(X(:),X(:),1)+2*sum(val);
function label = init(X, m)
[d,n] = size(X);
if numel(m) == 1                           % random initialization
    mu = X(:,randperm(n,m));
    [~,label] = min(dot(mu,mu,1)'/2-mu'*X,[],1); 
elseif all(size(m) == [1,n])               % init with labels
    label = m;
elseif size(m,1) == d                      % init with seeds (centers)
    [~,label] = min(dot(m,m,1)'/2-m'*X,[],1); 
end
end
function Y = normalize(X, dim)
    % Normalize the vectors to be summing to one
    %   By default dim = 1 (columns).
    % Written by Michael Chen (sth4nth@gmail.com).
    if nargin == 1
        % Determine which dimension sum will use
        dim = find(size(X)~=1,1);
        if isempty(dim), dim = 1; end
    end
    Y = X./sum(X,dim);
end
%}