function varargout = myKmeans(obj, x , varagin)
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
    addParameter(p_,'covariance_matirx', [],validationFcn_3_);
    parse(p_,varargin{:})
    p_ = p_.Results;
    minArgs=0;
    maxArgs=1;
    nargoutchk(minArgs,maxArgs)

    %% k-means
    if strcmp(p_distance, 'wishart')
        label = kmeans4Wishart(x, p_);
    else
        label = kmeans(im, p_.maxClusterNum,...
            'distance', p_.distance,...
            'Replicates', p_replicates,...
            'MaxIter', p_.max_iter,...
            'Options',statset('Display','final'));
    end
    %% Output arguments
    if nargout>=1, varargout{1} = label; end
end
function label = kmeans4Wishart(x, p_)
    label_sum = zeros(size(x,1));
    for 1 : p_.replicates 
        [label, c] = kmeans(x, p_.maxClusterNum,...
                        'distance', p_.distance,...
                        'Replicates', 1,...
                        'MaxIter', p_.max_iter,...
                        'Options',statset('Display','final'));
        counter = 1; 
        label = init(X, m);
        n = numel(label);
        idx = 1:n;
        last = zeros(1,n);
        while any(label ~= last) || counter>p_.max_iter
            [~,~,last(:)] = unique(label);                  % remove empty clusters
            mu = X*normalize(sparse(idx,last,1),1);         % compute cluster centers 
            [val,label] = min(dot(mu,mu,1)'/2-mu'*X,[],1);  % assign sample labels
            counter = counter+1;
        end
        label_sum = label_sum + label;
    end
end