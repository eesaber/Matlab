function varargout = mySVM(obj, x, y)
    % MYSVM use LIBSVM [1] to do the PolSAR classification
    %
    % Syntax:
	%  * 
	%
    % Inputs:
    %  * x
    %
    % Outputs:
    %  * 
    %
    % References:
    % [1] https://www.csie.ntu.edu.tw/~cjlin/libsvm/
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------
    
    dim = size(x,3);
    % Split Data
    train_data = reshape(x(:,1501:end,:),[],dim);
    train_label = reshape(y(:,1501:end),[],1);
    test_data = reshape(x(:,1:1500,:),[],dim);
    test_label =  reshape(y(:,1:1500),[],1);
    % sampling
    rng(1)
    p = randperm(size(train_data,1),80000);
    
    % Radial Kernel
    disp('SVM training...')
    tic
    model_linear = svmtrain(train_label(p), train_data(p,:), '-s 0 -t 2 -g 10 -c 100 -e 0.1');
    [predict_label_L, accuracy_L, dec_values_L] = ...
    svmpredict(test_label, test_data, model_linear);
    toc
    fprintf('Test case accuracy: %f \n', accuracy_L) % Display the accuracy using linear kernel
    % Output results
    disp('SVM predicting...')
    label = svmpredict([test_label; train_label], [test_data; train_data], model_linear);
    disp('Session over..')
    if nargout>=1
        varargout{1} = label;
    end

end