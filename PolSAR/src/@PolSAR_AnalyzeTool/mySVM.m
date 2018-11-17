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
    
    % Split Data
    rng(1)
    dim = size(x,3);
    p = randperm(numel(y));
    data_spilt = int32(numel(y)*0.9);
    x = reshape(x,[],dim);
    y = y(:);
    train_data = x(1:data_spilt,:);
    train_label = y(1:data_spilt);
    test_data = x(data_spilt+1:end, :);
    test_label =  y(data_spilt+1:end);
    
    % SVM
    disp('SVM training...')
    tic
    model_linear = svmtrain(train_label, train_data, '-s 0 -t 2 -g 1 -c 10 -e 0.1 -m 8000');
    [predict_label_L, accuracy_L, dec_values_L] = svmpredict(test_label, test_data, model_linear);
    toc
    fprintf('Test case accuracy: %f \n', accuracy_L) % Display the accuracy using linear kernel
    % Output results
    %yes_predict = input('Do you want to predict? [1/0] \n');
    yes_predict = 1;
    if yes_predict
        disp('SVM predicting...')
        label = svmpredict(y, x, model_linear);
    else
        label = [];
    end
    disp('Session over..')
    if nargout>=1, varargout{1} = label; 
        if nargout >= 2,  varargout{2} = model_linear;
        end
    end
end