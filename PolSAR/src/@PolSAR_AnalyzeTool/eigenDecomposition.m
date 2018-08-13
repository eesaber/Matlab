function [H, alpha_bar] = eigenDecomposition(obj, Calculate, Filename, varargin)
    % EIGENDECOMPOSITION implement eigen-decomposition.
    %
    % Syntax:
	%	[H,alpha_bar] = EIGENDECOMPOSITION(Calculate, Filename)
	%
    % Description:
    %   EIGENDECOMPOSITION(Calculate, Filename)
    %	* If Calculate is 1, the entropy data will be computed and saved
    %	with the Filename.
    %	* If Calculate is 0, the entropy of data will be loaded from the
    %	file Filename.
    % Name-Value Pair Arguments:
    %   * SaveResult - logical
    %     If Calculate is 1, then the calculated H, alpha_bar, ... etc 
    %     will be saved as <Filename>.m 
    %     If Calculate is 0, it will not do anything.
    %
    % Outputs:
    %    H - Image of entropy of the PolSAR data.
    %    alpha_bar - Image of average alpha angle of the PolSAR data.
    %
    % Other m-files required: none
    % Subfunctions: none
    % MAT-files required: none
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------
    parse_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'logical'},{});
    addParameter(parse_, 'SaveResults', 1, validationFcn_1_);
    
    parse(parse_,varargin{:})

    num= obj.IMAGE_SIZE(1)*obj.IMAGE_SIZE(2);

    if Calculate
		T = cat(1,cat(2, reshape(obj.T_11,[1,1,num]), reshape(obj.T_12,[1,1,num]), reshape(obj.T_13,[1,1,num])), ...
             cat(2,reshape(conj(obj.T_12),[1,1,num]), reshape(obj.T_22,[1,1,num]), reshape(obj.T_23,[1,1,num])),....
             cat(2,reshape(conj(obj.T_13),[1,1,num]), reshape(conj(obj.T_23),[1,1,num]), reshape(obj.T_33,[1,1,num])));
		if obj.IS_BIGFILE
			setCov2Zero(obj)
            setCoh2Zero(obj)
		end
        H = zeros(1, num); A_1 = zeros(1, num); A_2 = zeros(1, num); 
		lambda = zeros(3,num); alpha_bar = zeros(1, num);
        seg = [1, num/4, num/2, 3*num/4, num];
		disp('Calculating....')
        for b = 1 : 4
            parfor r = ceil(seg(b)) : floor(seg(b+1))
                if sum(sum(T(:,:,r)== zeros(3,3))) == 9
                    continue;
                end
                [U, L] = eig(T(:,:,r));
                P = diag(L)/sum(diag(L));
                H(r) = -sum(P.*log(P)/log(3));
                A_temp = sort(diag(L),'descend');
                A_1(r) = (A_temp(2)-A_temp(3))/(A_temp(2)+A_temp(3));
                A_2(r) = (A_temp(1)-A_temp(2))/(A_temp(1)+A_temp(2));
                lambda(:,r) = A_temp;
                alpha_bar(r) = sum(P'.*acosd(abs(U(1,:))));
            end
        end
		H = single(real(reshape(H, obj.IMAGE_SIZE)));
		A_1 = single(reshape(A_1, obj.IMAGE_SIZE));
		A_2 = single(reshape(A_2, obj.IMAGE_SIZE));
		alpha_bar = single(reshape(alpha_bar, obj.IMAGE_SIZE));
        lambda = single(reshape(lambda.', [obj.IMAGE_SIZE 3]));
        if parse_.Results.SaveResults
		    disp('Saving results....')
            save([obj.INPUT_PATH Filename '.mat'],'-v7.3', 'H', 'A_1', 'A_2', 'lambda','alpha_bar');
        end
    else 
        load([obj.INPUT_PATH Filename '.mat'])
    end
    figure
        imagesc(H)
        obj.plotSetting([0 1])
        plot_para('Filename','Entropy', 'Maximize',true)
    figure
        imagesc(10*log10(abs(lambda(:,:,1))))
        obj.plotSetting([-35 -5],'Colorbar_unit',"(dB)")
        plot_para('Filename','lambda_1', 'Maximize',true)
    figure
        imagesc(alpha_bar)
        obj.plotSetting([0 60],'Colorbar_unit',"(deg)")
        plot_para('Filename','alpha', 'Maximize',true)
    figure
        imagesc(A_1)
        obj.plotSetting([0 0.5])
        plot_para('Filename','Anisotropy1', 'Maximize',true)
    figure
        imagesc(A_2)
        obj.plotSetting([0.6 1])
        plot_para('Filename','Anisotropy2', 'Maximize',true) 
	if obj.IS_BIGFILE
		clear A_1 A_2 lambda
		obj.readPolsarData();
	end
end