function [H, alpha_bar] = eigenDecomposition(obj, varargin)
    % EIGENDECOMPOSITION implement eigen-decomposition.
    %
    % Syntax:  [H,alpha_bar] = eigenDecomp(NAME, VALUE)
    % Specifies function properties using one or more Name,Value pair arguments. 
    % Name-value pair settings apply to all the lines plotted. 
    %
    % Inputs:
    %    ('Calculate', c) - 
    %    ('Filename', fname) - 
    %    ('T', T) - 
    %    ('Plotsetting', plotf) -
    %
    % Outputs:
    %    H - Image of entropy of the PolSAR data.
    %    alpha_bar - Image of average alpha angle of the PolSAR data.
    %
    % Example: 
    %    Line 1 of example
    %    Line 2 of example
    %    Line 3 of example
    %
    % Other m-files required: none
    % Subfunctions: none
    % MAT-files required: none
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------
    
    %% H_alpha plot
    parse_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'logical'},{});
    validationFcn_2_ = @(x) validateattributes(x,{'char'},{'nonempty'});
    validationFcn_3_ = @(x) validateattributes(x,{'numeric'},{'3d'});
    validationFcn_4_ = @(x) validateattributes(x,{'function_handle'},{});
    addParameter(parse_,'Calculate',0,validationFcn_1_);
    addParameter(parse_,'Filename','',validationFcn_2_);
    addParameter(parse_,'T',[],validationFcn_3_);
    addParameter(parse_,'Plotsetting',@(x) 0,validationFcn_4_);
    parse(parse_,varargin{:})
    
    if parse_.Results.Calculate
        t = cputime;
        [~,~,num]= size(parse_.Results.T);
        T = parse_.Results.T;
        H = zeros(1, num); A_1 = zeros(1, num); A_2 = zeros(1, num); lambda = zeros(3,num);
        alpha_bar = zeros(1, num);
        seg = [1, num/4, num/2, 3*num/4, num];
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
            e = cputime -t;
            disp(['Execution time:' num2str(e)])
            %%
            H = single(real(reshape(H, obj.IMAGE_SIZE)));
            A_1 = single(reshape(A_1, obj.IMAGE_SIZE));
            A_2 = single(reshape(A_2, obj.IMAGE_SIZE));
            alpha_bar = single(reshape(alpha_bar, obj.IMAGE_SIZE));
            lambda = single(reshape(lambda.', [obj.IMAGE_SIZE 3]));
            disp('Saving results....')
            save([dir parse_.Results.Filename '.mat'],'-v7.3', 'H', 'A_1', 'A_2', 'lambda','alpha_bar');
    else 
        load([dir parse_.Results.Filename '.mat'])
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
end