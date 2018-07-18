function [H, lambda] = Eigen_decomp(varargin)
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
    
    global dir im_size
    if parse_.Results.Calculate
        %ind_ = randperm(im_size(1)*im_size(2));
        %figure
        %    H_Alpha(parse_.Results.T(:,:,ind_(1:2000)))
        %clear ind_
        % Plot each entropy of each pixel
        t = cputime;
        [~,~,num]= size(parse_.Results.T);
        T = parse_.Results.T;
        H = zeros(1, num); A_1 = zeros(1, num); A_2 = zeros(1, num); lambda = zeros(3,num);
        alpha_bar = zeros(1, num);
        
        parfor r = 1 : num
            %if sum(sum(parse_.Results.T(:,:,r)== zeros(3,3))) == 9
            %    continue;
            %end
            [U, L] = eig(T(:,:,r));
            P = diag(L)/sum(diag(L));
            H(r) = -sum(P.*log(P)/log(3));
            A_temp = sort(diag(L),'descend');
            A_1(r) = (A_temp(2)-A_temp(3))/(A_temp(2)+A_temp(3));
            A_2(r) = (A_temp(1)-A_temp(2))/(A_temp(1)+A_temp(2));
            lambda(:,r) = A_temp;
            alpha_bar(r) = sum(P'.*acosd(abs(U(1,:))));
        end
        e = cputime -t;
        disp(['Execution time:' num2str(e)])
        %%
        H = single(real(reshape(H, im_size)));
        A_1 = single(reshape(A_1, im_size));
        A_2 = single(reshape(A_2, im_size));
        alpha_bar = single(reshape(alpha_bar, im_size));
        lambda = single(reshape(lambda, [im_size 3]));
        disp('Saving results....')
        save([dir parse_.Results.Filename '.mat'],'-v7.3', 'H', 'A_1', 'A_2', 'lambda','alpha_bar');
    else 
        load([dir parse_.Results.Filename '.mat'])
    end
    %%
    x = [1900, 4200, 5100];
    y = [3100, 1600, 500];
    figure
        imagesc(H)
        parse_.Results.Plotsetting([0 1],1)
        GOM2()
        plot_para('Filename','Entropy_am', 'Maximize',true)
    %%
    figure
        imagesc(10*log10(abs(lambda(:,:,1))))
        parse_.Results.Plotsetting([-35 -5],1,'Colorbar_unit',"(dB)")
        GOM2()
        plot_para('Filename','lambda_1', 'Maximize',true)
    figure
        imagesc(alpha_bar)
        parse_.Results.Plotsetting([0 60],1,'Colorbar_unit',"(deg)")
        GOM2()
        plot_para('Filename','alpha', 'Maximize',true)
    figure
        imagesc(A_1)
        parse_.Results.Plotsetting([0 0.5],1)
        GOM2()
        plot_para('Filename','Anisotropy1', 'Maximize',true)
    figure
        imagesc(A_2)
        parse_.Results.Plotsetting([0.6 1],1)
        GOM2()
        plot_para('Filename','Anisotropy2', 'Maximize',true)
end
function ann_GOM1()
    annotation('rectangle',[0.125 0.5 0.04 0.425],'Color','k','Linewidth',2)
        annotation('textbox',[0.17 0.71 0.1 0.1],'String','$A_{e 1}$','Linestyle','none','Fontsize',40,'Color','w','Interpreter', 'latex')
    annotation('rectangle',[0.7 0.525 0.06 0.4],'Color','k','Linewidth',2)   
        annotation('textbox',[0.7 0.4 0.1 0.1],'String','$A_{e 2}$','Linestyle','none','Fontsize',40,'Color','w','Interpreter', 'latex')
    annotation('rectangle',[0.78 0.5 0.06 0.425],'Color','k','Linewidth',2)
        annotation('textbox',[0.78 0.4 0.1 0.1],'String','$A_{e 3}$','Linestyle','none','Fontsize',40,'Color','w','Interpreter', 'latex')
end
function GOM2()
    global im_size
    y =  [3100, 1600, 500];
    linew = 2;
    hold on
    plot([1, im_size(2)], y(1)*ones(2,1),'k','Linewidth', linew)
    plot([1, im_size(2)], y(2)*ones(2,1), 'k','Linewidth', linew)
    plot([1, im_size(2)], y(3)*ones(2,1), 'k','Linewidth', linew)
    text(100, y(1)-150,'$L_1$','FontSize',36,'FontName','CMU Serif Roman', 'Interpreter', 'latex')
    text(100, y(2)-150,'$L_2$','FontSize',36,'FontName','CMU Serif Roman', 'Interpreter', 'latex')
    text(100, y(3)-150,'$L_3$','FontSize',36,'FontName','CMU Serif Roman', 'Interpreter', 'latex')
    hold off
end    