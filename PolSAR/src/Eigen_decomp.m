function [H, alpha_bar] = Eigen_decomp(varargin)
    %% H_alpha plot
    parse_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'logical'},{});
    validationFcn_2_ = @(x) validateattributes(x,{'char'},{'nonempty'});
    validationFcn_3_ = @(x) validateattributes(x,{'numeric'},{'3d'});
    addParameter(parse_,'Calculate',0,validationFcn_1_);
    addParameter(parse_,'FileName','',validationFcn_2_);
    addParameter(parse_,'T',[],validationFcn_3_);
    parse(parse_,varargin{:})
    
    global dir im_size
    if parse_.Results.Calculate
        ind_ = randperm(im_size(1)*im_size(2));
        figure
            H_Alpha(parse_.Results.T(:,:,ind_(1:2000)))
        clear ind_
        % Plot each entropy of each pixel
        t = cputime;
        [~,~,num]= size(parse_.Results.T);
        H = zeros(1, num); A_1 = zeros(1, num); A_2 = zeros(1, num);
        alpha_bar = zeros(1, num);
        for r = 1 : num
            if sum(sum(parse_.Results.FileNameT(:,:,r)== zeros(3,3))) == 9
                continue;
            end
            [U, L] = eig(T(:,:,r));
            P = diag(L)/sum(diag(L));
            H(r) = -sum(P.*log(P)/log(3));
            A_temp = sort(diag(L),'descend');
            A_1(r) = (A_temp(2)-A_temp(3))/(A_temp(2)+A_temp(3));
            A_2(r) = (A_temp(1)-A_temp(2))/(A_temp(1)+A_temp(2));
            alpha_bar(r) = sum(P'.*acosd(abs(U(1,:))));
        end
        e = cputime -t;
        disp(e)
        %%
        H = real(reshape(H, im_size));
        A_1 = reshape(A_1, im_size);
        A_2 = reshape(A_2, im_size);
        alpha_bar = reshape(alpha_bar, im_size);
        save([dir parse_.Results.FileName '.mat'],'-v7.3', 'H', 'A_1', 'A_2', 'alpha_bar');
    else 
        load([dir parse_.Results.FileName '.mat'])
    end

    %%
    figure
        imagesc(H)
        %Plotsetting_GOM1([0 1])
        %ann_GOM1()
        ann_GOM2()
        Plotsetting_GOM2([0 1],1)
        plot_para('Filename','Entropy', 'Maximize',true)
    figure
        imagesc(alpha_bar)
        %Plotsetting_GOM1([0 90])
        Plotsetting_GOM2([0 90],1)
        plot_para('Filename','alpha', 'Maximize',true)
    figure
        imagesc(A_1)
        %Plotsetting_GOM1([0 0.4])
        Plotsetting_GOM2([0 0.8],1)
        plot_para('Filename','Anisotropy1', 'Maximize',true)
    figure
        imagesc(A_2)
        %Plotsetting_GOM1([0 1])
        Plotsetting_GOM2([0.6 1],1)
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
function ann_GOM2()
     annotation('rectangle',[0.3 0.65 0.04 0.25],'Color','k','Linewidth',2)
        annotation('textbox',[0.25 0.75 0.1 0.1],'String','$A_{e 1}$','Linestyle','none','Fontsize',40,'Interpreter', 'latex','Color','w')
    annotation('rectangle',[0.4 0.55 0.18 0.375],'Color','k','Linewidth',2)
        annotation('textbox',[0.6 0.75 0.1 0.1],'String','$A_{e 2}$','Linestyle','none','Fontsize',40,'Interpreter', 'latex','Color','w')
end