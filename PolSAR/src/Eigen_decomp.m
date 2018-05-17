function [H, alpha_bar] = Eigen_decomp(T, span)
    %% H_alpha plot
    global dir im_size
    if(0)
        ind_ = randperm(numel(span));
        figure
            H_Alpha(T(:,:,ind_(1:2000)))
        clear ind_
        % Plot each entropy of each pixel
        t = cputime;
        [~,~,num]= size(T);
        H = zeros(1, num); A_1 = zeros(1, num); A_2 = zeros(1, num);
        alpha_bar = zeros(1, num);
        for r = 1 : num
            if sum(sum(T(:,:,r)== zeros(3,3))) == 9
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
        save([dir 'eigen.mat'],'-v7.3', 'H', 'A_1', 'A_2', 'alpha_bar');
    else 
        load([dir 'eigen.mat'])
    end
    %{
    close all
    line = zeros(1,size(H,1));
    global dir;
    load([dir area '_hgt.mat']);
    for a=1 : size(H,1)
        line(a) = find(hgt(a,:)+10000,1);
    end
    %}
   %%
    figure
        imagesc(H)
        Plotsetting_1([0 1])
        annotation('rectangle',[0.125 0.5 0.04 0.425],'Color','k','Linewidth',2)
            annotation('textbox',[0.17 0.71 0.1 0.1],'String','$A_{E1}$','Linestyle','none','Fontsize',40,'Color','w','Interpreter', 'latex')
        annotation('rectangle',[0.7 0.525 0.06 0.4],'Color','k','Linewidth',2)   
            annotation('textbox',[0.7 0.4 0.1 0.1],'String','$A_{E2}$','Linestyle','none','Fontsize',40,'Color','w','Interpreter', 'latex')
        annotation('rectangle',[0.78 0.5 0.06 0.425],'Color','k','Linewidth',2)
            annotation('textbox',[0.78 0.4 0.1 0.1],'String','$A_{E3}$','Linestyle','none','Fontsize',40,'Color','w','Interpreter', 'latex')
        xlabel('Azimuth (km)')
        ylabel('Range (km)')  
        plot_para('Filename','output/Entropy', 'Maximize',true)
   %%
    figure
        imagesc(alpha_bar)
        Plotsetting_1([0 90])    
        xlabel('Azimuth (km)')
        ylabel('Range (km)')  
        plot_para('Filename','output/alpha', 'Maximize',true)
    figure
        imagesc(A_1)
        Plotsetting_1([0 0.6])
        xlabel('Azimuth (km)')
        ylabel('Range (km)')  
        plot_para('Filename','output/Anisotropy1', 'Maximize',true)
    figure
        imagesc(A_2)
        Plotsetting_1([0 1])
        xlabel('Azimuth (km)')
        ylabel('Range (km)')  
        plot_para('Filename','output/Anisotropy2', 'Maximize',true)
end