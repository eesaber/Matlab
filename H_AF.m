function [af] = H_AF(s, order)
% This function implement high-order ambiguity function (forth order
% function).
% Usage: H_AF(s, order),where "s" is a 1xN vector and "order == 2 v 4", 2 is the regular 
% ambguity function(A.F.). 4 is the forth order A.F. The function will return a
% 2-D array 

%{
clear
    % Testing parameters
    K_a = -100;
    t = linspace(-1 , 1, 1000 );
    s = exp(j* 2 * pi * 00 * t + j*pi * K_a *t.^2);
    %
    close all
    
    ambgfun(s ,1000,100)
%}
%% Order == 4
    if order == 4
       %xi_0 = floor(0.089 * length(s));
       xi_0 = 0 ;
       s_len = length(s) ;
       temp = zeros(1,s_len); % s_zp -> zeropadding
       size(s_zp)
       for n = -s_len : s_len
           s_1 = [temp(1:n), s(n, :), temp(n :s_len)];
           
           xi(n,:) = s_1 .* s_2 .* s_3 .* s_4; 
       end
       xi = fftshift(fft(xi), 1) ;
       mesh(xi)
    end
%%  Order == 2 
    if order == 2
        N = length(s);
        M = length(s) *0.089;
        K = nextpow2(length(s)) ; 
        % tril(ones(N,N),1) higer triangle matrix
        % tril(ones(N,N),-1) lower triangle matrix
        s_ = repmat(s, 2*N - 1, 1);
        s_star = ones(2 * N - 1, N);
        for i = 1 : N
            s_star(i,:) = [zeros(1,N-i) conj(s(N+1-i: N)) ] ;
        end 
        for i = 1 : N -1
            s_star(i + N,:) = [conj(s(1: N-i)) zeros(1,i)] ;
        end 
        %imagesc(real(s_star))
        K = s_ .* s_star;
        af = fft(K , 2^K, 2) ; % Top: \tau = T, Down: \tau = -T, Dim = 2: take FFT of row vector.
        if 1
            af = fftshift(af,2);
        end
        figure 
        mesh(abs(af))
        xlabel('$f_t$','Interpreter','latex')
        ylabel('$\xi$','Interpreter','latex')
        set(gca,'FontSize',24,'Fontname','Noto Sans CJK JP')
        colormap('Jet')
        colorbar
        %pbaspect([4 3 1])
    end
    
end