function [xi] = H_AF(s, order)
% This function implement high-order ambiguity function (forth order
% function).
% Usage: H_AF(s, order),where "s" is a 1xN vector and "order == 2 v 4", 2 is the regular 
% ambguity function(A.F.). 4 is the forth order A.F. The function will return a
% 2-D array 

%% Order == 4
    if order == 4
       len_s = length(s) ;
       xi_0 = floor(0.089 * len_s);
       temp = zeros(1, 2 * len_s); % s_zp -> zeropadding
       s_1 = zeros(len_s, len_s * 2 + xi_0 * 2);
       s_2 = s_1 ; s_3 = s_1; s_4 = s_1;
       for n = -len_s / 2 : len_s / 2
            s_1(n + len_s / 2 + 1, :) = [temp(1: -n + len_s / 2 ), s, temp(1: n + len_s / 2 + 2 * xi_0 )];
            %s_1 = [temp(-xi_0 - xi + 1 + xi_0   ) , s , temp(1 + xi_0 - (-xi_0 - xi + len_s/2)  ) ]
            s_2(n + len_s / 2 + 1, :) = [temp(1: n + len_s / 2 + 2* xi_0), s, temp(1: -n + len_s / 2) ];
            s_3(n + len_s / 2 + 1, :) = [temp(1: n + len_s / 2), s, temp(1: -n + len_s /2 + 2 * xi_0) ];
            s_4(n + len_s / 2 + 1, :) = [temp(1: -n + len_s /2 + 2 * xi_0), s, temp(1: n + len_s / 2) ];
       end
       xi = s_1 .* s_2 .* conj(s_3) .* conj(s_4); 
       xi = fftshift(fft(xi, 2^nextpow2(len_s)), 1) ;
       %xi = fft(xi, 2^nextpow2(len_s));
       return 
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