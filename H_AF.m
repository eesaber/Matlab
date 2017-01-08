%function [af] = H_AF(s)
%%
clear
    % Testing parameters
    K_a = -100;
    t = linspace(-1 , 1, 1000 );
    s = exp(j* 2 * pi * 00 * t + j*pi * K_a *t.^2);
    %
    close all
%%    
    N = length(s);
    M = nextpow2(length(s)) ; 
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
    af = fft(K , 2^M, 2) ; % Top: \tau = T, Down: \tau = -T, Dim = 2: take FFT of row vector.
    if 1
        af = fftshift(af,2);
    end
    %%
    figure 
    mesh(abs(af))
    xlabel('$f_t$','Interpreter','latex')
	ylabel('$\xi$','Interpreter','latex')
    set(gca,'FontSize',24,'Fontname','Noto Sans CJK JP')
	colormap('Jet')
    colorbar
    %pbaspect([4 3 1])
%end