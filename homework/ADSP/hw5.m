t = 0:.01:2*pi;
x = t;
y = sin(2*t) ;
comet(x,y);


%{
[a,b] = NTTm(4,5);

function [A, B] = NTTm(N, M)
	prime = primes(M-1);
	alpha = 0;
	for i = 1:length(prime)
		if (mod((M-1), prime(i))==0)
			all = (M-1)./prime(i);
			if (prime(i)^all~=1)
				alpha = prime(i);
			end
		else
			continue;
		end
	end
	col = (0:N-1)';
	rows = col';
	nk = col*rows;
	A = mod(alpha.^nk, M);
	B = mod((M-1)* alpha.^(-nk), M);
end
%}