function [AF] = GAF(s, total_ord, n, k)
% This function implement the "General Ambiguity function (GAF)" 
% Usage: GAF(s, total_ord, n), where "s" is the 1xN array, "total_ord" is
% the order, at begin, of the polynomial phase signal. "n" is the order
% "NOW". ?�k??is used in FFT which do 2^(nextpow2(length(s))+k) zero padding.
% Return an n-th order ambiguity function.
	t_mov = floor(length(s) / total_ord) ;
	for qq = 1 : n - 1
		s =  s(t_mov + 1: end) .* conj(s(1: length(s) - t_mov) ) ;
	end
	%figure
	%plot(real(s))
	AF = fftshift(fft(s, 2^(nextpow2(t_mov * total_ord)+k)) );
end