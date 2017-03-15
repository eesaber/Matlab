function [AF] = GAF(s, total_ord, n)
% This function implement the "General Ambiguity function (GAF)" 
% Usage: GAF(s, total_ord, n), where "s" is the 1xN array, "total_ord" is
% the order, at begin, of the polynomial phase signal. "n" is the order
% "NOW". Return an n-th order ambiguity function.
	if mod(length(s), total_ord * 2)
		s = [s zeros(1, total_ord * 2 - mod(length(s), total_ord * 2))];
	end
	t_mov = length(s) / total_ord / 2;
	for qq = 1 : n - 1
		s = [zeros(1, t_mov) s(1: length(s) - t_mov) ] .* conj( [ s(t_mov + 1:end) zeros(1, t_mov)] ) ;
	end
	AF = fftshift(fft(s, 2^nextpow2(length(s))));
end