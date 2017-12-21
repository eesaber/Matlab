clear, close all
	trans = 1;
	k = 10;
	if trans == 0
		F = linspace(-0.5,0.5,2*k+1);
	else 
		F = linspace(-0.5,0.5,2*k+1);
		F(k/2+1:3*k/2+1) = 0.5./(1+exp(-6*F(k/2+1:3*k/2+1)))-0.25;
	end
	H = 1j * 2 * pi * F;
	h_ip = fftshift(ifft(H));
	%h_ip()
	figure
	stem(linspace(-0.5,0.5,2*k+1),imag(H),'k','Linewidth',2)
	xlim([-0.5 0.5])
	xlabel('F')
	ylabel('Imag(H(\omega))')
	title('frequency responce')
	%plot_para('Maximize',true,'Filename','1')
	
	figure
	stem(-k:k,real(h_ip),'k','Linewidth',2)
	title('impulse responce')
	%plot_para('Maximize',true,'Filename','2')