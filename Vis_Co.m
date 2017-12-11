function Vis_Co()
	close all
	clear 
	s_p = 300;
	N = 20;
	thetavec = linspace(0,pi/2,N);
	phivec = linspace(0,2*pi,2*N);
	[th, ph] = meshgrid(thetavec,phivec);
	R = ones(size(th)); % should be your R(theta,phi) surface in general

	x = R.*sin(th).*cos(ph);
	y = R.*sin(th).*sin(ph);
	z = R.*cos(th);


	%{
	mu = [.25 .5];
	sig = [0.01, 0; 0, 0.02];
	F = mvnpdf([x(:) y(:)],mu, sig);
	F = reshape(F,s_p/2+1, s_p+1);
	F = ind2rgb(uint8((F/max(max(F)))*255),jet);

	q = normpdf(linspace(-1,1,s_p), x_mu, x_sig).'*normpdf(linspace(-1,1,s_p), y_mu, y_sig);
	q = fliplr(uint8(q/max(max(q))*255)).';
	C = ind2rgb(q,jet);
	figure
	imagesc(q)
	colormap jet
	colorbar
	%}

	figure
	%surf(x,y,z,F,'EdgeColor','none')
	surf(x,y,z,'FaceColor','w')
	%warp(x,y,z,C)
	view(45,30);
	xlabel('$|s_{2d}| \cos (2 \psi)$','Interpreter', 'latex')
	ylabel('$|s_{3d}| \sin (2 \psi)$','Interpreter', 'latex')
	zlabel('$|s_{1d}|$','Interpreter', 'latex')
	plot_para('Ratio',[2 2 1],'Maximize', true,'Filename', 'VisSph')

end