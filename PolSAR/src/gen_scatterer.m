function gen_scatterer()
%function [S, Sigma] = gen_scatterer()
% GEN_SCATTERER is used to generate dictionary matrix for simulation.
% The return of GEN_SCATTERER is an 3D matrix of size [row, column, # of atom].
	rng(180);
	n_atom = 8;
	range = 10;
	phi = linspace(0,pi,180);
	A = zeros(3,8);
	I = [1, 0; 0, 1];
    Sigma = zeros(8,9);
	% surface 3 atoms
	for k = 1 : 3
		S = abs(range*rand(2).*I);
		theta = rand(1)*pi;
		R_k = [1, 0, 0; 0, cos(2*theta), -sin(2*theta); 0, sin(2*theta), cos(2*theta)];
        Sigma(k,:) = 1/2*reshape([S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)]*[S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)]', [1 9]);
		A(:,k) = 1/sqrt(2)*R_k*[S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)];
	end
	% dihedral 2 big sigma 2 small sigma
	for k = 4 : 7
		S = range*(1-randn(2) + 1j-1j*randn(2)).*I;
		theta = rand(1)*pi;
		R_k = [1, 0, 0; 0, cos(2*theta), -sin(2*theta); 0, sin(2*theta), cos(2*theta)];
        Sigma(k,:) = 1/2*reshape([S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)]*[S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)]', [1 9]);
		A(:,k) = 1/sqrt(2)*R_k*[S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)];
	end
	% volume
    Base = [1, 0; 0, 0];
    S = zeros(2,2);
    for  k = 1 : numel(phi)
        R_i = [cos(phi(k)), sin(phi(k)); -sin(phi(k)), cos(phi(k))];
        S = S + sin(phi(k))*R_i*Base*R_i.'*(phi(2)-phi(1));
    end
    A(:,8) = 1/sqrt(2)*[S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)];
	Sigma(8,:) = 1/2*reshape([S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)]*[S(1,1)+S(2,2); S(1,1)-S(2,2); 2*S(1,2)]', [1 9]);
    temp = (abs(Sigma) < 0.0001).*repmat([1 0 0 0 1 0 0 0 1],[8,1]);
    Sigma = Sigma + temp;
 
   
	%%
    for qq = 1 : n_atom
		Vis_Co(A(:,qq),'Sigma',0.01*reshape(abs(Sigma(qq,:)),[3 3]),'Subplot',true,'Filename',num2str(qq))
	end
	%%
    label = ['(a)';'(b)';'(c)';'(d)';'(e)';'(f)';'(g)';'(h)'];
    for qq = 1 : n_atom
		subplot(2, n_atom/2, qq)
		I = imread([num2str(qq),'.jpg']);
		delete([num2str(qq),'.jpg'])
		imshow(I)
		q = size(I);
		%set(gca,'XTick',linspace(1,q(1),3), 'YTick',linspace(1,q(1),3),'XTickLabel',{'-1','0','1'},'YTickLabel',{'-1','0','1'})
		set(gca,'TickDir','in','XTick',linspace(1,q(1),3), 'YTick',linspace(1,q(1),5),'XTickLabel',{'-1','0','1'},'YTickLabel',{'-1','-0.5','0','0.5','1'})
		xlabel(label(qq,:),'Interpreter', 'latex')
		plot_para('Fontsize', 24)
	end
	plot_para('Fontsize', 24,'Ratio', [1 1 1], 'Maximize',true, 'Filename','SimAtom')
	
end