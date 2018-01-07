function [map] = gen_map(varargin)
% GEN_MAP is used to generate a spatial distribution of scatterers for simulation.
% The return of GEN_MAP is an 3D matrix of size [row, column, channel].
% How to call? Variable are in order of (row size, coloumn size, number of channel).
    switch nargin 
    	case 0
	        n_map = 8; 
	        col_size = 100; 
	        row_size = 100;
	    case 3
	    	optargs = {eps 17 @magic};
        	optargs(1:nargin) = varargin;
        	[row_size, col_size, n_map] = optargs{:};
        otherwise
        	disp('NUMBER OF INPUT VARIABLE ERROR')
    end

	map = zeros(row_size, col_size, n_map);
	rng(1); % control random seed
	
	% 1. First kind of Building.  
    map(1:row_size/2, 1:col_size/2, 2) = reshape(rand(row_size/2), [row_size/2, col_size/2, 1])>0.3;
	% 2. Second kind of Building.
    map(row_size/2+1:end, 1:col_size/2, 3) = reshape(rand(row_size/2), [row_size/2, col_size/2, 1])>0.3;
	% 3. Road 
    map(:, 1:5:col_size/2, 4) = 1;
    map(:, :, 4) = map(:, :, 4).*reshape(rand(row_size), [row_size, col_size, 1])>0.3;
	% 4. First kind of tree
	cen = [row_size/2 3*col_size/4];
	[x_, y_] = meshgrid(1:row_size, 1:col_size);
    map(:, :, 5) = (x_-cen(2)).^2+(y_-cen(1)).^2 < min(row_size/4, col_size/4).^2;
    map(:, :, 5) = map(:, :, 5).*reshape(rand(row_size), [row_size, col_size, 1])>0.3;
	% 5. Second kind of tree
	map(:, col_size/2+1:end, 6) = 1;
	map(:, :, 6) = map(:, :, 6).*reshape(rand(row_size), [row_size, col_size, 1])>0.3;
	% 6. Whatever
	map(:, :, 7) = (row_size/col_size*x_ + y_ <= row_size).*rand(row_size)>0.3;
	% 7. Whatever2
	map(:, :, 1) = (row_size/col_size*x_ + y_ >= row_size).*rand(row_size)>0.3;
	% 8. Random asymmetric	
	map(:,:, 8) = reshape(sum(map,3)<1, [row_size, col_size, 1]);
	% Plot the sptial distribution 
	figure(957)
	for rr = 1 : n_map
		subplot(2, n_map/2, rr)
		imagesc(~map(:, :, rr))
        %title(['# of ' num2str(rr)])
		set(gca, 'YDir','normal')
        colormap gray
        plot_para('Fontsize', 24,'Ratio', [1 1 1])
	end
	axis on
    plot_para('Fontsize', 24,'Ratio', [1 1 1], 'Maximize',true, 'Filename','RR')
end
