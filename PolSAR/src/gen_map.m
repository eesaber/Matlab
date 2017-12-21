function [map] = gen_map(n_map, col_size, row_size)
	mask = zeros(row_size, col_size, 1);
	map = zeros(row_size, col_size, n_map);
	rng(1024);
	temp = [ones(col_size*row_size*0.3) zeros(col_size*row_size*0.7)];
	% 1. Random asymmetric.
	map(:, :, 1) = map(:, :, 1) + reshape(temp(randperm(numel(temp))), [row_size, col_size, 1]); 
	% 2. First kind of Building.

	% 3. Second kind of Building.

	% 4. Road 

	% 5. First kind of tree

	% 6. Second kind of tree

	% 7. 

	% Plot the sptial distribution 
	for rr = 1 : k
    figure(k)
        imagesc(map(:, :, rr))
        %xlabel('range')
        %ylabel('azimuth')
        title(['# of ' num2str(k)])
	end
	pause
	close all

end