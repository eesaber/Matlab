function [map] = gen_map(n_map, col_size, row_size)
	mask = zeros(row_size, col_size, 1);
	map = zeros(row_size, col_size, n_map);
	

	% Plot the sptial distribution 
	for rr = 1 : k
    figure(k)
        imagesc(map(:,:,rr))
        %xlabel('range')
        %ylabel('azimuth')
        title(['# of ' num2str(k)])
	end
	pause
	close all

end