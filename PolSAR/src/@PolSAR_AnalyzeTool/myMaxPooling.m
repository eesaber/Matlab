function y = myMaxPooling(obj, x, r_size, c_size)
    switch nargin
        case 2 
            r_size = 10;
            c_size = 20;
        case 3
            c_size = 20;            
    end
    y = zeros(ceil(size(x,1)/r_size), ceil(size(x,2)/c_size));
    pad_im = padarray(x,[r_size*ceil(size(x,1)/r_size)-size(x,1),...
        c_size*ceil(size(x,2)/c_size)-size(x,2)],'replicate','post');
    for c_it = 1 : size(y,2)
        for r_it = 1 : size(y,1)
            y(r_it, c_it) = max(max(pad_im(1+r_size*(r_it-1):r_size*r_it,...
                1+c_size*(c_it-1):c_size*c_it)));
        end
    end
end