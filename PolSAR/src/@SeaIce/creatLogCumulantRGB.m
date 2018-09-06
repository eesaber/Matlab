function creatLogCumulantRGB(obj)
    if 1
        obj.im = cat(3, reshape(intensityMapping(obj.kai_1,'Bit','uint16'), size(obj.kai_1)),...
            reshape(intensityMapping(obj.kai_2,'Bit','uint16'), size(obj.kai_2)), ...
            reshape(intensityMapping(obj.kai_3,'Bit','uint16'), size(obj.kai_3)));
        obj.im = obj.im(2:end-1,3:end-2,:);
    else
        R = obj.kai_1(2:end-1,3:end-2);
        G = obj.kai_2(2:end-1,3:end-2);
        B = obj.kai_3(2:end-1,3:end-2);
        obj.im = cat(3, reshape(intensityMapping(R,'Bit','uint16'), size(R)),...
            reshape(intensityMapping(G,'Bit','uint16'), size(G)), ...
            reshape(intensityMapping(B,'Bit','uint16'), size(B)));
        clear R G B
    end
end