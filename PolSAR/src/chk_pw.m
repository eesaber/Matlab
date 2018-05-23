function chk_pw(sub_p)
    switch nargin
        case 0
            sub_p = '';
    end
            
	path = {'/home/akb/Code/Matlab/PolSAR/';'D:\Code\Simu\PolSAR'};
    for a = 1 : 2
        if exist(char(path(a))) == 7
            cd([char(path(a)), sub_p])
            return 
        end
    end
    pwd = input('Input the absolute path of "PolSAR"');
    cd(pwd)
end