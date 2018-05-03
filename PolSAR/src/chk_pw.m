function chk_pw()
	path = {'/home/akb/Code/Matlab/PolSAR';'D:\Code\Simu\PolSAR'};
    for a = 1 : 2
        if exist(char(path(a))) == 7
            cd(char(path(a)))
            return 
        end
    end
    pwd = input('Input the absolute path of "PolSAR"');
    cd(pwd)
end