function setCoh2Zero(obj)
    disp('Coherency matrix is CLEARED')
    disp('Remember to RELOAD the coherency matrix!')
    obj.T_11 = []; obj.T_22 = []; 
    obj.T_33 = []; obj.T_12 = [];
	obj.T_13 = []; obj.T_23 = [];
end