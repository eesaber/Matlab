function setCov2Zero(obj)
    disp('Covriance matrix is CLEARED')
    disp('Remember to RELOAD the covriance matrix!')
    obj.hh_hh = []; obj.hv_hv = [];
    obj.vv_vv = []; obj.hh_hv = [];
    obj.hh_vv = []; obj.hv_vv = [];
end