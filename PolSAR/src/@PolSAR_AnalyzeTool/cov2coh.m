function cov2coh(obj)
	% COV2COH transform coherency matrix \bar{\bar{T}} to
	% covariance matrix \bar{\bar{C}}.
    %
    % Syntax:
	%	COV2COH()
	%
    % Description:
    %   COV2COH()
	%	* If IS_BIGFILE is 1, the elements of covariance matrix will
	%	  set to empty matrix.
    %	* If IS_BIGFILE is 0, the elements of covariance matrix is kept.
	%	  
    % Outputs:
    %
    % Other m-files required: none
    % Subfunctions: none
    % MAT-files required: none
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
	%------------- <<<<< >>>>>--------------
	obj.T_11 = (obj.hh_hh + obj.vv_vv + obj.hh_vv + conj(obj.hh_vv))/2;
	obj.T_22 = (obj.hh_hh + obj.vv_vv - obj.hh_vv - conj(obj.hh_vv))/2;
	obj.T_33 = 2*obj.hv_hv;
	obj.T_12 = (obj.hh_hh - obj.vv_vv - obj.hh_vv + obj.conj(hh_vv))/2;
	obj.T_13 = obj.hh_hv + conj(obj.hv_vv);
	obj.T_23 = obj.hh_hv - conj(obj.hv_vv);
	if obj.IS_BIGFILE
		obj.hh_hh = []; obj.hv_hv = []; obj.vv_vv = []; obj.hh_hv = [];
		obj.hh_vv = []; obj.hv_vv = [];
	end
end