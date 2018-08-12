function coh2cov(obj)
	% COH2COV transform coherency matrix \bar{\bar{T}} to
	% covariance matrix \bar{\bar{C}}.
    %
    % Syntax:
	%	COH2COV()
	%
    % Description:
    %   COH2COV()
	%	* If IS_BIGFILE is 1, the elements of coherency matrix will
	%	  set to empty matrix.
    %	* If IS_BIGFILE is 0, the elements of coherency matrix is kept.
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
	obj.hh_hh = (obj.T_11 + obj.T_12 + obj.T_22 + conj(obj.T_12))/2;
	obj.hv_hv = obj.T_33/2;
	obj.vv_vv = (obj.T_11 - obj.T_12 + obj.T_22 - conj(obj.T_12))/2;
	obj.hh_hv = (obj.T_13 + obj.T_23)/2;
	obj.hh_vv = (obj.T_11 - obj.T_12 + conj(obj.T_12) - obj.T_22)/2;
	obj.hv_vv = (conj(obj.T_13) - conj(obj.T_23))/2;
	if obj.IS_BIGFILE
		obj.T_11 = []; obj.T_22 = []; obj.T_33 = []; obj.T_12 = [];
		obj.T_13 = []; obj.T_23 = [];
	end
end