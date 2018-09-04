function pauliDecomposition(obj, varargin)
    % PAULIDECOMPOSITION implement Pauli-decomposition.
    %
    % Syntax:
	%	PAULIDECOMPOSITION(Name, Value)
	%
    % Description:
    %   PAULIDECOMPOSITION(Name, Value)
    %	* Specify the filename and image scaleing by name-value pair.
    %
    % Outputs:
    %
    % Other m-files required: plot_para
    % Subfunctions: none
    % MAT-files required: none
    %
    %------------- <<<<< >>>>>--------------
    % Author: K.S. Yang
    % email: fbookzone@gmail.com
    %------------- <<<<< >>>>>--------------

    parse_ = inputParser;
    validationFcn_1_ = @(x) validateattributes(x,{'char'},{});
    addParameter(parse_,'Filename','',validationFcn_1_);
    addParameter(parse_,'Method','2sigma',validationFcn_1_);
    parse(parse_,varargin{:})

    if numel(obj.T_11) == 0
        B = (obj.hh_hh + obj.vv_vv + obj.hh_vv + conj(obj.hh_vv))/2;
	    R = (obj.hh_hh + obj.vv_vv - obj.hh_vv - conj(obj.hh_vv))/2;
        G = 2*obj.hv_hv;
    else
        B = obj.T_11;
	    R = obj.T_22;
        G = obj.T_33;
    end

    %% Pauli-decomposition
    figure
    if strcmp(parse_.Results.Method, '2sigma')
        % Show the Pauli vector with scaleing each component to [0, 255].
        image(cat(3, reshape(Stat(10*log10(R)), size(R)), ......
            reshape(Stat(10*log10(G)), size(G)), reshape(Stat(10*log10(B)), size(B))))
    else
        % Show the Pauli vector based on its power without scaleing each
        % component to [0, 255].
        up_ = 10; low_ = -35;
        Pauli = zeros([size(R), 3]);	
        % |S_vv - S_hh|^2 -> double bounce scattering 
        t_p = 10*log10(R);	
        t_p(t_p < low_) = low_;
        t_p(t_p > up_ ) = up_;
        Pauli(:,:,1) = (t_p-low_)/(up_-low_);	
        % |S_hv|^2 -> volume scattering
        t_p= 10*log10(G);
        t_p(t_p < low_) = low_;
        t_p(t_p > up_ ) = up_;
        Pauli(:,:,2) = (t_p-low_)/(up_-low_);
        % |S_vv + S_hh|^2 -> single scattering
        t_p = 10*log10(B);
        t_p(t_p < low_) = low_;
        t_p(t_p > up_ ) = up_;
        Pauli(:,:,3) = (t_p-low_)/(up_-low_);
        image(Pauli)
    end 
    set(gca,'Ydir','normal')   
    if numel(parse_.Results.Filename)
        plot_para('Maximize',true,'Filename', [obj.OUTPUT_PATH '/' parse_.Results.Filename]);
    end    
end
