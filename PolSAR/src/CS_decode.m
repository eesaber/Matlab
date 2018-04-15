function [hh_hh, hv_hv, vv_vv, hh_hv, hh_vv, hv_vv, dem] =  CS_decode()
    % CS_DECODE is used to decompress the "AIRSAR compressed Stokes
    % matrix". The decompress algorithm is in "AIRSAR data format". The 
    % return value of CS_DECODE is [hh_hh, hv_hv, vv_vv, hh_hv, hh_vv,
    % hv_vv].
    
    chk_pw()
    % mission
	mission_num = 1;
    switch mission_num
		case 1
			disp('AIRSAR Mission: Camp Roberts')
			im_size = [4028 1151];
			%im_size_c = im_size.*[2 ,1];
			dir = '/media/akb/2026EF9426EF696C/raw_data/Camp_Roberts/';
			task = 'ts554_';
		otherwise 
			error('You need to select a mission')
    end
    clear temp
    disp('Parsing input file...')
    
    fid = fopen([dir task 'l.datgr'],'r','n'); 
    hd_1 = fscanf(fid, '%c',[1, 11510]); % Read first line header.
    hd_1 = reshape(hd_1(1:1000), [50,20]).';
    disp(hd_1)
    hd_2 = fscanf(fid, '%c',[1, 11510]); % Read second line header.
    hd_2 = reshape(hd_2(1:5000), [50,100]).';
    hd_3 = fscanf(fid, '%c',[1, 11510]); % Read third line header.
    hd_3 = reshape(hd_3(1:11500), [50,230]).';
    for k = 1 : 3 
        [~] = fscanf(fid, '%c',[1, 11510]); % Neglect the rest of the header.
    end
    cs = fread(fid,[10*1150*4,1007],'int8');
    gen_fac = 1;
    M_11 = (cs(2:10:end, :)/254 + 1.5).*2.^(cs(1:10:end, :))*gen_fac;
    M_11(M_11 < 1e-30) = 0; % Set the number which is closed to 0 to 0.
    M_12 = cs(3:10:end, :).*M_11/127;
    M_13 = (-1+2*(cs(4:10:end, :)>=0)).*(cs(4:10:end, :)/127).^2.*M_11;
    M_14 = (-1+2*(cs(5:10:end, :)>=0)).*(cs(5:10:end, :)/127).^2.*M_11;
    M_23 = (-1+2*(cs(6:10:end, :)>=0)).*(cs(6:10:end, :)/127).^2.*M_11;
    M_24 = (-1+2*(cs(7:10:end, :)>=0)).*(cs(7:10:end, :)/127).^2.*M_11;
    M_33 = cs(8:10:end, :).*M_11/127;
    M_34 = cs(9:10:end, :).*M_11/127;
    M_44 = cs(10:10:end, :).*M_11/127;
    M_22 = M_11-M_33-M_44;
   
    hh_hh = M_11 + M_22 + 2*M_12;
    hv_hv = (M_33 + M_44);
    vv_vv = (hh_hh - 4*M_12);
    %hh_hh = rot90(hh_hh);
    hh_hv = ((M_23 -1j*M_24) + (M_13-1j*M_14));
    hv_vv = ((M_13-1j*M_14) - (M_23 -1j*M_24));
    hh_vv = (M_33 - M_44 - 2j*M_34);
    Pauli_decomp((hh_hh+vv_vv-hh_vv-conj(hh_vv))/2, 2*hv_hv,(hh_hh+vv_vv+hh_vv+conj(hh_vv))/2,'Filename','qq')
    dem = AIRSAR_DEM([dir task 'c.demi2']);
end
function [dem] = AIRSAR_DEM(dir)
    fid = fopen(dir,'r','n'); 
    hd_1 = fscanf(fid, '%c',[1, 2302]); % Read first line header.
    hd_1 = reshape(hd_1(1:1000), [50,20]).';
    [~] = fscanf(fid, '%c',[1, 2302*3]); % Skip 2~4 line header.
    hd_3 = fscanf(fid, '%c',[1, 2302]); % Read third line header.
    hd_3 = reshape(hd_3(1:1000), [50,20]).';
    
    DN = fread(fid,[1151, 4083],'int16');
    dem = 0.04112*(DN) + 16.5;
end