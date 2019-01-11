%%
f2()
f3()

function f1()
fid = fopen('/home/akb/下載/Haywrd_23501_18039_014_180801_L090HH_CX_01.slc','r','ieee-le'); 
OFFSET = 2*9900*95000*4;
fseek(fid, OFFSET, 'bof');
hh = fread(fid, [9900 95000],'real*4=>single',4);
fclose(fid);
fid = fopen('/home/akb/下載/Haywrd_23501_18039_014_180801_L090HH_CX_01.slc','r','ieee-le'); 
fseek(fid, OFFSET+4, 'bof');
ang = fread(fid, [9900 95000],'real*4=>single',4);
hh = hh.*exp(1j*ang);
clear ang;


figure
imagesc(10*log10(abs(hh(:,40000:70000))))
set(gca,'clim',[-30 20],'Colormap',colormap('jet'))
end

function f2()
fid = fopen('/home/akb/下載/Haywrd_23501_18039_014_180801_L090HH_CX_01.slc','r','ieee-le'); 
OFFSET = 2*9900*95000*4;
fseek(fid, OFFSET, 'bof');
hh = fread(fid, [9900*2 95000],'real*4=>single');
hh = hh(1:2:end, :).*exp(1j*hh(2:2:end, :));

figure
imagesc(10*log10(abs(hh(:,40000:70000))))
set(gca,'clim',[-30 20],'Colormap',colormap('jet'))
end
function f3()
fid = fopen('/home/akb/下載/Haywrd_23501_18039_014_180801_L090HH_CX_01.slc','r','ieee-le'); 
OFFSET = 2*9900*95000*4;
fseek(fid, OFFSET, 'bof');
hh = fread(fid, [9900*2 95000],'real*4=>single');
a = hh(1:2:end, :).*exp(1j*hh(2:2:end, :));

figure
imagesc(10*log10(abs(a(:,40000:70000))))
set(gca,'clim',[-30 20],'Colormap',colormap('jet'))
end
