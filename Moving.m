%% Parameters - C Band airborne SAR parameters
control  = 1;
% Target
x =0;
y = 19300;
dy = 400; % Length of the target area 
vx = 10; ax = 0; %azimuth
vy = 10; ay = 0; % range 
% Platform
h = 5000;
La = 2;     
vp = 150; 
f0 = 5.3e9; c = 3e8 ; lambda = c/f0 ;
dur = 4 ; %
PRF = 400; %%2.35e3
Kr = 20e12 ;
Tp = 2.5e-6; % Pulse width
B =  Kr*Tp;
Ta = 0.886 * (lambda * y)/(La *vp *1);
ts = linspace(-dur/2,dur/2,PRF*dur) ;  % Slow time space
%% Signal Space
x1 = x+ vx*ts + 0.5*ax*(ts).^2 ; y1 = (y+200) + vy*ts + 0.5*ay* (ts).^2;
R1 = sqrt(h^2 + (x1 - vp*ts).^2 + ( y1.^2 ));
x2 = x ; y2 =y;
R2 = sqrt(h^2 + (x2 - vp*ts).^2 + (y2 ).^2 );
upran = sqrt( (y+dy)^2 + h^2) ; downran =  sqrt( (y-dy)^2 + h^2);
t =  2*(downran)/c: 1/2/B : 2*(upran)/c + Tp  ;  % fast time space
s = zeros(length(ts) , length(t)) ;
% Moving 
%{%}
   wa = sinc(atan( (x1-vp*ts) ./y1 ) *La /lambda/0.886).^2;
for k = 1 : length(ts)
    td = 2* R1(k)/c ;
    s(k,:) = s(k,:) + (wa(k) * ( exp(-j*4*pi*f0*R1(k)/c + j*pi*Kr*(t- td ).^2).*(t>=td & t-td<=Tp)));
end

 
% Static
%{
    wa = sinc(atan( (x2-vp*ts) ./ y2 )  /lambda/0.886).^2;
for k = 1 : length(ts)
    td = 2* R2(k)/c ;
    s(k,:) = s(k,:) + ( wa(k) * (exp(-j*4*pi*f0*R2(k)/c + j*pi*Kr*(t- td).^2 ).*(t>=td & t-td<=Tp)));
end
%}


figure
imagesc(abs(s))
set(gca,'Ydir','normal')
%title('Raw data')
xlabel('Range (£G£n)')
ylabel('Azimuth (£G£b)')
%{%}
%% Range Compression 
% Reference sugnal 
R0 = sqrt( (y)^2+h^2);
td0 = 2*R0/c;
sr0 = exp(- j * pi*Kr*( - t -td0 ).^2).*(-t<=-td0 & -t+td0>= -Tp);

% Fourier transform
fsr0 = 1/2/B*fftshift(fft(fftshift(sr0))) ; 
for k = 1: length(ts)
    fs(k,:) = 1/2/B*fftshift(fft(fftshift(s(k,:)))) ;
    fsr(k,:) =  fsr0 .*fs(k,:); 
    sr(k,:) = fftshift(ifft(fftshift(fsr(k,:))));
end
for k = 1: length(t)
    fssr(:,k) = 1/2/B*fftshift(fft(fftshift(sr(:,k)))) ;
end


figure
imagesc(abs(sr))
set(gca,'Ydir','normal')
xlabel('Range (£G£n)')
ylabel('Azimuth (£G£b)')
%title('After range compression')
%{%}

%% Range cell migration - Keystone transform for range walk , Filter for range curveture
%Filter for range curveture
ftau = linspace(-B,B,length(t)); %1/2/B : 2*B/length(t) : 2*B 
Hc = exp( -j*2*pi*vp^2/c/R0* (ts.^2)'*ftau );
%Hc = 1;

%Filter for range walk
Hw = exp(j* 4*pi* vy*cos((y+200)/R0)/c * ts'*ftau);
%Hw =1;
f2  = fsr .* (Hc)  ;
f3 = f2 .* (Hw);

for k = 1 : length(ts)
srcmc(k,:) =  fftshift(ifft(fftshift(f3(k,:))));
temp(k,:) =   fftshift(ifft(fftshift(f2(k,:))));
end
figure
imagesc(abs((srcmc)))
set(gca,'Ydir','normal')
xlabel('Range (£G£n)')
ylabel('Azimuth (£G£b)')
%{
figure
imagesc(abs((temp)))
set(gca,'Ydir','normal')
xlabel('Range (£G£n)')
ylabel('Azimuth (£G£b)')
%}


%%  Azimuth compression 
%  cite from new approach for SAR imaging of ground moving
%  imaging and parameters estimation.....
%Ka = ( 2 *(vx - vp)^2 - ay*y ) *f0/c/R0 ;
 Ka =( 2 *(vx - vp)^2 - ay*y ) *f0/c/R0 ;
sa0 = exp(j*pi*(-ts).^2*Ka);
fsa0 = 1/PRF*fftshift(fft(fftshift(sa0))) ;
%feta = linspace(-PRF,PRF,length(ts)) ; 
%fsa0 = exp(j* pi * feta.^2 / Ka);
fdc = 2*vy*cos((y+200)/R0)/lambda;148;
Hdc = fftshift(fft(fftshift(exp(j* fdc * ts))));
%.*Hdc 
for k = 1 : length(t)
fsrcmc(: , k) = fftshift(fft(fftshift( srcmc(:,k)))) ;
fsac (: , k ) =  conj(fsa0) .*fsrcmc(:,k)';
sac (: , k) =  fftshift(ifft(fftshift( fsac(:,k) )));
end

%% Image and motion parameter
figure
imagesc(abs((sac )))
set(gca,'Ydir','normal')
xlabel('Range (£G£n)')
ylabel('Azimuth (£G£b)')

vr = c*tan( 26/1596)*PRF/2/2/B

%{
figure
subplot(2,2,1);imagesc(abs((s))) 
set(gca,'Ydir','normal')
title('(a)')
xlabel('Range (£G£n)')
ylabel('Azimuth (£G£b)')
subplot(2,2,2);imagesc(abs((sr))) 
title('(b)')
set(gca,'Ydir','normal')
xlabel('Range (£G£n)')
ylabel('Azimuth (£G£b)')
subplot(2,2,3);imagesc(abs((srcmc)))
title('(c)')
set(gca,'Ydir','normal')
xlabel('Range (£G£n)')
ylabel('Azimuth (£G£b)')
subplot(2,2,4);imagesc(abs((sac)))
title('(d)')
set(gca,'Ydir','normal')
xlabel('Range (£G£n)')
ylabel('Azimuth (£G£b)')
%}
%plot(real( exp(j*pi*ts.^2*Ka) ))