%% Q.1

[data,fs]= audioread('Audio3.wav');

yf = fft(data)/fs;



F_data = fftshift(fft(data,length(data))/fs);
S=F_data.*conj(F_data);
df = linspace(-fs/2,fs/2,length(F_data));


dt=linspace(0,length(data)/fs,length(data));
plot(dt,data);
title('Time Domain of signal');

figure;
plot(df, abs(F_data));
title('Fourier Transform of signal')

S=S';
p=length(S);
p=p-1;

for i= p/2 +1:p
   
    a=sum(S(p/2:i));
    b=sum(S(p/2:end));
    
    if a>0.99*b
       k=i;
       break;
        
    end
    
    
end

BW=(k-p/2)*(fs/p);

%% 2

load('LPF1_HW3.mat');


x = filter(LPF1,data);

xf = fftshift(fft(x,length(data))/fs);
df = linspace(-fsn/2,fsn/2,length(xf));
plot(df, abs(xf));

title('Filtered signal')

fm = 6500;



ups = 10;
fsn = fs * ups;
ux = interp(x, ups);
uxf = fftshift(fft(ux,length(ux))/fsn);

df = linspace(-fsn/2,fsn/2,length(uxf));
figure
plot(df, abs(ux));

title('Upsampled signal')

%% mod

xt =ux;
b = 3;
delf = b * fm;
fdel = delf / max(xt);
fc = 50000;

q=length(xt);


dt=linspace(0,q/(fsn),q);
intx=cumtrapz(dt,xt);
intx=intx';

dt = 1/fsn:1/fsn:q/fsn;
xc = cos(2*pi*fc*dt + 2*pi*fdel*intx);


xcf = fftshift(fft(xc,length(xc))/fsn);

df = linspace(-fsn/2,fsn/2,length(xcf));
plot(df, abs(xcf));

title('Modulated signal')

%% Dem



xcp = zeros(1,length(xc));

for i = 2: length(xc)-1
    xcp(i) =(xc(i+1)-xc(i-1))*fsn/ 2;
end

xcpf = fftshift(fft(xcp,length(xcp))/fsn);

df = linspace(-fsn/2,fsn/2,length(xcpf));

plot(df, abs(xcf));

title('Derivated Signal')
%%
ex=abs(xcp);

ex=ex-mean(ex);

exf = fftshift(fft(ex,length(ex))/fsn);

df = linspace(-fsn/2,fsn/2,length(exf));

plot(df, abs(exf));

title('Demodulated signal')


%%

load('LPF2_HW3.mat')
r = filter(LPF2,ex);

rf = fftshift(fft(r,length(r))/fsn);

df = linspace(-fsn/2,fsn/2,length(r));

plot(df, abs(rf));

title('Final Filtered Signal')
%% Downsampling
ds = downsample(r, ups);

dsf = fftshift(fft(ds,length(ds))/fs);

df = linspace(-fs/2,fs/2,length(dsf));



plot(df,abs(dsf));

title('final Signal')

dsf=dsf';

mse= immse(abs(dsf)*(max(abs(fft(data)/fs))/max(abs(dsf))),abs(fft(data)/fs))


%% Q.2
clc
clear

fs = 50000;
fm = 10;
dt=0:1/fs:10;
x = cos(2*pi*fm*dt);
fdel = 100;
fc = 1000;

q=length(x);
t=linspace(0,q/(fs),q);
intx=cumtrapz(t,x);


v = cos(2*pi*fc*dt + 2*pi*fdel*intx);

vf = fftshift(fft(v,length(v))/fs);

df = linspace(-fs/2,fs/2,length(vf));

plot(df, abs(vf));

title('Fourier Transform of Modulated Signal')

figure;
plot(dt,x);
xlim([0 0.1])
title('Message Signal')

figure
plot(dt,v);
xlim([0 0.1])
title('Modulated Signal')
%%

S=vf.*conj(vf);
df = linspace(-fs/2,fs/2,length(S));

S=S';
p=length(S);
p=p-1;

for i= p/2 +1:p
   
    a=sum(S(p/2:i));
    b=sum(S(p/2:end));
    
    if a>0.99*b
       k=i;
       break;
        
    end
    
    
end

BW=(k-p/2)*(fs/p);
%% b

y = v .^ 3;

yf = fftshift(fft(y,length(v))/fs);

df = linspace(-fs/2,fs/2,length(yf));

plot(df, abs(yf))

title('Modulated Signal after nonlinear block')
%%
load('LPF3_HW3.mat');

yy =0.75*filter(LPF3,y);

de = mean(grpdelay(LPF3));
yy(1:de) = [];


yyf = fftshift(fft(yy,length(yy))/fs);
df = linspace(-fs/2,fs/2,length(yyf));


plot(df, abs(yyf));
title('Signal without distortion')
%% c
lv = sign(v);
lvv = filter(LPF3, lv);
lvv(1:de) = [];

lvf = fftshift(fft(lv,length(lv))/fs);
df = linspace(-fs/2,fs/2,length(lvf));
plot(df, abs(lvf));
title('Limited Signal')
%%

lvvf = fftshift(fft(lvv,length(lvv))/fs);
df = linspace(-fs/2,fs/2,length(lvvf));
figure
plot(df, abs(lvvf));
title('Limited Filtered Signal')
%%

lvvp = zeros(1,length(lvv));
for i = 2: length(lvv)-1
    lvvp(i) =(lvv(i+1)-lvv(i-1))*fs/ 2;
end


lvvpf = fftshift(fft(lvvp,length(lvvp))/fs);
df = linspace(-fs/2,fs/2,length(lvvpf));

plot(df, abs(lvvpf));
title('Derivated Signal')

%%
ex = abs(lvvp);
exf = fftshift(fft(ex,length(ex))/fs);
df = linspace(-fs/2,fs/2,length(exf));
plot(df, abs(exf));
title('Envelope')
%%
o=fft(ex);
o(1)=0;
dm=ifft(o);

dmf = fftshift(fft(dm,length(dm))/fs);
df = linspace(-fs/2,fs/2,length(dmf));

plot(df, abs(dmf));
title('Demodulated signal')
%%
load('LPF4_HW3.mat');

dmm = filter(LPF4,dm);
de = mean(grpdelay(LPF4));
dmm(1:de) = [];
dmmf = fftshift(fft(dmm,length(dmm))/fs);
df = linspace(-fs/2,fs/2,length(dmmf));

plot(df, abs(dmmf));
title('first filtering of output')

%%
load('LPF5_HW3.mat');
oo = filter(LPF5,dmm);
de = mean(grpdelay(LPF5));
oo(1:de) = [];
oof = fftshift(fft(oo,length(oo))/fs);
df = linspace(-fs/2,fs/2,length(oof));
plot(df, abs(oof));
title('final Output')
%%
dt=linspace(0,10,length(oo));
figure
plot(dt,oo)
xlim([0 0.2])
title('final Signal in Time')
%% a
dt=0:1/fs:10;
w = v + 0.1* cos(2*pi*(fc+100/(2*pi))*t+pi/3);

wp = zeros(1,length(w));
for i = 2: length(w)-1
    wp(i) =(w(i+1)-w(i-1))*fs/ 2;
end


wpf = fftshift(fft(wp,length(wp))/fs);
df = linspace(-fs/2,fs/2,length(wpf));



exw = abs(wp);
exwf = fftshift(fft(exw,length(exw))/fs);
df = linspace(-fs/2,fs/2,length(exwf));


ow=fft(exw);
ow(1)=0;
dmw=ifft(ow);

dmwf = fftshift(fft(dmw,length(dmw))/fs);
df = linspace(-fs/2,fs/2,length(dmwf));



dmmw = filter(LPF4,dmw);
de = mean(grpdelay(LPF4));
dmmw(1:de) = [];
dmmwf = fftshift(fft(dmmw,length(dmmw))/fs);
df = linspace(-fs/2,fs/2,length(dmmwf));


oow = filter(LPF5,dmmw);
de = mean(grpdelay(LPF5));
oow(1:de) = [];
oowf = fftshift(fft(oow,length(oow))/fs);
df = linspace(-fs/2,fs/2,length(oowf));


dt=linspace(0,10,length(oow));
figure
plot(dt,oow)
title('output Signal in Time')
xlim([0.1 0.2])
%% b

dt=0:1/fs:10;
[b, a]= butter(1, 10/(fs/2));
freqz(a,b)
vf =3*filter(a,b,v);
w = vf + 0.1* cos(2*pi*(fc+100/(2*pi))*dt+pi/3);
wf= 1/3*filter(b,a,w);
figure
freqz(b,a)
wfp = zeros(1,length(wf));
for i = 2: length(wf)-1
    wfp(i) =(wf(i+1)-wf(i-1))*fs/ 2;
end


wfpf = fftshift(fft(wfp,length(wfp))/fs);
df = linspace(-fs/2,fs/2,length(wfpf));



exww = abs(wfp);
exwwf = fftshift(fft(exww,length(exww))/fs);
df = linspace(-fs/2,fs/2,length(exwwf));


oww=fft(exww);
oww(1)=0;
dmww=ifft(oww);

dmwwf = fftshift(fft(dmww,length(dmww))/fs);
df = linspace(-fs/2,fs/2,length(dmwwf));



dmmww = filter(LPF4,dmww);
de = mean(grpdelay(LPF4));
dmmww(1:de) = [];
dmmwwf = fftshift(fft(dmmww,length(dmmww))/fs);
df = linspace(-fs/2,fs/2,length(dmmwwf));



ooww = filter(LPF5,dmmww);
de = mean(grpdelay(LPF5));
ooww(1:de) = [];
oowwf = fftshift(fft(ooww,length(ooww))/fs);
df = linspace(-fs/2,fs/2,length(oowwf));



dt=linspace(0,10,length(ooww));
figure
plot(dt,ooww)
title('output Signal in Time')

xlim([0.1 0.2])
