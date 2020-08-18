%Q.1
syms t w

fc=10;

x=exp(-t^2);
X=fourier(x);
c=cos(2*pi*fc*t);

v=x*c;
V=fourier(v);

r(t)=v*c;

R=fourier(r);

subplot(1,2,1)
fplot(x)
title('x in Time domain')

subplot(1,2,2)
fplot(X)
xlim([-20 20])

title('x in Frequency domain')

figure();
subplot(1,2,1)
fplot(v)
title('v in Time domain')

subplot(1,2,2)
fplot(V)
xlim([-80 80])
ylim([0 1])
title('v in Frequency domain')

figure();
subplot(1,2,1)
fplot(r)
title('r in Time domain')

subplot(1,2,2)
fplot(R)
xlim([-160 160])
ylim([0 1])
title('r in Frequency domain')

h(t)=20*sinc(10*t);
%%
d=int(r(w)*h(w-t),w,-inf,inf);

D(w)=fourier(d);
subplot(1,2,1)
fplot(d)
title('d in Time domain')

subplot(1,2,2)
fplot(X)
title('d in frequency domain')
%%
%Q.2

x='message.wav';

[data,fs]=audioread(x);

yf=fft(data,length(data))/fs;

Yf=fftshift(fft(data,length(data))/fs);

S=Yf.*conj(Yf);
p=length(S);
df=-fs/2:fs/length(data):fs/2;
df(p+1)=[];
plot(df,S);
title('Spectral density');
figure();
plot(df,abs(yf));

title('Frequency domain');

for i= p/2 +1:p
   
    a=sum(S(p/2:i));
    b=sum(S(p/2:end));
    
    if a>0.99*b
       k=i;
       break;
        
    end
    
    
end

BW=(k-p/2)*(fs/p);


%%

upr = 14;

fss=fs*upr;
qq = upsample(data, upr);

Q=fftshift(fft(qq,length(qq))/fss);
df=-fss/2:fss/length(Q):fss/2;

df(length(Q)+1)=[];

plot(df,abs(Q))
title('Frequency domain');
%%


load LPF_22

qu= filter(LPF,qq);

Qu=fftshift(fft(qu,length(qu))/fss);
df=-fss/2:fss/length(Qu):fss/2;
df(length(Qu)+1)=[];

plot(df,abs(Qu))
title('Frequency domain');
%%

dt = 1/fss:1/fss:length(qu)/(fss);
fc = 102000;

ca=2*cos(2*pi*fc*dt);


mqu=ca.*qu';

mQu=fftshift(fft(mqu,length(mqu))/fss);
df=-fss/2:fss/length(mQu):fss/2;
df(length(mQu)+1)=[];
figure()
plot(df,abs(mQu))
title('Frequency domain with carrier');
%%
load HPF_22

ussb= filter(HPF,mqu);

Ussb=fftshift(fft(ussb,length(ussb))/fss);
df=-fss/2:fss/length(Ussb):fss/2;
df(length(Ussb)+1)=[];

plot(df,abs(Ussb))
title('Frequency domain');
%%
dt = 1/fss:1/fss:length(qu)/(fss);
fc = 3000;

ca=2*cos(2*pi*fc*dt);


mqu=ca.*qu';

mQu=fftshift(fft(mqu,length(mqu))/fss);
df=-fss/2:fss/length(mQu):fss/2;
df(length(mQu)+1)=[];

load HPFU1_22

ussb= filter(HPFU1,mqu);

Ussb=fftshift(fft(ussb,length(ussb))/fss);
df=-fss/2:fss/length(Ussb):fss/2;
df(length(Ussb)+1)=[];

plot(df,abs(Ussb))
title('Frequency domain with carrier 3k');
%%
dt = 1/fss:1/fss:length(ussb)/(fss);
fc = 99000;

ca=2*cos(2*pi*fc*dt);
mqu=ca.*ussb;

load HPF_22

us= filter(HPF,mqu);

Ussb=fftshift(fft(us,length(us))/fss);
df=-fss/2:fss/length(Ussb):fss/2;
df(length(Ussb)+1)=[];

plot(df,abs(Ussb))
title('Frequency domain');



%%
%last section:using ideal hilbert

dt = 1/fss:1/fss:length(qu)/(fss);
fc = 102000;

ca1=cos(2*pi*fc*dt);
ca2=sin(2*pi*fc*dt);

x1=ca1.*qu';

hdata=imag(hilbert(data));

uhdata=upsample(hdata,upr);

hlpf=filter(LPF,uhdata);

x2=ca2.*hlpf';

out=x1-x2;

Ussb=fftshift(fft(out,length(out))/fss);
df=-fss/2:fss/length(Ussb):fss/2;
df(length(Ussb)+1)=[];

plot(df,abs(Ussb))
title('Frequency domain');



%%
%Q.3


t=-10:0.01:10;

fm=1;
fc=10;
x=cos(2*pi*fm*t);
figure();
plot(t,x);
title('signal');

c=cos(2*pi*fc*t);
xc=x.*c;
figure();
plot(t,xc);
title('modulated signal');

xr=2*cos(2*pi*fc*t).*cos(2*pi*fc*t).*cos(2*pi*fm*t);
figure();
plot(t,xr);
title('demodulated signal');

out=ifft(fft(2*sinc(2*t)).*fft(xr));

figure();
plot(t,out);
title('filtered signal');

%phase

p=2*pi*rand;
xrp=2*cos(2*pi*fc*t+p).*cos(2*pi*fc*t).*cos(2*pi*fm*t);

outp=ifft(fft(2*sinc(2*t)).*fft(xrp));

figure();
plot(t,outp);
title('output with delay phase');




