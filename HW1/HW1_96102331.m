%%Q.1
%random_pulse
len=1;
sigma=80;
fs=1000;
T=0.1;
n=2000;

[r0,r1]=random_pulse(1000,1,0.1,1);
plot(r0,r1);
title('random digital wave')



for i = 1 : n
    
    q0 = sin(400*pi*[0:1/fs:1]);
    [x,y] = random_pulse(1000, 1, 0.1, 1);
    v= q0 .*y;
    
    Yr(i,:)=xcorr(y);%autocorrelation
    Yr_s(i,:)=xcorr(v);
      
    
    Yf(i,:)=fftshift(fft(y,length([-1:0.001:1])));%fourie
    Yf_s(i,:)=fftshift(fft(v,length([-1:0.001:1])));
    
    S2y(i,:)=Yf(i,:).*conj(Yf(i,:));%spectral density
    S2y_s(i,:)=Yf_s(i,:).*conj(Yf_s(i,:));
    
    
end


yr_m=mean(Yr);
dx=-1:1/fs:1;
(dx,yr_m);
title('Auto correlation with T=0.3,sigma=1')


S1y_m=fftshift(fft(yr_m,length([-1:0.001:1])));
df=-1:0.001:1;
plot(df,abs(S1y_m));
title('Spectral density');

S2y_m=mean(S2y);
plot(df,S2y_m);
title('Spectral density');

MSE=immse(abs(S1y_m),S2y_m);

Yr_sm=mean(Yr_s);
plot(dx,Yr_sm);
title('Auto correlation of modulated signal')

S2y_sm=mean(S2y_s);
plot(df,S2y_sm);
title('Spectral density of modulated signal')


%%
%Q.2

%spectral density of Audio

[d,F]= audioread('clip.wav');

Yf=fftshift(fft(d,length([-1:0.001:1])));

Sy=Yf.*conj(Yf);

df=-1:0.001:1;

plot(df,Sy);

title('Spectral density');

%%
%LTI or not???

t=[0:0.2:2*pi];

X0=sin(t);
Y0=cos(t);

I1=Channel(2*X0-Y0);
O1=2*Channel(X0)-Channel(Y0);

plot(t,I1);
hold on
plot(t,O1);
title('LTI or not??')
%%
%LTI or not???

I2=Channel(X0+3*Y0);
O2=Channel(X0)+3*Channel(Y0);

plot(t,I2);
hold on
plot(t,O2);
title('LTI or not??')
%%
%LTI or not???

t=-2*pi:0.2:2*pi;

x1=cos(t);
y1=Channel(x1);

x2=cos(t-pi/2);
y2=Channel(x2);

plot(t,y1,'b');
hold on
plot(t,y2,'r--');
title('LTI or not??')
%%
%LTI or not???

x3=sin(t);
y3=Channel(x3);

x4=sin(t-pi/4);
y4=Channel(x4);
%%
%Impulse Response of channel

ii = dirac([0:1:1000]);
ii(1)=1;
IM=Channel(ii);%Impulse Response


H=fftshift(fft(IM,length([-1:0.001:1])));
df=-1:0.001:1;

plot(df,abs(H));
title('Impulse Response')
ylabel('Magnitude')
figure();
plot(df,angle(H));
title('Impulse Response')
ylabel('Phase')
%%
%group delay and phase delay

Ph=angle(H);

PhD= -Ph./(length([-1:0.001:1]));
GrD= -diff(Ph)/(0.001);

df=-1:0.001:1-0.001;
plot(df,GrD);
title('Group delay')


%%
%spectral density for channel output

y2=Channel(d);

Y2f=fftshift(fft(y2,length([-1:0.001:1])));
S2y=Y2f.*conj(Y2f);

plot(df,S2y);
title('Spectral density')
%%
%channel estimation and Equalizer

[d,fs]= audioread('clip.wav');
o=Channel(d);
data = iddata(o,d,1/fs); 
sys = tfest(data,2);
equ=1/sys;

pzmap(equ);
figure();
plot(-0.125:0.000125:0.125,1./abs(H));
title('Impulse Response of equalizer')
ylabel('Magnitude')
figure();
plot(df,-angle(H));
title('Impulse Response of equalizer')
ylabel('Phase')
%%
%non-causal Equalizer
a=(1./abs(H));
b=exp(-i.*angle(H));
Hequ=a.*b;
Of=fftshift(fft(o,length([-1:0.001:1]))).';
Dequ=Hequ.*Of;
d2=ifftshift(ifft(Dequ));


%%
%Q.3
%Impulse response
fs=22050;
t=(-0.1:1/fs:0.1);
imp=t==0;

IM=imag(hilbert(imp));
plot(t,IM)
%%
%normal_window

M=1000;
q=length(t);

wind=zeros(1,length(t));
wind(((q+1)/2)-M/2:((q+1)/2)+M/2)=1;
h1=IM.*wind;


H1(1:M+1)=h1(((q+1)/2)-M/2:((q+1)/2)+M/2);
Hf1=fftshift(fft(H1,length([-1:0.001:1])));

df=-1:0.001:1;
plot(df,abs(Hf1));
title('Impulse Response for M=1000')
ylabel('Magnitude')
figure();
plot(df,unwrap(angle(Hf1)));
title('Impulse Response for M=1000')
ylabel('Phase')
%%
%hamming_window

M=1000;
q=length(t);
dn2=[0:1:M];

wind2(1:M+1)=0.54-0.46*cos(2*pi*dn2/M);
H2(1:M+1)=IM(((q+1)/2)-M/2:((q+1)/2)+M/2);
H2=H2.*wind2;

H2f=fftshift(fft(H2,length([-1:0.001:1])));
df=-1:0.001:1;

plot(df,abs(H2f));
title('Impulse Response for M=1000')
ylabel('Magnitude')
figure();
plot(df,unwrap(angle(H2f)));
title('Impulse Response for M=1000')
ylabel('Phase')
%%
%kaiser window

M=1000;
q=length(t);
alpha = 0.1;

for i = 1:M+1
   wind3(i)= besseli(0, pi*alpha*sqrt(1-(2*i/M)^2))/besseli(0, pi*alpha);
end

H3(1:M+1)=IM(((q+1)/2)-M/2:((q+1)/2)+M/2);
H3=H3.*wind3;

H3f=fftshift(fft(H3,length([-1:0.001:1])));
df=-1:0.001:1;

plot(df,abs(H3f));
title('Impulse Response for M=1000')
ylabel('Magnitude')
figure();
plot(df,unwrap(angle(H3f)));
title('Impulse Response for M=1000')
ylabel('Phase')

%%
function [x,y] = random_pulse(fs, sigma, T, L)
D=rand*T;

k=ceil((L-D)/T);

x=0:1/fs:L;

y=zeros(1,length(x));

for i=1:k
    
    y(1,floor(fs*D)+floor((i-1)*T*fs)+1:floor(D*fs)+floor(i*T*fs))=normrnd(0,sqrt(sigma));
    
end
a=length(y);
y(length(x)+1:a)=[];

end