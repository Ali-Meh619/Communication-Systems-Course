%% Final Project Communication Systems


%% Phase two:Encoder
clear
clc
image=imread('12.gif');
resizedimage=imresize(image,1/8);
probability=zeros(1,256);
probability(2,1:256)=1:256;

for i=1:size(resizedimage,1)
    for j= 1:size(resizedimage,2)
        a=resizedimage(i,j);
        probability(1,a+1)=probability(1,a+1)+1;
    end
end

[temp, order] = sort(probability(1,:),'descend');
probabilitysorted =probability(:,order);
Zeroprobability=0;

for i=1:256

    probabilitysorted(3,i)=1;
end

tree=cell(1);
leaf=cell(1);
tree={probabilitysorted};


q=3;

nk=1;
qq=1;

for j=1:size(probabilitysorted,2)
    
    c=0;
    d=0;

   
    for i=1:size(tree,2) 
        
        x=tree{j,i};
        
        if size(x,2)>1
        
  
        [y,z]=spilit(tree{j,i});
        
        tree{j+1,nk}=y;
        tree{j+1,nk+1}=z;
        
        nk=nk+2;
  
        elseif size(x,2)==1
     
        c=c+1;
        
        leaf{1,qq}=x;
        
        qq=qq+1;
        
        elseif size(x,2)==0
            
            d=d+1;
        
        end
    
    end 
    
    if (d+c)==size(tree,2) && j~=1
    
        break;     
        
    end
    
    
    nk=1;
    
end   


for i=1:size(leaf,2)
   x=leaf{1,i};
   
  tt(1:2,i)=x(2:3,1);   
end
[temp, order] = sort(tt(1,:));
xx=tt(:,order);



%%
codeword='';
for i=1:size(resizedimage,2)
    for j=1:size(resizedimage,1)
        [r,c]=find(xx==resizedimage(j,i));
        a=num2str(xx(r+1,c));
        a=a(2:end);
        codeword=strcat(codeword,a);
    end
end

code=codeword;


%% Modulator
clc
codeword=code;
fs=10^5;
Ts=0.03;
fc=10^4;
M=2;
N=length(codeword);
Sequence=zeros(1,N); %%sequence made
for i=1:N
    if codeword(1,i)=='0'
        Sequence(1,i)=1;
    else
        Sequence(1,i)=2;
    end
end
[modulated,t]=modulator(M,N,Sequence,Ts,fc,fs);

Sm=baseband(M,Ts,fc,fs);
a=modulated;
%% Part t
clc 
noise=zeros(1,length(modulated));
var=800;
noise=randn(1,length(modulated))*sqrt(var);
modulatorafternoise=modulated+noise;
a=modulatorafternoise;

    %%
    clc
    load('BP2.mat');
    
   s0=filter(Bandpass2,1,Sm(1,:));
   s1=filter(Bandpass2,1,Sm(2,:));
   channel=filter(Bandpass2,1,a);
   idealchannel=filter(Bandpass2,1,modulated);
   
   %% part P
    clc
     fourier=fftshift(fft(modulated));
     f=linspace(-fs/2,fs/2,length(fourier));
     plot(f,abs(fourier));
      
     clc
     S=fourier.*conj(fourier);
     Etotal=sum(S);
     l=length(fourier);
     middle=(l)/2;
     Etest=S(middle);
     for i=0:(l/2)-1
         Etest=Etest+S(middle+i)+S(middle-i);
         if(Etest>0.99*Etotal)
             breakingpoint=i;
             break
         end
     end
 fm=fs*breakingpoint/l;
   %%
   clc
N=length(channel)/(Ts*fs);
st=channel(1,1:Ts*fs);
c=zeros(M,N);
detector=zeros(1,N);

    for n=1:N
        down=1+floor((n-1)*Ts*fs);
        up=length(s0)+floor((n-1)*Ts*fs);
        c(1,n)=sum(s0.*channel(1,down:up));
        c(2,n)=sum(s1.*channel(1,down:up));
    end
    


th=((c(1,1)+c(2,1))/2);
for i=1:N
     if c(1,i)>th
          detector(1,i)=0;
     else
         detector(1,i)=1;
     end
end

codenew='';
for i=1:N
    if detector(1,i)==0
        codenew=strcat(codenew,'0');
    else
        codenew=strcat(codenew,'1');
    end
end

%% Part G
figure(11);
scatter(c(1,:),c(2,:));

%% Error
error=sum(abs(code-codenew));
%% SNR
clc
Pt=sum(modulated.*modulated);
Pnt=sum(noise.*noise);
SNRt=Pt/Pnt
Pr=sum(idealchannel.*idealchannel);
Pnr=sum(channel.*channel)-Pr;
SNRb=Pr/Pnr
Sequencenew=detector+1;
[demodulated,t]=modulator(M,N,Sequencenew,Ts,fc,fs);
demodulated=filter(Bandpass2,1,demodulated);
Pl=sum(demodulated.*demodulated);
Pnde=Pr-Pl;
SNR=Pl/Pnde
%% Decoder
code=codenew;


nmb=cell(2,size(xx,2));

j=1;
for i=1:size(xx,2)
     
      x=num2str(xx(2,i));
      x=x(2:end);
      nmb{1,j}=x;
      nmb{2,j}=xx(1,i);
      j=j+1;
         
end
%% Decoder
code=codenew;


nmb=cell(2,size(xx,2));

j=1;
for i=1:size(xx,2)
     
      x=num2str(xx(2,i));
      x=x(2:end);
      nmb{1,j}=x;
      nmb{2,j}=xx(1,i);
      j=j+1;
         
end
%%

nn=string(code);

kk=nmb{1,1};

up=length(code);




o=length(code) ;  
    
yy = zeros(1,size(resizedimage,1)*size(resizedimage,2));


for i = 1: size(resizedimage,1)*size(resizedimage,2)
    for j = 1:size(nmb,2)
        
        if length(code) >= length(nmb{1,j})
        if string(code(1: length(nmb{1,j})))==string(nmb{1,j})
            yy(i) = nmb{2,j};
            code(1: length(nmb{1,j})) =[];
            break;
        end
        end
    end
end


     %%   
     clc
    imagekk=reshape(yy,[size(resizedimage,1),size(resizedimage,2)]);
    jop=uint8(imagekk);
    imshow(jop);
%%

function codeword=shannon_encoder(begin_point,end_point,p,code)
 for i=begin_point:end_point
     if sum(p(begin_point:i))>sum(p(i+1:end_point))
         break;
     end
 end
 code(begin_point:i)=10*code(begin_point:i)+0;
 code(i+1:end_point)=10*code(i+1:end_point)+1;
 high_point=i;
 low_point=i+1;
  if ((low_point+1)<end_point)
      codeword=code;
     code=shannon_encoder(low_point,end_point,p,code);
    if ((begin_point+1)<high_point)
     codeword=code;
     code=shannon_encoder(begin_point,high_point,p,code);
    end
  end
  if ((begin_point+1)<high_point)
     codeword=code;
     code=shannon_encoder(begin_point,high_point,p,code);
  end
  if ((begin_point+1)==high_point)
      code(begin_point)=10*code(begin_point)+0;
      code(high_point)=10*code(high_point)+1;
      codeword=code;
      return
  end
  if (end_point==(low_point+1))
     code(low_point)=10*code(low_point)+0;
     code(end_point)=10*code(end_point)+1;
     codeword=code;
     return
 end
 
 if ((end_point==low_point)&&(begin_point==high_point))
     codeword=code;
     return
 end
 
end
function [modulated,tTotal]=modulator(M,N,sequence,Ts,fc,fn)
n=floor(N*Ts*fn);
tTotal=linspace(0,N*Ts,N*Ts*fn);
modulated=zeros(1,n);
  for m=1:N
    u=sqrt(2*M/(sequence(m)*Ts));
    t=linspace(0,sequence(m)*Ts/M,sequence(m)*Ts*fn/M);
    modulated(1,(1+fix((m-1)*Ts*fn)):(length(t)+fix((m-1)*Ts*fn)))=u*cos(2*pi*fc*t);
  end 
end

function [modulated,tTotal]=phase(M,N,sequence,Ts,fc,fn,phi)
n=floor(N*Ts*fn);
tTotal=linspace(0,N*Ts,N*Ts*fn);
modulated=zeros(1,n);
  for m=1:N
    u=sqrt(2*M/(sequence(m)*Ts));
    t=linspace(0,sequence(m)*Ts/M,sequence(m)*Ts*fn/M);
    modulated(1,(1+fix((m-1)*Ts*fn)):(length(t)+fix((m-1)*Ts*fn)))=u*cos(2*pi*fc*t+phi);
  end 
end


function S=baseband(M,Ts,fc,fn)
n=Ts*fn;
S=zeros(M,n);
  for m=1:M
    u=sqrt(2*M/(m*Ts));
    t=linspace(0,m*Ts/M,m*Ts*fn/M);
    S(m,1:fn*m*Ts/M)=u*cos(2*pi*fc*t);
  end 
end


function [y,z]=spilit(x)

for i=1:size(x,2)
    
     if sum(x(1,1:i))>=sum(x(1,i+1:end))
         
         y=ones(3,i);
         y(1:2,:)=x(1:2,1:i);
         y(3,:)=10*x(3,1:i)+0;
         
         z=ones(3,size(x,2)-i);
         z(1:2,:)=x(1:2,i+1:end);
         z(3,:)=10*x(3,i+1:end)+1;
         
         break;
     end
     end
end