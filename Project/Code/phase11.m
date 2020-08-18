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
if probabilitysorted(1,i)==0
    probabilitysorted(3,i)=9;
    Zeroprobability=Zeroprobability+1;
else
    probabilitysorted(3,i)=1;
end
end

probabilitysorted(:,size(probabilitysorted,2)-Zeroprobability+1:end)=[];

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
[modulated1,t]=phased(M,N,Sequence,Ts,fc,fs,pi/3);
[modulated2,t]=phased(M,N,Sequence,Ts,fc,fs,pi/6);
[modulated3,t]=phased(M,N,Sequence,Ts,fc,fs,pi/12);

Sm=baseband(M,Ts,fc,fs);

%% Part t
clc 
noise=zeros(1,length(modulated1));
var=100;
noise=randn(1,length(modulated1))*sqrt(var);
modulatorafternoise1=modulated1+noise;
modulatorafternoise2=modulated2+noise;
modulatorafternoise3=modulated3+noise;

    %%
    clc
    load('BP2.mat');
    
   s0=filter(Bandpass2,1,Sm(1,:));
   s1=filter(Bandpass2,1,Sm(2,:));
   channel1=filter(Bandpass2,1,modulatorafternoise1);
   channel2=filter(Bandpass2,1,modulatorafternoise2);
   channel3=filter(Bandpass2,1,modulatorafternoise3);
   

   %%
   clc
N=length(channel1)/(Ts*fs);
st=channel1(1,1:Ts*fs);
c1=zeros(M,N);
c2=zeros(M,N);
c3=zeros(M,N);
detector1=zeros(1,N);
detector2=zeros(1,N);
detector3=zeros(1,N);

    for n=1:N
        down=1+floor((n-1)*Ts*fs);
        up=length(s0)+floor((n-1)*Ts*fs);
        c1(1,n)=sum(s0.*channel1(1,down:up));
        c1(2,n)=sum(s1.*channel1(1,down:up));
        c2(1,n)=sum(s0.*channel2(1,down:up));
        c2(2,n)=sum(s1.*channel2(1,down:up));
        c3(1,n)=sum(s0.*channel3(1,down:up));
        c3(2,n)=sum(s1.*channel3(1,down:up));
    end
    


th=((c1(1,1)+c1(2,1))/2);
for i=1:N
     if c1(1,i)>th
          detector1(1,i)=0;
     else
         detector1(1,i)=1;
     end
end

th=((c2(1,1)+c2(2,1))/2);
for i=1:N
     if c2(1,i)>th
          detector2(1,i)=0;
     else
         detector2(1,i)=1;
     end
end
th=((c3(1,1)+c3(2,1))/2);
for i=1:N
     if c3(1,i)>th
          detector3(1,i)=0;
     else
         detector3(1,i)=1;
     end
end
codenew1='';
for i=1:N
    if detector1(1,i)==0
        codenew1=strcat(codenew1,'0');
    else
        codenew1=strcat(codenew1,'1');
    end
end

%% Part G
figure(11);
scatter(c1(1,:),c1(2,:));
hold on
scatter(c2(1,:),c2(2,:));
hold on
scatter(c3(1,:),c3(2,:));
     %%   

    now=zeros(64,64);
for i=1:64
    for j=1:64
        n=(i-1)*64+j;
        if n<=length(yy(3,:))
        now(j,i)=yy(3,n);
        end
    end
end
now=mat2gray(now);
imshow(now);
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

function [modulated,tTotal]=phased(M,N,sequence,Ts,fc,fn,phi)
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
