%%Final Project Communication Systems


%% Phase one/Shannon Fano 
clear
clc
image=imread('8.gif');
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

%%
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

%%
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

%%

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



p=1;
o=length(code) ;  
    
i=2;
j=1;

while(o>0)
     t=i;   
        

  for k=1:size(nmb,2)
     
      if string(code(j:i))==string(nmb{1,k})
         
          
         yy(1,p)=i;
         yy(2,p)=j;
         yy(3,p)=nmb{2,k};
         o=o-length(code(j:i));
         
         p=p+1;
         j=i+1;
         i=j+1;
         
         break;
         
         
      end
      
  end
  
  
  if t==i
      
      i=i+1;
  end
  
    end
         

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
 clear
clc
M=5;
N=5;
Sequence=zeros(1,N); %%sequence made
Sequence(1)=1;
Sequence(2)=2;
Sequence(3)=4;
Sequence(4)=3;
Sequence(5)=5;
Ts=2;  %% sending rate
fc=5;
fs=10^5;
[modulated,t]=modulator(M,N,Sequence,Ts,fc,fs);
Sm=baseband(M,Ts,fc,fs);
figure(14);
plot(t,modulated);
N=length(modulated)/(Ts*fs);
st=modulated(1,1:Ts*fs);
c=zeros(M,N);
detector=zeros(1,N);
for m=1:M
    for n=1:N
        c(m,n)=sum(Sm(m,:).*modulated(1,1+(n-1)*Ts*fs:n*Ts*fs));
    end
    
end

for i=1:N
     max_num=max(c(1,i));
     max_index=1;
     for m=2:M
         if max(c(m,i))>max_num
             max_num=max(c(m,i));
             max_index=m;
         end
     end
    
    detector(1,i)=max_index;
end




domain = linspace(-pi,pi,length(st));
f = (domain/(2*pi))*fs;



%% Detector

th=((c(1,1)+c(2,1))/2);
for i=1:N
     if c(1,i)>th
          detector(1,i)=0;
     else
         detector(1,i)=1;
     end
end
    
    
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