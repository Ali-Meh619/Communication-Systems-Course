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