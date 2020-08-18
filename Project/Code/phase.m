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