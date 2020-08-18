function S=baseband(M,Ts,fc,fn)
n=Ts*fn;
S=zeros(M,n);
  for m=1:M
    u=sqrt(2*M/(m*Ts));
    t=linspace(0,m*Ts/M,m*Ts*fn/M);
    S(m,1:fn*m*Ts/M)=u*cos(2*pi*fc*t);
  end 
end