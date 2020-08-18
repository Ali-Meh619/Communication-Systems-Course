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