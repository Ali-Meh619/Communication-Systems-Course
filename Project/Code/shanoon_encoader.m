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