function [ac,flag]=computer_current_weight(y_XPw,wModel,Size)
 result=zeros(1,y); 
  
    resultpre=zeros(1,y);
    for i=1:y
         wp = zeros(1,Size+1);
         for j=1:(Size+1)
             wp(j)=y_XPw(i,j);
         end            
         result(i)=CompareTwoModels(wp,wModel,Size+1);
         
         if(result(i)<1-result(i))
             resultpre(i)=1-result(i);
         else
             resultpre(i)=result(i);
         end
    end 
    [~,index]=max(resultpre);
    if(result(index)~=resultpre(index))
        if(prediction_array(index)<resultpre(index))
            prediction_array(index)=resultpre(index);
            flag(index)=1;
            wModel=transpose(wModel);
            for j=1:(Size+1)
              wModellist(index,j)=-wModel(j);
            end  
        end
    else
         if(prediction_array(index)<resultpre(index))
            prediction_array(index)=resultpre(index);
            flag(index)=0;
            wModel=transpose(wModel);
            for j=1:(Size+1)
              wModellist(index,j)=wModel(j);
            end  
        end
        
    end
end