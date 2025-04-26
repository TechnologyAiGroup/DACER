function f = modelAccHa(w,APhi,Size,nRows,flag)

    %epsilion = 0;
    %epsilion = w(end);
    %w = w(1:end-1);
    Res = 0;
    f=zeros(1,nRows);
    if(flag==0)
        for i=1:nRows
        %Compute Delay diffference
        for j=1:Size
          Res = Res + w(j)*APhi(i,j);
        end
        
        %Compute response
        if (Res>0)
            f(i)=0;
        else
            f(i)=1;
        end 
        end
    else
        for i=1:nRows
        %Compute Delay diffference
        for j=1:Size
          Res = Res + w(j)*APhi(i,j);
        end
        
        %Compute response
        if (Res>0)
            f(i)=1;
        else
            f(i)=0;
        end 
        end
    end
        
end