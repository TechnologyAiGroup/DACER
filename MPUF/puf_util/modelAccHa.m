function f = modelAccHa(w,APhi,AResponse,Size,nRows)

    %epsilion = 0;
    %epsilion = w(end);
    %w = w(1:end-1);
    Res = 0;
    f=0;
    
    for i=1:nRows
        %Compute Delay diffference
        for j=1:Size
          Res = Res + w(j)*APhi(i,j);
        end
        
        %Compute response
        if (Res>0)
            Res=0;
        else
            Res=1;
        end 
        
        %Compare responese 
        if (Res==AResponse(i))
            f=f+1;
        end 
        
    end
    f = f/nRows;
end