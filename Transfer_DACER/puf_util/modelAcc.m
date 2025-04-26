function f = modelAcc(w,APhi,AResponse,Size,nRows,nXOR)

    %epsilion = 0;
    %epsilion = w(end);
    %w = w(1:end-1);
    Res = 0;
    f=0;
    R=0;
    for i=1:nRows
        %Compute Delay diffference
        for k=1:nXOR
            for j=1:Size
                Res = Res + w(k,j)*APhi(i,j);
            end

            if (Res>0)
                Res=0;
            else
                Res=1;
            end 
                R=R+Res;
        end
        R=mod(R,2);
        %Compute response

        %Compare responese 
        if (R==AResponse(i))
            f=f+1;
        end 
        R=0;
    end

    f = f/nRows;
end