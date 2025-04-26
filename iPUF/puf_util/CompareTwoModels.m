function count = CompareTwoModels(w1,w2,Size)
    
    nTest  = 10000;
    chalSize = Size-1;
    Test= randi([0 1], nTest, chalSize);
    APhi = Transform(Test, nTest, chalSize);
    
    Res1 = ComputeResponseAPUF(w1,APhi,nTest,Size);
    Res2 = ComputeResponseAPUF(w2,APhi,nTest,Size);    
  
    count = 0;
    for i=1:nTest
      if(Res1(i)==Res2(i))
         count = count+1;
      end    
    end
    
    count = count/nTest;
    
end