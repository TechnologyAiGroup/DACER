function [testchallenge,testresponse]=gettest_set(x_xpw,chalSize,x)
   nTrS = 20000; %size of test set
   testchallenge= randi([0 1], nTrS, chalSize);
   Phi_TrSp = Transform(testchallenge, nTrS, chalSize);
   testresponse=ComputeResponseXOR(x_xpw,x,Phi_TrSp,nTrS,size(Phi_TrSp,2));
end