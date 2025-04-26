clear
addpath(genpath('./puf_util'));
chalSize = 64;
nhighXOR = 1;
nbelowXOR = 3;
nTrS = 20000;
nTeS = 5000;
sigmaNoise = 0.10;

flag_diag = 0; %sep cma-es 1-used
sum_repeat = 25;%repeat num
flag_fitloss = 0;% 1-used
loss_func = 'CDF_loss_p2';
loss_alpha = 1;
dis_mu = 0.20142668480353;
dis_sig = 0.0194236824685995;

flag_earlystop = 1;% 1-used

f_record = "./record/" + flag_diag + flag_fitloss + flag_earlystop ...
    + "_" + chalSize + "_" + nhighXOR+"_"+ nbelowXOR+ "_" + sigmaNoise + "_record.csv";
record = zeros(sum_repeat+1, 4);
for i = 1:sum_repeat
    [acc, ind_time, tot_time, try_num] = iPUF_ATTACK(chalSize,nhighXOR,nbelowXOR,nTrS,nTeS,sigmaNoise, ...
        flag_diag,flag_earlystop, ...
        flag_fitloss,loss_func,loss_alpha,dis_mu,dis_sig);
    record(i, :) = [acc, ind_time, tot_time, try_num];
end
mask = (record(:, 1) ~= 0);
record(end, :) = mean(record(mask, :), 1);

if ~exist('./record', 'dir')
    mkdir('./record');
end

writetable(array2table(record, 'VariableNames', {'Accuracy', 'Ind time', 'Tot time', 'Try num'}), f_record);
