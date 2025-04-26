function loss = CDF_loss_p2(coef, dis_mu, dis_sig)

if coef < dis_mu + 1.96 * dis_sig  && coef > dis_mu - 1.96 * dis_sig
    loss = 0;
else
    loss = (abs(coef - dis_mu) - 1.96 * dis_sig)^2;
end