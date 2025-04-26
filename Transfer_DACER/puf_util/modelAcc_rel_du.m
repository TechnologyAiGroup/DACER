

function f = modelAcc_rel_du(w,x,h)

    %epsilion = w(end);
    %w = w(1:end-1);
    
    P = x * w;
    P = abs(P);
    h_tilde = P;
    %h_tilde = P > abs(epsilion);
    %h_tilde = double(h_tilde);
    
    %r = corrcoef(h,h_tilde,'rows','complete');
    r = corrcoef(h,h_tilde,'rows','complete');
    %f = -1*r(1,2);
    f = r(1,2);
end
