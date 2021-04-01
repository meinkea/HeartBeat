% v 1.0.2
function P = CoverianceUpdate(K, P_bar_xx, P_bar_hh)

    tmp = K*P_bar_hh*K';
    P = P_bar_xx - K*P_bar_hh*K';
    %P = (P + P')/2;
    
    %P = (eye(n)-K*H')*P_bar_xx*(eye(n)-K*H')'-K*R*K';
    
end