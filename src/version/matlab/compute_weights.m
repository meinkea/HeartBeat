% v 1.0.2
function  [Wc,Wm]=  computeWeights(n,lambda,alpha,beta)

Wc = (1/(2*(n + lambda))) * ones(2 * n + 1, 1);
Wm = (1/(2 * (n + lambda))) * ones(2 * n + 1,1);
Wc(1) = lambda / (n + lambda) + (1. - alpha^2 + beta);
Wm(1) = lambda / (n + lambda);

Wc = Wc';
Wm = Wm';

end