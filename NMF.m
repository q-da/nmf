%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NMF算法 %%%%%%%%%%%%%%%%%%%%%%%%% V = WH
function [W,H]=NMF(V,k,maxiter)
[m,n]=size(V);
W=rand(m,k);
H=rand(k,n);
errs = zeros(maxiter,1);
myeps = 1e-10;
for t = 1:maxiter

   W = W .* ( (V*H') ./ max(W*(H*H'), myeps) ); 
%   W = normalize_W(W,1);
   H = H .* ( (W'*V) ./ max((W'*W)*H, myeps) );

   loss = sum((V-W*H).^2);
   errs(t) = sum(sum(loss));
end


end