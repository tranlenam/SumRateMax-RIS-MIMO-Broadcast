function Sp = proj_S(S,Nr,P,K)
d = [];
V = zeros(Nr,Nr,K);
for user = 1:K
    [V(:,:,user),D] = eig(S(:,:,user));
    d = [d; diag(D)];
end
% water filling
d_out = pow_proj(P,d);
Sp = zeros(Nr,Nr,K);
for user = 1:K
    idx = (user-1)*Nr+1:user*Nr;
    Sp(:,:,user) = V(:,:,user)*diag(d_out(idx))*(V(:,:,user))';
end
    
end