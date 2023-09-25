function [theta_out,T] = optimizeRIS(Hdirall,H1,H2all,theta,lind,Sin,Ndim,Kuser)
% solve (21) when l = lind
A = eye(Ndim);
B = 0;
theta_out = theta;

Htot = Hdirall+pagemtimes(H2all,H1.*theta);
h1kl = H1(lind,:);
for user = 1:Kuser
    S = Sin(:,:,user);
    h2kl = H2all(:,lind,user);  
    C = Htot(:,:,user)-h2kl*theta(lind)*h1kl;
    A = A + C'*S*C +h1kl'*h2kl'*S*h2kl*h1kl;
    B = B + C'*S*h2kl;
end
d = h1kl*(A\B);
theta_out(lind) = exp(-1i*angle(d));
T=0;
end

