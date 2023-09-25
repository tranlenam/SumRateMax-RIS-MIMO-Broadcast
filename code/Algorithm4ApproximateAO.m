function [Rao,t_hyb,theta,Itheta,I,T] = Algorithm4ApproximateAO(Nt,Nr,Pt,K,Hdir,H1,H2,theta,Sin,no_cycle)

% no_cycle = 30;  % normally 300, however for rate vs users 100

Lmax = K*Nt/Pt;
% Lmax = 2;

eta = 1e-6;

t_hyb = zeros(1,2*no_cycle+1);

L20 = 1e-4;
L2=L20;
R = zeros(1,2*no_cycle+1);
R(1) = computeRate(Hdir,H1,H2,theta,Sin,Nt);
Itheta = 0;
I = 0;
T = 0;
tic
for cycle = 1:2:2*no_cycle
    
    Htot = Hdir+pagemtimes(H2,diag(theta)*H1);
    Hsum = eye(Nt)+sum(pagemtimes(pagemtimes(Htot,'ctranspose',Sin,'none'),Htot),3);
    
    % updated covariance matrices
    [Sin] = Algorithm2(Pt,eta,Lmax,Hsum,Htot,Sin);
    t_hyb(cycle+1) = t_hyb(cycle+1)+toc;
    %trace(sum(Sin,3))
    R(cycle+1) = computeRate(Hdir,H1,H2,theta,Sin,Nt);
    
    % RIS Phase Shift optimization
    g_RIS = grad_theta(Htot,H1,H2,Sin,Nt);
    [theta,R(cycle+2),L2,I] = linesearchTheta(theta,g_RIS,R(cycle+1),L2,Hdir,H1,H2,Sin,Nt,30);
    t_hyb(cycle+2) = t_hyb(cycle+2)+toc;
end

Rao = R/log(2);