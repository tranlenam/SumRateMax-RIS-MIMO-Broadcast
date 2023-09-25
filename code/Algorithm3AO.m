function [Rao,t_ao,theta,Iiter,Titer] = Algorithm3AO(Nt,Nris,Pt,K,Hdir,H1,H2,theta,Sin,no_cycle)

Lmax = K*Nt/Pt;
% Lmax=2;
eta = 1e-6;

t_ao = zeros(1,2*no_cycle+1);
R = zeros(1,2*no_cycle+1);

R(1) = computeRate(Hdir,H1,H2,theta,Sin,Nt);
Iiter = 0;
Titer = 0;
tic
for cycle = 1:2:2*no_cycle
    % compute total user matrices and sum of products Htot{user}'*Sin{user}*Htot{user}
    Htot = Hdir+pagemtimes(H2,diag(theta)*H1);
    Hsum = eye(Nt)+sum(pagemtimes(pagemtimes(Htot,'ctranspose',Sin,'none'),Htot),3);
    
    % updated covariance matrices
    [Sin,~,deltaI,deltaT] = Algorithm2(Pt,eta,Lmax,Hsum,Htot,Sin);
    R(cycle+1) = computeRate(Hdir,H1,H2,theta,Sin,Nt);
    t_ao(cycle+1) = t_ao(cycle+1)+toc;
     
    % Optimize RIS phase shifts sequentially
    for ris_ind = 1:Nris
        [theta,~] = optimizeRIS(Hdir,H1,H2,theta,ris_ind,Sin,Nt,K);
        
    end
    R(cycle+2) = computeRate(Hdir,H1,H2,theta,Sin,Nt);
    t_ao(cycle+2) = t_ao(cycle+2)+toc;
    
    Iiter = Iiter+deltaI;
    Titer = Titer+deltaT;
end
Rao = R/log(2);
Iiter = Iiter/no_cycle;
Titer = Titer/no_cycle;
