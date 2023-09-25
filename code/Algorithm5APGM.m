function [Rpgm,tpgm,theta,Is,Itheta] = Algorithm5APGM(Nt,Nr,Pt,K,Hdir,H1,H2,theta,Sin,no_iter)
L10 = 1e-4;
L20 = 1e-4;

tpgm = zeros(1,2*no_iter+1);
Racc = zeros(1,2*no_iter+1);


S = Sin;

L1 = L10;
L2 = L20;
Racc(1) = computeRate(Hdir,H1,H2,theta,S,Nt);

Is = 0; 
Itheta = 0;
tic
for iter = 1:2:2*no_iter

    Htot = Hdir+pagemtimes(H2,diag(theta)*H1);
    Hsum = eye(Nt)+sum(pagemtimes(pagemtimes(Htot,'ctranspose',S,'none'),Htot),3);

    g_S = grad_S(Htot,Hsum);
    RS=computeRate(Hdir,H1,H2,theta,S,Nt);
    [S,Rprev,L1,I]=linesearchS(S,g_S,RS,L1,Hdir,H1,H2,theta,Nt,Nr,Pt,K,30);
    Racc(iter+1) = Rprev;
    tpgm(iter+1) = tpgm(iter+1)+toc;
    Is = Is+I;
    
    Htot = Hdir+pagemtimes(H2,diag(theta)*H1);
    % update RIS phase shifts   
    g_RIS = grad_theta(Htot,H1,H2,S,Nt);
    
    Rtheta = computeRate(Hdir,H1,H2,theta,S,Nt);
    [theta,Rprev,L2,I] = linesearchTheta(theta,g_RIS,Rtheta,L2,Hdir,H1,H2,S,Nt,30);
    Racc(iter+2) = Rprev;
    tpgm(iter+2) = tpgm(iter+2)+toc;
    Itheta = Itheta+I;
end
Rpgm = Racc/log(2);
Is = Is/no_iter;
Itheta = Itheta/no_iter;



