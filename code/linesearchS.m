function [Snew,Rnew,Lout,iLineSearch] = linesearchS(S,gradS,Rprev,Lin,Hdir,H1,H2,theta,Nt,Nr,Pt,K,maxIter)
Lout = Lin;
Lmax = 1e6; % the maximum value of Lipschitz constant of gradS
for iLineSearch = 1:maxIter
    Snew = S+1/Lout*gradS;
    Snew = proj_S(Snew,Nr,Pt,K);
    Rnew = computeRate(Hdir,H1,H2,theta,Snew,Nt);
    Sdiff = Snew-S;
    Sprod = pagemtimes(gradS,Sdiff);
    trsum = trace(sum(Sprod,3));
    Qs = real(Rprev + trsum - Lout/(2)*(norm(Sdiff(:)))^2);
    if (Rnew>Qs) || (Lout > Lmax)
        break;
    else
        Lout = Lout*2;
    end
end
end

