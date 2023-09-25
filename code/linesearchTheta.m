function [theta_new,Rnew,Lout,iLineSearch] = linesearchTheta(theta,gradRis,Rprev,Lin,Hdir,H1,H2,S,Nt,maxIter)
Lout = Lin;
for iLineSearch = 1:maxIter
    theta_new = theta+1/Lout*gradRis;
    theta_new = theta_new./abs(theta_new);
    Rnew = computeRate(Hdir,H1,H2,theta_new,S,Nt);

    Qr = Rprev + 2*real(gradRis'*(theta_new-theta)) - Lout*(norm(theta_new-theta))^2;
    if (Rnew>Qr) 
        break;
    else
        Lout = Lout*2;
    end
end

end

