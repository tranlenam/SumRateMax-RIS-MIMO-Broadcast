function R = computeRate(Hdir,H1,H2,theta,Sin,Ndim)
Htot = Hdir+pagemtimes(H2,diag(theta)*H1);
Hsum = eye(Ndim)+sum(pagemtimes(pagemtimes(Htot,'ctranspose',Sin,'none'),Htot),3);
R = real(log(det(Hsum)));
end