function grad_arr = grad_theta(Htot,H1,H2,Y,Nt)
sum_mat = sum(pagemtimes(pagemtimes(H2,'ctranspose',Y,'none'),Htot),3);
Hsum = eye(Nt)+sum(pagemtimes(pagemtimes(Htot,'ctranspose',Y,'none'),Htot),3);
grad_arr = diag(sum_mat*inv(Hsum)*H1');
end