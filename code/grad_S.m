function grad_mat = grad_S(Htot,Hsum)
% grad_mat = [];
% for user = 1:K
%     grad_mat = blkdiag(grad_mat, Htot{user}*inv(Hsum)*(Htot{user})'); 
% end 
grad_mat = pagemtimes(Htot,pagemtimes(inv(Hsum),'none',Htot,'ctranspose'));

end