function [Sout,obj_seq_GS,I,T] = Algorithm2(P,eta,Lmax,Hsum,Htot,Sin)
% This function solves (12)

Lmin = 0;
rep = 1;
K = size(Sin,3);
Sum_tr = 0;
for ind = 1:K
    Sum_tr = Sum_tr+trace(Sin(:,:,ind));
end

maxIter=10*K  ; % maximum number of iterations for each bisection step
obj_seq_GS = zeros(maxIter+1,1);
d=zeros(K,1);
for k=1:K
    d(k)=max(real((eig(Htot(:,:,k)*((Htot(:,:,k))')))));
end
T = 0; % number of outer iterations
I = 0; % number of inner iterations
while ((Lmax-Lmin)>=eta)
    L = (Lmin+Lmax)/2;
    % Gauss-Southwell method
    obj_seq_GS_current = log(real(det(Hsum))) - L*(Sum_tr-P);

    %% Algorithm 1 GBCM
    steplength = zeros(K,1);
    for n = 1:maxIter
        Hprod = pagemtimes(Htot,pagemtimes(inv(Hsum),'none',Htot,'ctranspose'));
        gradS = Hprod-L*eye(size(Htot(:,:,ind),1));
        for k=1:K
            steplength(k) = norm(Sin(:,:,k)-ProjPSD(Sin(:,:,k)+1/(d(k)^2)*gradS(:,:,k)),'fro');
        end
        [~,indx] = max(steplength);
        [Hsum,Sin,Sum_tr] = update_S(Hsum,Htot,Sin,L,indx,Sum_tr);
        obj_seq_GS_new = log(real(det(Hsum))) - L*(Sum_tr-P);
        if(abs(obj_seq_GS_new-obj_seq_GS_current)<1e-6) && (n>K)
            break
        end
        obj_seq_GS_current = obj_seq_GS_new;
    end
    %%
    I = I+n;
    T = T+1;
    s = real(trace(sum(Sin,3)));
    if P<s
        Lmin = L;
    else
        Lmax = L;
    end
end
Sout = Sin;
end

function Xplus = ProjPSD(X)
[U,D]=eig(X,'vector');
Xplus = U*diag(max(real((D)),0))*U';
end


function [Hsum_out,Sout,Str_out]= update_S(Hsum_in,Htot,Sin,lambda,ind,Str_in)
% take the channel and the covariance matrix
% solve (15) for k = ind
Hk = Htot(:,:,ind);
Sk = Sin(:,:,ind);
tr1 = sum(diag(Sk));
% calulate the covariance matrix
Hbar = Hsum_in-(Hk')*Sk*Hk;
if size(Hk,1)==1
    d = real(Hk*(Hbar\(Hk')));
    Sk = max(1/lambda-1./d,0);
else
    [V,D] = eig(Hk*(Hbar\(Hk')));
    d = real(diag(D));
    d = max(1/lambda-1./d,0);
    Sk = V*diag(d)*V';
end

% update the array of covariance matrices
Sout = Sin;
Sout(:,:,ind) = Sk;
Hsum_out = Hbar+(Hk')*Sk*Hk;
Str_out = real(Str_in-tr1+sum(d));
end
