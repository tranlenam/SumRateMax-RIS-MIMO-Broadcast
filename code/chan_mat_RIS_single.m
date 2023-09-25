function [Hdir,H1,H2] = chan_mat_RIS_singe(Nt,Nr,Nris,lt,ht,lk_all,dk_all,hk_all,no_mat,K,f,dist_ris,hris,N0,DIR,varargin)

lambda = 3e8/f;     % wavelength
dt = lambda/2;      % TX antenna space
dr = lambda/2;      % RX antenna space
dris = lambda/2;    % RIS element space
k = 2*pi/lambda;    % wavenumber
Kuser = size(lk_all,1); % numebr of users

% geometrical placement
% x, y and z axis 
% TX array
tx_arr = zeros(1,Nt); 
tx_arr(2,:) = (sort(0:Nt-1,'descend')-(Nt-1)/2)*dt+lt; 
tx_arr(3,:) = ht; 

% RIS
center = [dist_ris 0]; % RIS center position 
N1 = sqrt(Nris);
N2 = N1;
ris_pos = RISPosition(N1,N2,dris,center); % RIS elements coordinates 
a = repmat(ris_pos{1},N1,1);   % placing RIS elements in proper coordinates
ris_arr(1,:) = a(:)';        
ris_arr(2,:) = zeros(1,Nris);
ris_arr(3,:) = repmat(ris_pos{2},1,N2)+hris; 

if isempty(varargin)
    alpha = 2;
else
    alpha = varargin{1};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TX-RIS
for l1 = 1:Nris
    for t1 = 1:Nt
        d1(l1,t1) = norm(tx_arr(:,t1)-ris_arr(:,l1));
    end
end
% TX-RIS
dt_ris = sqrt(dist_ris^2+lt^2+(ht-hris)^2);  % TX-RIS distance
H1_los = exp(-1i*k*d1);
FSPL_1 = 1;
H1 = Rician(H1_los,FSPL_1,no_mat,K,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hdir = zeros(Nr,Nt,Kuser);
H2 = zeros(Nr,Nris,Kuser);
for user = 1:Kuser
    lk = lk_all(user);
    dk = dk_all(user);
    hk = hk_all(user);
    
    % RX array postion
    rx_arr(1,:) = dk*ones(1,Nr);
    rx_arr(2,:) = (sort(0:Nr-1,'descend')-(Nr-1)/2)*dr+lk;
    rx_arr(3,:) = hk;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % direct TX-RX paths/channel matrix
    for i1 = 1:Nr
        for j1 = 1:Nt
            d(i1,j1) = norm(rx_arr(:,i1)-tx_arr(:,j1));
        end
    end
    % Ditect link matrix
    Hdir_los = exp(-1i*k*d);
    dt_k = sqrt(dk^2+(lt-lk)^2+(ht-hk)^2);                    % TX-RX distance
    FSPL_dir = (lambda/(4*pi))^2/dt_k^alpha(1)*DIR;
    Hdir(:,:,user) = Rician(Hdir_los,FSPL_dir,1,K,N0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % indirect paths (TX-RIS-RX)
    for l1 = 1:Nris
        for r1 = 1:Nr
            d2(r1,l1) = norm(rx_arr(:,r1)-ris_arr(:,l1));
        end
    end
    
    dris_k = sqrt((dist_ris-dk)^2+lk^2+(hris-hk)^2);    % RIS-RX distance
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    Gt = 2; Gr = 2;
    cosThetaI = lt/dt_ris; cosThetaR = lk/dris_k;
    FSPLindir = Gt*Gr*lambda^4/(256*pi^2)*...
        cosThetaI*cosThetaR/(dt_ris*dris_k)^2;
%     if FSPLindir == 0
%         FSPLindir
%         pause       
%     end
    % RIS-RX
    H2_los = exp(-1i*k*d2);
    FSPL_2 = FSPLindir;
    H2(:,:,user) = Rician(H2_los,FSPL_2,1,K,N0);
end
end

function pos = RISPosition(N1,N2,dist,center)
d1 = (0:N1-1)-(N1-1)/2;
d2 = (0:N2-1)-(N2-1)/2;
pos{1} = center(1)+d1*dist;
pos{2} = center(2)+d2*dist;
end 

% Rician channel
function Hout = Rician(Hlos,FSPL,no_mat,K,N0)
Hlos = repmat(Hlos,no_mat,1);
Hnlos = sqrt(1/2)*(randn(size(Hlos))+1i*randn(size(Hlos)));
Htot = sqrt(FSPL)/sqrt(K+1)*(Hlos*sqrt(K)+Hnlos);
% adding the noise power
Htot = Htot/sqrt(N0);
%%%%%%%%%%%%%%%
dim = size(Hlos,1)/no_mat;
if no_mat>1
    for ind = 1:no_mat
        Hout{ind} = Htot((ind-1)*dim+1:ind*dim,:);
    end
else
    Hout = Htot;
end
end