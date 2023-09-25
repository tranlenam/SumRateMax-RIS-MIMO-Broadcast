function [Hdir,H1,H2] = generateChannels(Nt,Nr,Nris,lt,ht,no_mat,dist_ris,hris,N0,DIR,K,varargin)

Kr = 1;                  % Rician K factor 
f = 2000e6;              % frequency

alpha_dir = 3;
arr = 200:2:500;
arr = 600:2:1000; % added for simulation
dk_all = randi([1 length(arr)],K,no_mat);
dk_all = arr(dk_all)';
dk_all = [100; 100];
% dk_all=[100];
dk_all = [100;100; 100];

lk_all = randi([1 70],K,no_mat);
hk_all = randi([150 200],K,no_mat)/100;

if isempty(varargin)
%     [Hdir,H1,H2] = chan_mat_RIS_noise_new(Nt,Nr,Nris,lt,ht,lk_all,dk_all,hk_all,no_mat,Kr,f,dist_ris,hris,N0,DIR,alpha_dir);
    [Hdir,H1,H2] = chan_mat_RIS_single(Nt,Nr,Nris,lt,ht,lk_all,dk_all,hk_all,no_mat,Kr,f,dist_ris,hris,N0,DIR,alpha_dir);
else
    prob = varargin{1};
    [Hdir,H1,H2] = chan_mat_RIS_prob(Nt,Nr,Nris,lt,ht,lk_all,dk_all,hk_all,no_mat,Kr,f,dist_ris,hris,N0,DIR,prob,alpha_dir);
end 