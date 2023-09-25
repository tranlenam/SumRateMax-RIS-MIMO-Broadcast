clear; 
close all; 
clc;
rng(1)
DIR = 1;
Kuser = [3];
Nt = 4;
Nr = 4;
Nris = 100; %225;
lt = 20;
ht = 10;
Pt = 1;
dist_ris = 30;
hris = 5;
N0 = db2pow(-110);

no_mat = 1;

no_iter = 100;

Sin = zeros(Nr,Nr,Kuser);
for ind = 1:length(Kuser)
    tPGM = 0; RPGM = 0;
    tAAO = 0; RAAO = 0;
    tAO = 0; RAO = 0;
    K = Kuser(ind);
    userspec = sprintf('K_%d',K);
    for iChan = 1:no_mat
        
        [Hdir,H1,H2] = generateChannels(Nt,Nr,Nris,lt,ht,1,dist_ris,hris,N0,DIR,K);
        
        % initial RIS phase shifts
        theta0= rand(Nris,1)+1i*rand(Nris,1);
        theta = theta0./abs(theta0);
        
        % initial covariance matrices (randomly generated)
        for user = 1:K
            S0 = rand(Nr,Nr)+1i*rand(Nr,Nr);
            S0 = S0*S0';
            Sin(:,:,user) = S0/trace(S0)*Pt/K;
        end

        % Algorithm 3 (AO)
        [R,t] = Algorithm3AO(Nt,Nris,Pt,K,Hdir,H1,H2,theta,Sin,no_iter);
        tAO = (tAO*(iChan-1)+t)/iChan;
        RAO = (RAO*(iChan-1)+R)/iChan;

        % Algorithm 4 (Approximate AO)
        [R,t] = Algorithm4ApproximateAO(Nt,Nr,Pt,K,Hdir,H1,H2,theta,Sin,no_iter);
        tAAO = (tAAO*(iChan-1)+t)/iChan;
        RAAO = (RAAO*(iChan-1)+R)/iChan;

        % Algorithm 5 (APGM)
        [R,t] = Algorithm5APGM(Nt,Nr,Pt,K,Hdir,H1,H2,theta,Sin,no_iter);
        tPGM = (tPGM*(iChan-1)+t)/iChan;
        RPGM = (RPGM*(iChan-1)+R)/iChan;

    end
    figure
    plot(tPGM,RPGM,'r','DisplayName','PGM');
    hold on;
    plot(tAAO,RAAO,'b','DisplayName','Approx. AO');
    hold on;
    plot(tAO,RAO,'g','DisplayName','AO');
    legend('show','location','southeast');
    title('Convergence vs. wall-clock time');
    saveas(gcf,"../results/Convergence_clocktime.png")
    figure
    plot(RPGM,'r','DisplayName','Algorithm 5 (PGM)');
    hold on;
    plot(RAAO,'b','DisplayName','ALgorithm 4 (Approx. AO)');
    hold on;
    plot(RAO,'g','DisplayName','Algorithm 3 (AO)');
    legend('show','location','southeast');
    title('Convergence vs. iteration');
    saveas(gcf,"../results/Convergence_iterationcount.png")
    
end