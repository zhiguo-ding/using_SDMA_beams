clear all;
M = 10;
N = M;
snrdb = [0: 10: 30 ];
ct = 10000;
Rp = 0.5; epsP = 2^Rp-1;
Rs = 1;   epsS = 2^Rs-1;
Ps = 1; Pp=1;

for i = 1:length(snrdb)
    snr = 10^(snrdb(i)/10);
    sum1 = 0;sum2 = 0;sum3 = 0;sum4 = 0;sum5 = 0;sum6 = 0;
    for ict = 1 : ct
        %primary users' channel vectors
        G = complex(sqrt(0.5)*randn(N,M),sqrt(0.5)*randn(N,M));
        % make G orthogonal
        Gx = complex(sqrt(0.5)*randn(N,M),sqrt(0.5)*randn(N,M));
        [G,v,d] = svd(Gx);
        G = G*diag(complex(sqrt(0.5)*randn(1,M),sqrt(0.5)*randn(1,M)));
        
        D_power = diag(1./[diag(inv(G'*G))]);        
        F = G*inv(G'*G)*sqrt(D_power)/sqrt(M); %normalized digital beamforming 
     
        h = complex(sqrt(0.5)*randn(N,1),sqrt(0.5)*randn(N,1));
        
        gm = abs(diag(G'*F)).^2; %primary users' effective channel gains
        hm = abs(F'*h).^2; %secondary user's effective channel gains on beams
        
        alphamPI = min(1,epsP/snr./gm);
        alphaSIIx1 = max(0,(gm-epsP/snr)./gm/(1+epsP));        
        for m = 1 : M
            him = hm; him(m)=[];
            alphai = alphamPI; alphai(m) = [];
            temp1 = sum(him.*alphai);
            alphaSIIx2(m,1) = max(0, (hm(m)-epsP*temp1-epsP/snr)/hm(m)/(1+epsP));
        end
        alphaSII = min(alphaSIIx1,alphaSIIx2);
        alphaPII = 1 - alphaSII;
             
        %direct implmentation 
        temp2=[];temp3=[];
        for m = 1 : M
            him = hm; him(m)=[];
            alphai = alphamPI; alphai(m) = [];
            temp1 = sum(him.*alphai);
            temp2(m,1) = hm(m)*alphaPII(m)/(hm(m)*alphaSII(m)+temp1+1/snr);
            temp3(m,1) = hm(m)*alphaSII(m)/(temp1+1/snr);
        end
        stateI = alphaSII>0;%temp2>=epsP;
        stateII = temp3>=epsS;
        state = stateI.*stateII;
        if max(state)>0 %one of the beams is successful
            sum1 = sum1 +1;
        end
        
        temp4=[];
        %probability way, it provides the same outage probaiblity as sum1
        %but easy for analysis as well as ergodic rate
        for m = 1 : M
            him = hm; him(m)=[];
            alphai = alphamPI; alphai(m) = [];
            temp1 = sum(him.*alphai);
            if alphaSII(m)>0
                temp4 = [temp4 hm(m)*alphaSII(m)/(temp1+1/snr)];
            end
        end
        if isempty(temp4)
            sum2 = sum2 +1;
        elseif max(temp4)<epsS
            sum2 = sum2 + 1;
        end
        
        %rate
        if isempty(temp4) %no beam  can be used or no power can be assingd
            ratex1(ict) = 0;
        else
            ratex1(ict) = log2(1+max(temp4));
        end
        
     end
    pout1(i) = 1-sum1/ct;
    poutx1(i) = 1-sum3/ct;
    pout2(i) = sum2/ct;
    Rate1(i) = mean(ratex1);

end

%semilogy(snrdb,pout1)%,snrdb,pout4)
plot(snrdb,Rate1);%,snrdb,Rate2)
