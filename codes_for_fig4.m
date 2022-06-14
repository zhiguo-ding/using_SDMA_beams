clear all;
M = 10;
N = M;
snrdb = [0: 10: 30];
ct = 1000;
Rp = 0.5; epsP = 2^Rp-1;
Rs = 1;   epsS = 2^Rs-1;
Ps = 1; Pp=1;

for i = 1: length(snrdb)
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
        %order two vectors according to hm, descending
        [temp1,hmindex]=sort(hm,'descend');
        gm = gm(hmindex);
        hm = hm(hmindex);

        for Dsize = 1 : M
                   
            for m = 1 : M
                alP(m,1) = min(1,epsP/gm(m)/snr);
                alS(m,1) = 1-alP(m,1);
                etam(m,1) = epsP*(gm(m)+1/snr)/gm(m)/(1+epsP);
            end
            tauD = 0;
            for j = Dsize+1:M
                tauD = tauD + hm(j)*alP(j);
            end
            tauD = tauD + 1/snr;

            A = []; % No other constraints
            b = [];
            Aeq = [];%-eye(2*Dsize);
            beq = [];%zeros(2*Dsize,1);
            lb = [];
            ub = [];
            x0 = zeros(2*Dsize,1);
            options = optimoptions('fmincon','Display', 'off','OptimalityTolerance',10^(-20), 'StepTolerance', 10^(-20), 'FunctionTolerance', 10^(-20));
            x = fmincon(@(x) -sum(sqrt(hm(1:Dsize)).*x(1:Dsize)),x0,A,b,Aeq,beq,lb,ub,@(x) mycons(x,hm,Dsize,tauD,epsP,etam),options);

            if max(mycons(x,hm,Dsize,tauD,epsP,etam))>=0
                ratesz(Dsize)=0;
            else

                alS(1:Dsize) = (x(1:Dsize)).^2;
                alP(1:Dsize) = x(Dsize+1:end);

                for m = 1 : Dsize
                    Rmptilde(m) = log2(1+ hm(m)*alP(m)/...
                        (sum(hm(m+1:end).*alP(m+1:end)) ...
                        + sum(sqrt(hm(1:Dsize).*alS(1:Dsize)))^2 + 1/snr));
                end

                ratesz(Dsize) = log2(1+sum(sqrt(hm(1:Dsize).*alS(1:Dsize)))^2/...
                    (sum(hm(Dsize+1:end).*alP(Dsize+1:end)) + 1/snr));
            end
        end
        rates(ict) = max(ratesz);
 
        %behchmark %%%%%%%
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
        
        temp4=[];
        %probability way
        for m = 1 : 1
            him = hm; him(m)=[];
            alphai = alphamPI; alphai(m) = [];
            temp1 = sum(him.*alphai);
            if alphaSII(m)>0
                temp4 = [temp4 hm(m)*alphaSII(m)/(temp1+1/snr)];
            end
        end 
        
        %rate
        if isempty(temp4) %no beam  can be used or no power can be assingd
            ratex1(ict) = 0;
        else
            ratex1(ict) = log2(1+max(temp4));
        end
         
    end  
    Ratex(i) = mean(ratex1);
    Rates(i) = mean(rates);
end

plot(snrdb,Ratex,'-o',snrdb,Rates);%,snrdb,Rate2)


function [c,ceq] = mycons(x,hm,Dsize,tauD,epsP,etam)
%first Dsize element of x is x, the second half is alpha^P
for i = 1: Dsize         
    c(3*(i-1)+1) =  epsP*(sum(sqrt(hm(1:Dsize)).*x(1:Dsize)))^2 +...
        epsP*sum(hm(i+1:Dsize).*x(i+1:Dsize)) - hm(i)*x(Dsize+i)+epsP*tauD;
    c(3*(i-1)+2) = -x(Dsize+i)+etam(i);
    c(3*(i-1)+3) = x(i)^2 + x(Dsize+i) -1;
end
    c(3*Dsize+1:3*Dsize+2*Dsize) = -x;
    ceq = [];
 
end
