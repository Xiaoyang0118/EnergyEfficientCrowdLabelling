clear all;
clc;

K = 20;%number of annotators
M = 1:1:10;%number of objects
N = 3;%type of clusters
B = 10^4;%bandwith (KHz)
D = [6*10^6,4*10^6,2*10^6];%data size (MB)
dK = [1,3,5];%annotator consumption
sigma = 10^(-5);%noise power (uW)
T = 1000;%transmission duration (s)
gamma = exp(D./(B*T))-1;%SNR threshold
dP = [sigma*T*(gamma(1)+gamma(2)*(1+gamma(1))+gamma(3)*(1+gamma(2))*(1+gamma(1))),sigma*T*(gamma(1)+gamma(2)*(1+gamma(1))),sigma*T*gamma(1)];%power consumption
P_BB = zeros(1,length(M));%initial output of the proposed design
P_1 = zeros(1,length(M));%initial output of the whole-layers design
P_2 = zeros(1,length(M));%initial output of the random clustering design
P_3 = zeros(1,length(M));%initial output of the first-layer design
h = ones(1,K);%channel gain
g = sort(h,'descend');%sorting channels with the descending gains
sumP_1 = 0;%initial power consumption of the whole-layers design
sumP_3 = 0;%initial power consumption of the first-layer design
    
for i = 1:1:length(M)
    k = 0;
    p = 0;
    kh = 0;
    ph = 0;
    Start = 1;
    End = 1;
    Starth = 1;
    Endh = 1; 
    knew = K+1;
    knewh = K+1;
    
    for m = 1:1:M(i)
        Start = End;
        End = {};
        kvect = {};
        pvect = {};
        for s = 1:1:length(Start)
            for n = 1:1:N
                knew = k(s) + dK(n);
                if knew <= K
                    pnew = p(s) + dP(n)/g(knew);
                    End = [End,length(End)+1];
                    kvect = [kvect,knew];
                    pvect = [pvect,pnew];
                end
            end
        end
        k = cell2mat(kvect);
        p = cell2mat(pvect);
    end    
    P_BB(i) = min(p);

    for m = 1:1:M(i)
        Starth = Endh;
        Endh = {};
        kvecth = {};
        pvecth = {};
        for s = 1:1:length(Starth)
            for n = 1:1:N
                knewh = kh(s) + dK(n);
                if knewh <= K
                    pnewh = ph(s) + dP(n)/h(knewh);
                    Endh = [Endh,length(Endh)+1];
                    kvecth = [kvecth,knewh];
                    pvecth = [pvecth,pnewh];
                end
            end
        end
        kh = cell2mat(kvecth);
        ph = cell2mat(pvecth);
    end    
    P_2(i) = min(ph);
    
    sumK_1 = dK(1)*M(i);
    if sumK_1 > K
        sumP_1 = 0;
    else
        sumP_1 = sumP_1 + dP(1)/g(sumK_1);
    end
    P_1(i) = sumP_1;
    
    sumK_3 = dK(3)*M(i);
    if sumK_3 > K
        sumP_3 = 0;
    else
        sumP_3 = sumP_3 + dP(3)/g(sumK_3);
    end
    P_3(i) = sumP_3;
end

semilogy(M,P_BB,'r',M,P_1,'b',M,P_2,'g',M,P_3,'k');
legend('Our proposed design', 'Whole-layers design', 'Random clustering design', 'First-layer design');
xlabel('Number of Objects');
ylabel('Power Consumption')