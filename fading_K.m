clear all;
clc;

K = 0:1:50;%number of annotators
P = 0.1;%power budget(W)
N = 3;%type of clusters
B = 10^4;%bandwith (KHz)
D = [6*10^6,4*10^6,2*10^6];%data size (MB)
dK = [1,3,5];%annotator consumption
sigma = 10^(-5);%noise power (uW)
T = 1000;%transmission duration (s)
gamma = exp(D./(B*T))-1;%SNR threshold
dP = [sigma*T*(gamma(1)+gamma(2)*(1+gamma(1))+gamma(3)*(1+gamma(2))*(1+gamma(1))),sigma*T*(gamma(1)+gamma(2)*(1+gamma(1))),sigma*T*gamma(1)];%power consumption
M_BB = zeros(1,length(K));%initial output of the proposed design
M_1 = zeros(1,length(K));%initial output of the whole-layers design
M_2 = zeros(1,length(K));%initial output of the random clustering design
M_3 = zeros(1,length(K));%initial output of the first-layer design
h = rand(1,length(K))%channel gain
g = sort(h,'descend');%sorting channels with the descending gains
% load('g_L.mat');

for i = 1:1:length(K)
    k = 0;
    p = 0;
    kh = 0;
    ph = 0;
    m = 0;
    Start = 1;
    End = 1;
    Starth = 1;
    Endh = 1;    
    m_1 = 0;
    m_2 = 0;
    m_3 = 0;
    sumP_1 = 0;
    sumP_3 = 0;
    sumK_1 = 0;
    sumK_3 = 0;
    knew = max(K)+1;
    pnew = P+1;
    knewh = max(K)+1;
    pnewh = P+1;

    while ~isempty(End)
        m = m + 1;
        Start = End;
        End = {};
        kvect = {};
        pvect = {};
        for s = 1:1:length(Start)
            for n = 1:1:N
                knew = k(s) + dK(n);
                if knew <= K(i) 
                    pnew = p(s) + dP(n)/g(knew);
                    if pnew <= P
                        End = [End,length(End)+1];
                        kvect = [kvect,knew];
                        pvect = [pvect,pnew];
                    end
                end
            end
        end
        k = cell2mat(kvect);
        p = cell2mat(pvect);
    end
    M_BB(i) = m - 1
    
    while ~isempty(Endh)
        m_2 = m_2 + 1;
        Starth = Endh;
        Endh = {};
        kvecth = {};
        pvecth = {};
        for s = 1:1:length(Starth)
            for n = 1:1:N
                knewh = kh(s) + dK(n);
                if knewh <= K(i) 
                    pnewh = ph(s) + dP(n)/h(knewh);
                    if pnewh <= P
                        Endh = [Endh,length(Endh)+1];
                        kvecth = [kvecth,knewh];
                        pvecth = [pvecth,pnewh];
                    end
                end
            end
        end
        kh = cell2mat(kvecth);
        ph = cell2mat(pvecth);
    end
    M_2(i) = m_2 - 1

    while sumP_1 <= P && sumK_1 <= K(i)
        m_1 = m_1 + 1;
        sumK_1 = dK(1)*m_1;
        if sumK_1 <= K(i)
            sumP_1 = sumP_1 + dP(1)/g(sumK_1);
        end
    end
    M_1(i) = m_1-1
    
    while sumP_3 <= P && sumK_3 <= K(i)
        m_3 = m_3 + 1;
        sumK_3 = dK(3)*m_3;
        if sumK_3 <= K(i)
            sumP_3 = sumP_3 + dP(3)/g(sumK_3);
        end
    end
    M_3(i) = m_3-1
end
plot(K,M_BB,'r',K,M_1,'b',K,M_2,'g',K,M_3,'k');
legend('Our proposed design', 'Whole-layers design', 'Random clustering design', 'First-layer design');
xlabel('Number of Annotators');
ylabel('Throughput (Number of Objects)')