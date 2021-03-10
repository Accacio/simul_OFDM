% parte 1
%% (a) Transmissao

% Criando o vetor s(n) dos tipos QAM4 e BPSK
n_s=100;

qam4=[sqrt(2)/2+sqrt(2)/2j  sqrt(2)/2-sqrt(2)/2j ... 
    -sqrt(2)/2+sqrt(2)/2j -sqrt(2)/2-sqrt(2)/2j];

bpsk=[-1 1]; 

% Cria vetor de tamanho n com valores aleatorios com probabilidade uniforme do vetor qam4
Sbpsk=randsample(bpsk,n_s,true)';

% Cria vetor de tamanho n com valores aleatorios com probabilidade uniforme do vetor bpsk
Sqam4=randsample(qam4,n_s,true)';                                           

%% (b) Canais


% Criando os canais aleatórios
n=1000;  % n = 10000
L=41;

H = zeros(n, L);
for k = 1:n
    H(k,:) = fir1(L-1,rand);
end

% Fourier inversa dos dois s(n)

Sbpsk2 = ifft(Sbpsk);
Sqam42 = ifft(Sqam4);

%% (c) Equalizacao



% parte 2