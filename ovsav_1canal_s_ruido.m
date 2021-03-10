%% Constantes
M=50
L=41
err=1e-10
%Tamanho da sequencia
n_s=1000;

%Valores possíveis na constelação qam4 
qam4=[sqrt(2)/2+sqrt(2)/2j  sqrt(2)/2-sqrt(2)/2j ...
     -sqrt(2)/2+sqrt(2)/2j -sqrt(2)/2-sqrt(2)/2j];
 
%Valores possíveis no BPSK
bpsk=[-1 1];
 

%Cria vetor de tamanho n com valores 
% pegos aleatoriamente com probabilidade uniforme do vetor qam4
Sqam4=randsample(qam4,n_s,true)


%Cria vetor de tamanho n com valores 
% pegos aleatoriamente com probabilidade uniforme do vetor bpsk
Sbpsk=randsample(bpsk,n_s,true)

%Fazendo bkp para posterior consulta
Sqam4_init=Sqam4;
Sbpsk_init=Sbpsk;

%% Serial para paralelo
% Completa com zeros para ser múltiplo de M para não ter problema na matriz

Sqam4=[repmat(0,1,M-1) Sqam4_init...
       repmat(0,1,(M-mod((length(Sqam4_init)+M-1),M))*(0~=mod((length(Sqam4_init)+M-1),M)))]
 
Sbpsk=[repmat(0,1,M-1) Sbpsk_init ...
       repmat(0,1,(M-mod((length(Sbpsk_init)+M-1),M))*(0~=mod((length(Sbpsk_init)+M-1),M)))];
 

Sqam4_par=reshape(Sqam4,[M length(Sqam4)/M])
Sbpsk_par=reshape(Sbpsk,[M length(Sbpsk)/M]);
%% IDFT

Sqam4_ifft=ifft(Sqam4_par)
Sbpsk_ifft=ifft(Sbpsk_par);

%% Prefixo Cíclico

Sqam4_ifft_cicl=[repmat(Sqam4_ifft(end-L+2:end,:),[1 1]);repmat(Sqam4_ifft,[1 1])]


Sbpsk_ifft_cicl=[repmat(Sbpsk_ifft(end-L+2:end,:),[1 1]);repmat(Sbpsk_ifft,[1 1])]

%% Paralelo para Serial
X_qam4=reshape(Sqam4_ifft_cicl,[1 numel(Sqam4_ifft_cicl)]);
X_bpsk=reshape(Sbpsk_ifft_cicl,[1 numel(Sbpsk_ifft_cicl)]);

%% (b) Canais

F=rand(1,L)

new_F=[F zeros(1,M-length(F))]
C=gallery('circul',new_F')


%% Passar pelo canal


Yk_qam4=filter(F,1,X_qam4)
Yk_bpsk=filter(F,1,X_bpsk);

%% Serial para paralelo

Y_qam4_cicl=reshape(Yk_qam4,[M+L-1 length(Yk_qam4)/(M+L-1)])
Y_bpsk_cicl=reshape(Yk_bpsk,[M+L-1 length(Yk_bpsk)/(M+L-1)])

%% Remover prefixo cíclico

Y_qam4=Y_qam4_cicl(end-M+1:end,:)
Y_bpsk=Y_bpsk_cicl(end-M+1:end,:)

%% DFT

Sbar_qam4=fft(Y_qam4).'
Sbar_bpsk=fft(Y_bpsk).'

%% (c) Equalizacao
%Diagonal 

Lambda=ifft(fft(C).').'

S_qam4_par_tilde=(inv(Lambda)*Sbar_qam4')
S_bpsk_par_tilde=(inv(Lambda)*Sbar_bpsk')

% Recuperar num de elem

S_qam4_tilde=reshape(S_qam4_par_tilde,[1 numel(S_qam4_par_tilde)]);
S_qam4_tilde=conj(S_qam4_par_tilde(M:end-(M-mod((length(Sqam4_init)+M-1),M))*(0~=mod((length(Sqam4_init)+M-1),M))))


S_bpsk_tilde=reshape(S_bpsk_par_tilde,[1 numel(S_bpsk_par_tilde)]);
S_bpsk_tilde=conj(S_bpsk_par_tilde(M:end-(M-mod((length(Sqam4_init)+M-1),M))*(0~=mod((length(Sqam4_init)+M-1),M))))

% Decisao

S_qam4_par_invtilde=sqrt(2)/2*((real(S_qam4_tilde)>0&imag(S_qam4_tilde)>0)*(1+j)...
   +(real(S_qam4_tilde)>0&imag(S_qam4_tilde)<0)*(1-j)...
  +(real(S_qam4_tilde)<0&imag(S_qam4_tilde)>0)*(-1+j)...
  +(real(S_qam4_tilde)<0&imag(S_qam4_tilde)<0)*(-1-j))


S_bpsk_par_invtilde=(real(S_bpsk_tilde)>0)*(1)...
  +(real(S_bpsk_tilde)<0)*(-1)

%Calculo de erro


qam4_err=(Sqam4_init-S_qam4_par_invtilde)<err;
qam4_ber=100-sum(qam4_err)*100/n_s

bpsk_err=(Sbpsk_init-S_bpsk_par_invtilde)<err;
bpsk_ber=100-sum(bpsk_err)*100/n_s



