%% Constantes
M=5
L=3
delta=(L-1)/2
err=1e-10

%Tamanho da sequencia
n_s=10;

%Valores possíveis na constelação qam4 
qam4=[sqrt(2)/2+sqrt(2)/2j  sqrt(2)/2-sqrt(2)/2j ...
     -sqrt(2)/2+sqrt(2)/2j -sqrt(2)/2-sqrt(2)/2j];
 
%Valores possíveis no BPSK
bpsk=[-1 1];

%Cria vetor de tamanho n_s com valores 
% pegos aleatoriamente com probabilidade uniforme do vetor qam4
Sqam4=randsample(qam4,n_s,true)
%Cria vetor de tamanho n_s com valores 
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

%% Encher zeros
Sqam4_ifft_pad=[zeros(delta,length(Sqam4)/M); Sqam4_ifft; zeros(delta,length(Sqam4)/M)]
Sbpsk_ifft_pad=[zeros(delta,length(Sbpsk)/M); Sbpsk_ifft; zeros(delta,length(Sbpsk)/M)]

%% (b) Canais

F=rand(1,L)

H_col=[F(1) zeros(1,M+2*delta-L)];
H_lin=[F zeros(1,M-L+2*delta)];
H=toeplitz(H_col,H_lin)

%Matriz de Transmissão de zeros
X=[zeros(delta,M);eye(M);zeros(delta,M)]

%Modelo do canal com transmissão de zeros
H0=H*X

%% Passar pelo canal
Yk_qam4=H*Sqam4_ifft_pad
Yk_bpsk=H*Sbpsk_ifft_pad

%%Inversão do canal

S_qam4_par_tilde_ifft=inv(H0)*Yk_qam4
S_bpsk_par_tilde_ifft=inv(H0)*Yk_bpsk
%% DFT

S_qam4_par_tilde=fft(S_qam4_par_tilde_ifft)
S_bpsk_par_tilde=fft(S_bpsk_par_tilde_ifft)


%%Paralelo para serial

S_qam4_tilde=reshape(S_qam4_par_tilde,[1 numel(S_qam4_par_tilde)])
S_bpsk_tilde=reshape(S_bpsk_par_tilde,[1 numel(S_bpsk_par_tilde)])

%% Recuperar num de elem

S_qam4_tilde=S_qam4_par_tilde(M:end-(M-mod((length(Sqam4_init)+M-1),M))*(0~=mod((length(Sqam4_init)+M-1),M)))
S_bpsk_tilde=S_bpsk_par_tilde(M:end-(M-mod((length(Sbpsk_init)+M-1),M))*(0~=mod((length(Sbpsk_init)+M-1),M)))


% Decisao

S_qam4_tilde=sqrt(2)/2*((real(S_qam4_tilde)>0&imag(S_qam4_tilde)>0)*(1+j)...
   +(real(S_qam4_tilde)>0&imag(S_qam4_tilde)<0)*(1-j)...
  +(real(S_qam4_tilde)<0&imag(S_qam4_tilde)>0)*(-1+j)...
  +(real(S_qam4_tilde)<0&imag(S_qam4_tilde)<0)*(-1-j))


S_bpsk_par_invtilde=(real(S_bpsk_tilde)>0)*(1)...
  +(real(S_bpsk_tilde)<0)*(-1)

%Calculo de erro

qam4_err=(Sqam4_init-S_qam4_tilde)<err;
qam4_ber=100-sum(qam4_err)*100/n_s

bpsk_err=(Sbpsk_init-S_bpsk_par_invtilde)<err;
bpsk_ber=100-sum(bpsk_err)*100/n_s
