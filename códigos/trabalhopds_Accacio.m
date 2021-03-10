%%% Parte 1 - OFDM


%% Constantes

M=80
L=40

%Tamanho da sequencia
n_s=1000;   
                                                             
%Valores possíveis na constelação qam4 
qam4=[sqrt(2)/2+sqrt(2)/2j  sqrt(2)/2-sqrt(2)/2j ...
     -sqrt(2)/2+sqrt(2)/2j -sqrt(2)/2-sqrt(2)/2j];
 
%Valores possíveis no BPSK
bpsk=[-1 1];

%Cria vetor de tamanho n com valores 
% pegos aleatoriamente com probabilidade uniforme do vetor qam4
Sbpsk=randsample(bpsk,n_s,true);
%Cria vetor de tamanho n com valores 
% pegos aleatoriamente com probabilidade uniforme do vetor bpsk
Sqam4=randsample(qam4,n_s,true);

Sqam4_init=Sqam4
Sbpsk_init=Sbpsk
%% Serial para paralelo

% Completa com zeros para ser múltiplo de M
Sbpsk=[repmat(0,1,M-1) Sbpsk_init ...
    repmat(0,1,(M-mod((length(Sbpsk_init)+M-1),M))*(0~=mod((length(Sbpsk_init)+M-1),M)))];
  
Sqam4=[repmat(0,1,M-1) Sqam4_init ...
    repmat(0,1,(M-mod((length(Sqam4_init)+M-1),M))*(0~=mod((length(Sqam4_init)+M-1),M)))];


Sbpsk_par=reshape(Sbpsk,[M length(Sbpsk)/M]);
Sqam4_par=reshape(Sqam4,[M length(Sqam4)/M]);

%% IDFT
Sbar_bpsk=ifft(Sbpsk_par);
Sbar_qam4=ifft(Sqam4_par)

%% Prefixo Cíclico

X_bpsk_par=[repmat(Sbar_bpsk,[1 1]);repmat(Sbar_bpsk(1:L,:),[1 1])];
X_qam4_par=[repmat(Sbar_qam4,[1 1]);repmat(Sbar_qam4(1:L,:),[1 1])]

%% Paralelo para Serial
X_bpsk=reshape(X_bpsk_par,[1 numel(X_bpsk_par)]);
X_qam4=reshape(X_qam4_par,[1 numel(X_qam4_par)]);

%% (b) Canais

n_canais=1000;
L_canais=41;

F=rand(1,L_canais)

%Circular tipo Add Save
%Ajustar com N amostras de X
new_F=[F zeros(1,M-length(F))]
%
C=gallery('circul',new_F');

%% Passar pelo canal
Y_bpsk=filter(F,1,X_bpsk);
Y_qam4=filter(F,1,X_qam4);

%% Serial para paralelo

% Completa com zeros para ser múltiplo de M+L
new_M=M+L;

%Y_bpsk=[repmat(0,1,new_M-1) Y_bpsk ...
%     repmat(0,1,(new_M-mod((length(Y_bpsk)+new_M-1),new_M))*(0~=mod((length(Y_bpsk)+new_M-1),new_M)))];

%Y_qam4=[repmat(0,1,new_M-1) Y_qam4 ...
%     repmat(0,1,(new_M-mod((length(Y_qam4)+new_M-1),new_M))*(0~=mod((length(Y_qam4)+new_M-1),new_M)))];

Y_bpsk_par=reshape(Y_bpsk,[new_M length(Y_bpsk)/new_M]);
Y_qam4_par=reshape(Y_qam4,[new_M length(Y_qam4)/new_M])

%% Remover prefixo cíclico
Ylin_bpsk_par=Y_bpsk_par(1:end-L,:);
Ylin_qam4_par=Y_qam4_par(1:end-L,:)

%% DFT
Y_bpsk_fft=fft(Ylin_bpsk_par);
Y_qam4_fft=fft(Ylin_qam4_par)
%% (c) Equalizacao

%Diagonal 
Lambda=ifft(fft(C).').';
S_bpsk_par_tilde=inv(Lambda)*Y_bpsk_fft;
S_qam4_par_tilde=inv(Lambda)*Y_qam4_fft;

S_bpsk_par_tilde = S_bpsk_par_tilde./sqrt(abs(S_bpsk_par_tilde).^2);
S_qam4_par_tilde = S_qam4_par_tilde./sqrt(abs(S_qam4_par_tilde).^2);
Sqam4_par

% Decision
S_qam4_par_invtilde=sqrt(2)/2*((real(S_qam4_par_tilde)>0&imag(S_qam4_par_tilde)>0)*(1+j)...
   +(real(S_qam4_par_tilde)>0&imag(S_qam4_par_tilde)<0)*(1-j)...
  +(real(S_qam4_par_tilde)<0&imag(S_qam4_par_tilde)>0)*(-1+j)...
  +(real(S_qam4_par_tilde)<0&imag(S_qam4_par_tilde)<0)*(-1-j))

S_bpsk_par_invtilde=(real(S_bpsk_par_tilde)>0)*(1)...
  +(real(S_bpsk_par_tilde)<0)*(-1)
%% Paralelo para Serial
S_qam4_invtilde=reshape(S_qam4_par_invtilde,[1 numel(S_qam4_par_invtilde)])
S_bpsk_invtilde=reshape(S_bpsk_par_invtilde,[1 numel(S_bpsk_par_invtilde)])
S_qam4_invtilde=S_qam4_invtilde(M+1:end+1-(M-mod((length(Sqam4_init)+M-1),M))*(0~=mod((length(Sqam4_init)+M-1),M)))
S_bpsk_invtilde=S_bpsk_invtilde(M+1:end+1-(M-mod((length(Sbpsk_init)+M-1),M))*(0~=mod((length(Sbpsk_init)+M-1),M)))
% S_qam4_par_invtilde = S_qam4_par_invtilde(S_qam4_par_invtilde~=0);
qam4_err=Sqam4_init==S_qam4_invtilde
bpsk_err=Sbpsk_init==S_bpsk_invtilde
qam4_ber=100-sum(qam4_err)*100/n_s
bpsk_ber=100-sum(bpsk_err)*100/n_s

% parte 2