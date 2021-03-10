close all
clear
%% Constantes
M=50;
L=41;
delta=(L-1);
err=1e-10;
SNR=-50:2:50;

%Tamanho da sequencia
n_s=1000;
N_canais=10000;

%Valores possíveis na constelação qam4 
qam4=[sqrt(2)/2+sqrt(2)/2j  sqrt(2)/2-sqrt(2)/2j ...
     -sqrt(2)/2+sqrt(2)/2j -sqrt(2)/2-sqrt(2)/2j];
 
%Valores possíveis no BPSK
bpsk=[-1 1];

%Cria vetor de tamanho n_s com valores 
% pegos aleatoriamente com probabilidade uniforme do vetor qam4
Sqam4=randsample(qam4,n_s,true);
%Cria vetor de tamanho n_s com valores 
% pegos aleatoriamente com probabilidade uniforme do vetor bpsk
Sbpsk=randsample(bpsk,n_s,true);


%Fazendo bkp para posterior consulta
Sqam4_init=Sqam4;
Sbpsk_init=Sbpsk;

%% Serial para paralelo
% Completa com zeros para ser múltiplo de M para não ter problema na matriz

Sqam4=[repmat(0,1,M-1) Sqam4_init...
       repmat(0,1,(M-mod((length(Sqam4_init)+M-1),M))*(0~=mod((length(Sqam4_init)+M-1),M)))];
 
Sbpsk=[repmat(0,1,M-1) Sbpsk_init ...
       repmat(0,1,(M-mod((length(Sbpsk_init)+M-1),M))*(0~=mod((length(Sbpsk_init)+M-1),M)))];
 

Sqam4_par=reshape(Sqam4,[M length(Sqam4)/M]);
Sbpsk_par=reshape(Sbpsk,[M length(Sbpsk)/M]);

%% IDFT

Sqam4_ifft=ifft(Sqam4_par);
Sbpsk_ifft=ifft(Sbpsk_par);

%% Encher zeros
Sqam4_ifft_pad=[zeros(delta,length(Sqam4)/M); Sqam4_ifft; zeros(delta,length(Sqam4)/M)];
Sbpsk_ifft_pad=[zeros(delta,length(Sbpsk)/M); Sbpsk_ifft; zeros(delta,length(Sbpsk)/M)];


wait=waitbar(0,'Calculando...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
for Snr_atual=1:length(SNR)
%% (b) Canais
for canal_atual=1:N_canais
    if getappdata(wait,'canceling')
        break
    end
F=rand(1,L);

H_col=[F(1) zeros(1,M+2*delta-L)];
H_lin=[F zeros(1,M-L+2*delta)];
H=toeplitz(H_col,H_lin);

%Matriz de Transmissão de zeros
X=[zeros(delta,M);eye(M);zeros(delta,M)];

%Modelo do canal com transmissão de zeros
H0=H*X;

%% Passar pelo canal
Yk_qam4=H*Sqam4_ifft_pad;
Yk_bpsk=H*Sbpsk_ifft_pad;

%%Ruído
Yk_qam4_noise=awgn(Yk_qam4,SNR(Snr_atual),'measured');
Yk_bpsk_noise=awgn(Yk_bpsk,SNR(Snr_atual),'measured');

X_cicl=horzcat([eye(M)],[eye(L-1) ;zeros(M-L+1,L-1)]);

%% Inclusão de "sufixo" cíclico
Yk_qam4_cicl=X_cicl*Yk_qam4_noise;
Yk_bpsk_cicl=X_cicl*Yk_bpsk_noise;

%%"Inversão" do canal
inv_canal=inv(X_cicl*H0);
S_qam4_par_tilde_ifft=inv_canal*Yk_qam4_cicl;
S_bpsk_par_tilde_ifft=inv_canal*Yk_bpsk_cicl;


%% DFT

S_qam4_par_tilde=fft(S_qam4_par_tilde_ifft);
S_bpsk_par_tilde=fft(S_bpsk_par_tilde_ifft);

%%Paralelo para serial

S_qam4_tilde=reshape(S_qam4_par_tilde,[1 numel(S_qam4_par_tilde)]);
S_bpsk_tilde=reshape(S_bpsk_par_tilde,[1 numel(S_bpsk_par_tilde)]);

%% Recuperar num de elem

S_qam4_tilde=S_qam4_par_tilde(M:end-(M-mod((length(Sqam4_init)+M-1),M))*(0~=mod((length(Sqam4_init)+M-1),M)));
S_bpsk_tilde=S_bpsk_par_tilde(M:end-(M-mod((length(Sbpsk_init)+M-1),M))*(0~=mod((length(Sbpsk_init)+M-1),M)));


% Decisao

S_qam4_tilde=sqrt(2)/2*((real(S_qam4_tilde)>0&imag(S_qam4_tilde)>0)*(1+j)...
   +(real(S_qam4_tilde)>0&imag(S_qam4_tilde)<0)*(1-j)...
  +(real(S_qam4_tilde)<0&imag(S_qam4_tilde)>0)*(-1+j)...
  +(real(S_qam4_tilde)<0&imag(S_qam4_tilde)<0)*(-1-j));


S_bpsk_par_invtilde=(real(S_bpsk_tilde)>0)*(1)...
  +(real(S_bpsk_tilde)<0)*(-1);

%Calculo de erro

qam4_err=(Sqam4_init-S_qam4_tilde)<err;
qam4_ser(canal_atual)=100-sum(qam4_err)*100/n_s;

bpsk_err=(Sbpsk_init-S_bpsk_par_invtilde)<err;
bpsk_ser(canal_atual)=100-sum(bpsk_err)*100/n_s;
end %Loop dos Canais

qam4_mean_ser(Snr_atual)=mean2(qam4_ser);
bpsk_mean_ser(Snr_atual)=mean2(bpsk_ser);
waitbar((Snr_atual)/(length(SNR)),wait,sprintf(['%1.2f' '%% Concluído' ],100*(Snr_atual)/(length(SNR))));

if getappdata(wait,'canceling')
        break
    end
end %Loop dos ruidos
delete(wait)

plot1=figure;
set(gcf, 'PaperPosition', [0 0 8 6])
set(gcf, 'PaperSize', [8 6])
plot(SNR(1:Snr_atual),qam4_mean_ser)
title('QAM4 Overlap-Add A')
xlabel('SNR')
ylabel('Probabilidade de Erro (%)')
saveas(plot1,['ovadd_invcicl-' 'QAM4-n_s-' num2str(n_s) '-N_canais_' num2str(N_canais) '.pdf'])

plot2=figure;
set(gcf, 'PaperPosition', [0 0 8 6])
set(gcf, 'PaperSize', [8 6])
plot(SNR(1:Snr_atual),bpsk_mean_ser)
title('BPSK Overlap-Add A')
xlabel('SNR')
ylabel('Probabilidade de Erro (%)')
saveas(plot2,['ovadd_invcicl-' 'BPSK-n_s-' num2str(n_s) '-N_canais_' num2str(N_canais) '.pdf'])
