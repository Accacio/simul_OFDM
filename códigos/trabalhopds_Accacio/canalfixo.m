close all

%% CANAL FIXO Overlap-Save Utilizando DFT e inversão diagonal F=[0.133957201146203 0.210366058971429 0.765448616816247 0.128977964632708 0.885880555211058 0.745044777555957 0.169805477542131 0.702511920757993 0.0453414709241675 0.400278773040264 0.348597710539630 0.566410419719651 0.773814795695134 0.808922861924320 0.548979215418482 0.998182621012049 0.453870397584661 0.417274496867465 0.598585668522936 0.0435876542682528 0.600003351340846 0.617996213550836 0.0376328814111715 0.123033078034647 0.0561570905855576 0.736035461066522 0.924548068225402 0.120688016672555 0.974744723554751 0.571772370857491 0.863063007350897 0.635420076254051 0.430152201072417 0.0505642926558801 0.953063335249642 0.540852461486257 0.794477743189057 0.101275932582940 0.181425941266032 0.457374357379481 0.441434273038271]

clear

% Constantes
M=50;
L=41;
delta=0;
err=1e-10;
cicl=L-1;
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

% Serial para paralelo
% Completa com zeros para ser múltiplo de M para não ter problema na matriz

Sqam4=[repmat(0,1,M-1) Sqam4_init...
       repmat(0,1,(M-mod((length(Sqam4_init)+M-1),M))*(0~=mod((length(Sqam4_init)+M-1),M)))];
 
Sbpsk=[repmat(0,1,M-1) Sbpsk_init ...
       repmat(0,1,(M-mod((length(Sbpsk_init)+M-1),M))*(0~=mod((length(Sbpsk_init)+M-1),M)))];
 

Sqam4_par=reshape(Sqam4,[M length(Sqam4)/M]);
Sbpsk_par=reshape(Sbpsk,[M length(Sbpsk)/M]);

% IDFT

Sqam4_ifft=ifft(Sqam4_par);
Sbpsk_ifft=ifft(Sbpsk_par);

%Coloca Prefixo

X_cicl=horzcat([eye(M)],[eye(cicl);zeros(M-(cicl),cicl)])';

Sqam4_ifft_cicl=X_cicl*Sqam4_ifft;

Sbpsk_ifft_cicl=X_cicl*Sbpsk_ifft;

wait=waitbar(0,'Calculando...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
for Snr_atual=1:length(SNR)

% (b) Canais
    
F=[0.133957201146203 0.210366058971429 0.765448616816247 0.128977964632708 0.885880555211058 0.745044777555957 0.169805477542131 0.702511920757993 0.0453414709241675 0.400278773040264 0.348597710539630 0.566410419719651 0.773814795695134 0.808922861924320 0.548979215418482 0.998182621012049 0.453870397584661 0.417274496867465 0.598585668522936 0.0435876542682528 0.600003351340846 0.617996213550836 0.0376328814111715 0.123033078034647 0.0561570905855576 0.736035461066522 0.924548068225402 0.120688016672555 0.974744723554751 0.571772370857491 0.863063007350897 0.635420076254051 0.430152201072417 0.0505642926558801 0.953063335249642 0.540852461486257 0.794477743189057 0.101275932582940 0.181425941266032 0.457374357379481 0.441434273038271];

H_col=[F(1) zeros(1,M+cicl+2*delta-L)];
H_lin=[F zeros(1,M+cicl-L+2*delta)];
H=toeplitz(H_col,H_lin);

%Matriz de Transmissão de zeros
X=[zeros(delta,M);eye(M);zeros(delta,M)];


%Modelo do canal com transmissão de prefixo cíclico
H0=H*X_cicl;
C=H0;
Lambda=ifft(fft(C).').';

% Passar pelo canal
Yk_qam4=H*Sqam4_ifft_cicl;
Yk_bpsk=H*Sbpsk_ifft_cicl;

% Ruído

Yk_qam4_noise=awgn(Yk_qam4,SNR(Snr_atual),'measured');
Yk_bpsk_noise=awgn(Yk_bpsk,SNR(Snr_atual),'measured');

% DFT

S_qam4_par_tilde_ifft=fft(Yk_qam4_noise);
S_bpsk_par_tilde_ifft=fft(Yk_bpsk_noise);

%"Inversão" do canal
inv_canal=inv(Lambda);
S_qam4_par_tilde=inv_canal*S_qam4_par_tilde_ifft;
S_bpsk_par_tilde=inv_canal*S_bpsk_par_tilde_ifft;

%Paralelo para serial

S_qam4_tilde=reshape(S_qam4_par_tilde,[1 numel(S_qam4_par_tilde)]);
S_bpsk_tilde=reshape(S_bpsk_par_tilde,[1 numel(S_bpsk_par_tilde)]);

% Recuperar num de elem

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
qam4_ser=100-sum(qam4_err)*100/n_s;

bpsk_err=(Sbpsk_init-S_bpsk_par_invtilde)<err;
bpsk_ser=100-sum(bpsk_err)*100/n_s;

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
title('QAM4 Overlap-Save')
xlabel('SNR')
ylabel('Razão de Erro de Símbolo')
saveas(plot1,['canalfixo_ovsav-' 'QAM4-n_s-' num2str(n_s) '.pdf'])


plot2=figure;
set(gcf, 'PaperPosition', [0 0 8 6])
set(gcf, 'PaperSize', [8 6])
plot(SNR(1:Snr_atual),bpsk_mean_ser)
title('BPSK Overlap-Save')
xlabel('SNR')
ylabel('Razão de Erro de Símbolo')
saveas(plot2,['canalfixo_ovsav-' 'BPSK-n_s-' num2str(n_s) '.pdf'])









%% CANAL FIXO Overlap-Add Invertendo Matriz círculante F=[0.133957201146203 0.210366058971429 0.765448616816247 0.128977964632708 0.885880555211058 0.745044777555957 0.169805477542131 0.702511920757993 0.0453414709241675 0.400278773040264 0.348597710539630 0.566410419719651 0.773814795695134 0.808922861924320 0.548979215418482 0.998182621012049 0.453870397584661 0.417274496867465 0.598585668522936 0.0435876542682528 0.600003351340846 0.617996213550836 0.0376328814111715 0.123033078034647 0.0561570905855576 0.736035461066522 0.924548068225402 0.120688016672555 0.974744723554751 0.571772370857491 0.863063007350897 0.635420076254051 0.430152201072417 0.0505642926558801 0.953063335249642 0.540852461486257 0.794477743189057 0.101275932582940 0.181425941266032 0.457374357379481 0.441434273038271]

clear
% Constantes
M=50;
L=41;
delta=(L-1);
err=1e-10;
SNR=-50:2:50;

%Tamanho da sequencia
n_s=1000;

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

% Serial para paralelo
% Completa com zeros para ser múltiplo de M para não ter problema na matriz

Sqam4=[repmat(0,1,M-1) Sqam4_init...
       repmat(0,1,(M-mod((length(Sqam4_init)+M-1),M))*(0~=mod((length(Sqam4_init)+M-1),M)))];
 
Sbpsk=[repmat(0,1,M-1) Sbpsk_init ...
       repmat(0,1,(M-mod((length(Sbpsk_init)+M-1),M))*(0~=mod((length(Sbpsk_init)+M-1),M)))];
 

Sqam4_par=reshape(Sqam4,[M length(Sqam4)/M]);
Sbpsk_par=reshape(Sbpsk,[M length(Sbpsk)/M]);

% IDFT

Sqam4_ifft=ifft(Sqam4_par);
Sbpsk_ifft=ifft(Sbpsk_par);

% Encher zeros
Sqam4_ifft_pad=[zeros(delta,length(Sqam4)/M); Sqam4_ifft; zeros(delta,length(Sqam4)/M)];
Sbpsk_ifft_pad=[zeros(delta,length(Sbpsk)/M); Sbpsk_ifft; zeros(delta,length(Sbpsk)/M)];


wait=waitbar(0,'Calculando...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
for Snr_atual=1:length(SNR)
% (b) Canais


F=[0.133957201146203 0.210366058971429 0.765448616816247 0.128977964632708 0.885880555211058 0.745044777555957 0.169805477542131 0.702511920757993 0.0453414709241675 0.400278773040264 0.348597710539630 0.566410419719651 0.773814795695134 0.808922861924320 0.548979215418482 0.998182621012049 0.453870397584661 0.417274496867465 0.598585668522936 0.0435876542682528 0.600003351340846 0.617996213550836 0.0376328814111715 0.123033078034647 0.0561570905855576 0.736035461066522 0.924548068225402 0.120688016672555 0.974744723554751 0.571772370857491 0.863063007350897 0.635420076254051 0.430152201072417 0.0505642926558801 0.953063335249642 0.540852461486257 0.794477743189057 0.101275932582940 0.181425941266032 0.457374357379481 0.441434273038271];


H_col=[F(1) zeros(1,M+2*delta-L)];
H_lin=[F zeros(1,M-L+2*delta)];
H=toeplitz(H_col,H_lin);

%Matriz de Transmissão de zeros
X=[zeros(delta,M);eye(M);zeros(delta,M)];

%Modelo do canal com transmissão de zeros
H0=H*X;

% Passar pelo canal
Yk_qam4=H*Sqam4_ifft_pad;
Yk_bpsk=H*Sbpsk_ifft_pad;

%Ruído
Yk_qam4_noise=awgn(Yk_qam4,SNR(Snr_atual),'measured');
Yk_bpsk_noise=awgn(Yk_bpsk,SNR(Snr_atual),'measured');

X_cicl=horzcat([eye(M)],[eye(L-1) ;zeros(M-L+1,L-1)]);

% Inclusão de "sufixo" cíclico
Yk_qam4_cicl=X_cicl*Yk_qam4_noise;
Yk_bpsk_cicl=X_cicl*Yk_bpsk_noise;

%"Inversão" do canal
inv_canal=inv(X_cicl*H0);
S_qam4_par_tilde_ifft=inv_canal*Yk_qam4_cicl;
S_bpsk_par_tilde_ifft=inv_canal*Yk_bpsk_cicl;


% DFT

S_qam4_par_tilde=fft(S_qam4_par_tilde_ifft);
S_bpsk_par_tilde=fft(S_bpsk_par_tilde_ifft);

%Paralelo para serial

S_qam4_tilde=reshape(S_qam4_par_tilde,[1 numel(S_qam4_par_tilde)]);
S_bpsk_tilde=reshape(S_bpsk_par_tilde,[1 numel(S_bpsk_par_tilde)]);

% Recuperar num de elem

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
qam4_ser=100-sum(qam4_err)*100/n_s;

bpsk_err=(Sbpsk_init-S_bpsk_par_invtilde)<err;
bpsk_ser=100-sum(bpsk_err)*100/n_s;


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
ylabel('Razão de Erro de Símbolo')
saveas(plot1,['canalfixo_ovadd_invcicl-' 'QAM4-n_s-' num2str(n_s) '.pdf'])

plot2=figure;
set(gcf, 'PaperPosition', [0 0 8 6])
set(gcf, 'PaperSize', [8 6])
plot(SNR(1:Snr_atual),bpsk_mean_ser)
title('BPSK Overlap-Add A')
xlabel('SNR')
ylabel('Razão de Erro de Símbolo')
saveas(plot2,['canalfixo_ovadd_invcicl-' 'BPSK-n_s-' num2str(n_s) '.pdf'])

%% CANAL FIXO Overlap-Add Utilizando DFT e inversão diagonal F=[0.133957201146203 0.210366058971429 0.765448616816247 0.128977964632708 0.885880555211058 0.745044777555957 0.169805477542131 0.702511920757993 0.0453414709241675 0.400278773040264 0.348597710539630 0.566410419719651 0.773814795695134 0.808922861924320 0.548979215418482 0.998182621012049 0.453870397584661 0.417274496867465 0.598585668522936 0.0435876542682528 0.600003351340846 0.617996213550836 0.0376328814111715 0.123033078034647 0.0561570905855576 0.736035461066522 0.924548068225402 0.120688016672555 0.974744723554751 0.571772370857491 0.863063007350897 0.635420076254051 0.430152201072417 0.0505642926558801 0.953063335249642 0.540852461486257 0.794477743189057 0.101275932582940 0.181425941266032 0.457374357379481 0.441434273038271]

clear
% Constantes
M=50;
L=41;
delta=(L-1);
err=1e-10;
SNR=-50:2:50;

%Tamanho da sequencia
n_s=1000;

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

% Serial para paralelo
% Completa com zeros para ser múltiplo de M para não ter problema na matriz

Sqam4=[repmat(0,1,M-1) Sqam4_init...
       repmat(0,1,(M-mod((length(Sqam4_init)+M-1),M))*(0~=mod((length(Sqam4_init)+M-1),M)))];
 
Sbpsk=[repmat(0,1,M-1) Sbpsk_init ...
       repmat(0,1,(M-mod((length(Sbpsk_init)+M-1),M))*(0~=mod((length(Sbpsk_init)+M-1),M)))];
 

Sqam4_par=reshape(Sqam4,[M length(Sqam4)/M]);
Sbpsk_par=reshape(Sbpsk,[M length(Sbpsk)/M]);

% IDFT

Sqam4_ifft=ifft(Sqam4_par);
Sbpsk_ifft=ifft(Sbpsk_par);

% Encher zeros
Sqam4_ifft_pad=[zeros(delta,length(Sqam4)/M); Sqam4_ifft; zeros(delta,length(Sqam4)/M)];
Sbpsk_ifft_pad=[zeros(delta,length(Sbpsk)/M); Sbpsk_ifft; zeros(delta,length(Sbpsk)/M)];


wait=waitbar(0,'Calculando...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
for Snr_atual=1:length(SNR)
% (b) Canais

    
F=[0.133957201146203 0.210366058971429 0.765448616816247 0.128977964632708 0.885880555211058 0.745044777555957 0.169805477542131 0.702511920757993 0.0453414709241675 0.400278773040264 0.348597710539630 0.566410419719651 0.773814795695134 0.808922861924320 0.548979215418482 0.998182621012049 0.453870397584661 0.417274496867465 0.598585668522936 0.0435876542682528 0.600003351340846 0.617996213550836 0.0376328814111715 0.123033078034647 0.0561570905855576 0.736035461066522 0.924548068225402 0.120688016672555 0.974744723554751 0.571772370857491 0.863063007350897 0.635420076254051 0.430152201072417 0.0505642926558801 0.953063335249642 0.540852461486257 0.794477743189057 0.101275932582940 0.181425941266032 0.457374357379481 0.441434273038271];

H_col=[F(1) zeros(1,M+2*delta-L)];
H_lin=[F zeros(1,M-L+2*delta)];
H=toeplitz(H_col,H_lin);

%Matriz de Transmissão de zeros
X=[zeros(delta,M);eye(M);zeros(delta,M)];

%Modelo do canal com transmissão de zeros
H0=H*X;



% Passar pelo canal
Yk_qam4=H*Sqam4_ifft_pad;
Yk_bpsk=H*Sbpsk_ifft_pad;

%Ruído
Yk_qam4_noise=awgn(Yk_qam4,SNR(Snr_atual),'measured');
Yk_bpsk_noise=awgn(Yk_bpsk,SNR(Snr_atual),'measured');

X_cicl=horzcat([eye(M)],[eye(L-1) ;zeros(M-L+1,L-1)]);

% Inclusão de "sufixo" cíclico
Yk_qam4_cicl=X_cicl*Yk_qam4_noise;
Yk_bpsk_cicl=X_cicl*Yk_bpsk_noise;

C=X_cicl*H0;
Lambda=ifft(fft(C).').';


% DFT

S_qam4_par_tilde=fft(Yk_qam4_cicl);
S_bpsk_par_tilde=fft(Yk_bpsk_cicl);

%"Inversão" do canal

inv_canal=inv(Lambda);
S_qam4_par_tilde=inv_canal*S_qam4_par_tilde;
S_bpsk_par_tilde=inv_canal*S_bpsk_par_tilde;


%Paralelo para serial

S_qam4_tilde=reshape(S_qam4_par_tilde,[1 numel(S_qam4_par_tilde)]);
S_bpsk_tilde=reshape(S_bpsk_par_tilde,[1 numel(S_bpsk_par_tilde)]);

% Recuperar num de elem

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
qam4_ser=100-sum(qam4_err)*100/n_s;

bpsk_err=(Sbpsk_init-S_bpsk_par_invtilde)<err;
bpsk_ser=100-sum(bpsk_err)*100/n_s;

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
title('QAM4 Overlap-Add B')
xlabel('SNR')
ylabel('Razão de Erro de Símbolo')
saveas(plot1,['canalfixo_ovadd_ofdm-' 'QAM4-n_s-' num2str(n_s) '.pdf'])

plot2=figure;
set(gcf, 'PaperPosition', [0 0 8 6])
set(gcf, 'PaperSize', [8 6])
plot(SNR(1:Snr_atual),bpsk_mean_ser)
title('BPSK Overlap-Add B')
xlabel('SNR')
ylabel('Razão de Erro de Símbolo')
saveas(plot2,['canalfixo_ovadd_ofdm-' 'BPSK-n_s-' num2str(n_s) '.pdf'])

%% CANAL FIXO Redundancia Mínima F=[0.133957201146203 0.210366058971429 0.765448616816247 0.128977964632708 0.885880555211058 0.745044777555957 0.169805477542131 0.702511920757993 0.0453414709241675 0.400278773040264 0.348597710539630 0.566410419719651 0.773814795695134 0.808922861924320 0.548979215418482 0.998182621012049 0.453870397584661 0.417274496867465 0.598585668522936 0.0435876542682528 0.600003351340846 0.617996213550836 0.0376328814111715 0.123033078034647 0.0561570905855576 0.736035461066522 0.924548068225402 0.120688016672555 0.974744723554751 0.571772370857491 0.863063007350897 0.635420076254051 0.430152201072417 0.0505642926558801 0.953063335249642 0.540852461486257 0.794477743189057 0.101275932582940 0.181425941266032 0.457374357379481 0.441434273038271]

clear
% Constantes
M=50;
L=41;
delta=(L-1)/2;
err=1e-10;
SNR=-50:2:50;

%Tamanho da sequencia
n_s=1000;

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

% Serial para paralelo
% Completa com zeros para ser múltiplo de M para não ter problema na matriz

Sqam4=[repmat(0,1,M-1) Sqam4_init...
       repmat(0,1,(M-mod((length(Sqam4_init)+M-1),M))*(0~=mod((length(Sqam4_init)+M-1),M)))];
 
Sbpsk=[repmat(0,1,M-1) Sbpsk_init ...
       repmat(0,1,(M-mod((length(Sbpsk_init)+M-1),M))*(0~=mod((length(Sbpsk_init)+M-1),M)))];
 

Sqam4_par=reshape(Sqam4,[M length(Sqam4)/M]);
Sbpsk_par=reshape(Sbpsk,[M length(Sbpsk)/M]);

% IDFT

Sqam4_ifft=ifft(Sqam4_par);
Sbpsk_ifft=ifft(Sbpsk_par);

% Encher zeros
Sqam4_ifft_pad=[zeros(delta,length(Sqam4)/M); Sqam4_ifft; zeros(delta,length(Sqam4)/M)];
Sbpsk_ifft_pad=[zeros(delta,length(Sbpsk)/M); Sbpsk_ifft; zeros(delta,length(Sbpsk)/M)];



wait=waitbar(0,'Calculando...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
for Snr_atual=1:length(SNR)

% (b) Canal


    
F=[0.133957201146203 0.210366058971429 0.765448616816247 0.128977964632708 0.885880555211058 0.745044777555957 0.169805477542131 0.702511920757993 0.0453414709241675 0.400278773040264 0.348597710539630 0.566410419719651 0.773814795695134 0.808922861924320 0.548979215418482 0.998182621012049 0.453870397584661 0.417274496867465 0.598585668522936 0.0435876542682528 0.600003351340846 0.617996213550836 0.0376328814111715 0.123033078034647 0.0561570905855576 0.736035461066522 0.924548068225402 0.120688016672555 0.974744723554751 0.571772370857491 0.863063007350897 0.635420076254051 0.430152201072417 0.0505642926558801 0.953063335249642 0.540852461486257 0.794477743189057 0.101275932582940 0.181425941266032 0.457374357379481 0.441434273038271];

H_col=[F(1) zeros(1,M+2*delta-L)];
H_lin=[F zeros(1,M-L+2*delta)];
H=toeplitz(H_col,H_lin);

%Matriz de Transmissão de zeros
X=[zeros(delta,M);eye(M);zeros(delta,M)];

%Modelo do canal com transmissão de zeros
H0=H*X;

% Passar pelo canal
Yk_qam4=H*Sqam4_ifft_pad;
Yk_bpsk=H*Sbpsk_ifft_pad;

%Ruído
Yk_qam4_noise=awgn(Yk_qam4,SNR(Snr_atual),'measured');
Yk_bpsk_noise=awgn(Yk_bpsk,SNR(Snr_atual),'measured');

%"Inversão" do canal

inv_canal=inv(H0);

S_qam4_par_tilde_ifft=inv_canal*Yk_qam4_noise;
S_bpsk_par_tilde_ifft=inv_canal*Yk_bpsk_noise;
% DFT

S_qam4_par_tilde=fft(S_qam4_par_tilde_ifft);
S_bpsk_par_tilde=fft(S_bpsk_par_tilde_ifft);


%Paralelo para serial

S_qam4_tilde=reshape(S_qam4_par_tilde,[1 numel(S_qam4_par_tilde)]);
S_bpsk_tilde=reshape(S_bpsk_par_tilde,[1 numel(S_bpsk_par_tilde)]);

% Recuperar num de elem

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
qam4_ser=100-sum(qam4_err)*100/n_s;

bpsk_err=(Sbpsk_init-S_bpsk_par_invtilde)<err;
bpsk_ser=100-sum(bpsk_err)*100/n_s;

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
title('QAM4 Redundancia Mínima')
xlabel('SNR')
ylabel('Probabilidade de Erro (%)')
saveas(plot1,['canalfixo_redunmin-' 'QAM4-n_s-' num2str(n_s) '.pdf'])

plot2=figure;
set(gcf, 'PaperPosition', [0 0 8 6])
set(gcf, 'PaperSize', [8 6])
plot(SNR(1:Snr_atual),bpsk_mean_ser)
title('BPSK Redundancia Mínima')
xlabel('SNR')
ylabel('Probabilidade de Erro (%)')
saveas(plot2,['canalfixo_redunmin-' 'BPSK-n_s-' num2str(n_s) '.pdf'])