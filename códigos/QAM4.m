% QAM-4 - Modula��o de Amplitude em Quadratura

M = 4;                                                                     % Tamanho da Modula��o
k = log2(M);                                                               % Bits necess�rios por valor
n = 200;                                                                   % N�mero de bits processados

vetorBinario = randi([0 1],n,1);                                           % Gerando vetor aleat�rio de n bits

vetorInput = reshape(vetorBinario,length(vetorBinario)/k,k);               % Empacotar os dados em k bits
x = bi2de(vetorInput);                                                     % Converter para inteiros

%Essas linhas s�o reservadas para voc� observar o formato dos valores
%gerados pelo programa e sua frequ�ncia de apari��o.
%-------------------------------------------------------------------------%
% figure;                                                 
% stem(x(1:n/2));
% title('Entrada Aleat�ria');
% xlabel('Index');
% ylabel('Valor Quadratura');
%-------------------------------------------------------------------------%

y = qammod(x,M,0);                                                         % Resultado modulado, mas com pot�ncia n�o unit�ria

potenciaMedia = mean(abs(y).^2);
S1n = y/sqrt(potenciaMedia);                                                % Resultado final, com pot�ncia unit�ria

%-------------------------------------------------------------------------%
%BPSK - Binary Phase Shift Keying

n = 100;                                                                   % N�mero de bits processados

vetorBinario = randi([0 1],n,1);                                           % Gerando vetor aleat�rio de n bits
S2n = 2*vetorBinario - 1;                                                   % Resultado final, com pot�ncia unit�ria