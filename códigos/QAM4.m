% QAM-4 - Modulação de Amplitude em Quadratura

M = 4;                                                                     % Tamanho da Modulação
k = log2(M);                                                               % Bits necessários por valor
n = 200;                                                                   % Número de bits processados

vetorBinario = randi([0 1],n,1);                                           % Gerando vetor aleatório de n bits

vetorInput = reshape(vetorBinario,length(vetorBinario)/k,k);               % Empacotar os dados em k bits
x = bi2de(vetorInput);                                                     % Converter para inteiros

%Essas linhas são reservadas para você observar o formato dos valores
%gerados pelo programa e sua frequência de aparição.
%-------------------------------------------------------------------------%
% figure;                                                 
% stem(x(1:n/2));
% title('Entrada Aleatória');
% xlabel('Index');
% ylabel('Valor Quadratura');
%-------------------------------------------------------------------------%

y = qammod(x,M,0);                                                         % Resultado modulado, mas com potência não unitária

potenciaMedia = mean(abs(y).^2);
S1n = y/sqrt(potenciaMedia);                                                % Resultado final, com potência unitária

%-------------------------------------------------------------------------%
%BPSK - Binary Phase Shift Keying

n = 100;                                                                   % Número de bits processados

vetorBinario = randi([0 1],n,1);                                           % Gerando vetor aleatório de n bits
S2n = 2*vetorBinario - 1;                                                   % Resultado final, com potência unitária