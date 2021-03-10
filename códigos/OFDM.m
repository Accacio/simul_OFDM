close all
clear
clc

%%
%--------------------------------------------------------------------------
% Variables
%--------------------------------------------------------------------------

numbBits = input('Number of bits [128 / 64] = ');                          % Number of bits to be transmitted
mOrder = input('Modulation order [4 / 2] = ');                             % Modulation order    
M = 64;                                                                    % Decimation Factor
L = 41;                                                                    % Size of channels
N = M + L - 1;
vSNR = -15:0.1:15;                                                         % SNR range values (dB)
n = 10000;                                                                 % Number of channels

%%
%--------------------------------------------------------------------------
% Transmitter
%--------------------------------------------------------------------------

data = randi([0 1],1,numbBits);                                            % Random data generation

%%
% Modulation --------------------------------------------------------------

if (mOrder==4) 
    inputData = qam4_mod(data,mOrder);                                     % 4-QAM modulation
else
    inputData = 2*data - 1;                                                % BPSK modulation
end

%%
% IDFT --------------------------------------------------------------------

DFT = dftmtx(M);
IDFT = (conj(DFT))/M;

idftData = IDFT*(inputData.');                                             % IDFT

%%
% Add Cyclic Prefix -------------------------------------------------------

cPrefixData = [eye(M); [eye(L-1) zeros(L-1, M-L+1)]]*idftData;             % Add cyclic prefix

%%
%--------------------------------------------------------------------------
% Channel
%--------------------------------------------------------------------------

nBER = 0;
totalBER = zeros(length(vSNR));

for noise = 1:length(vSNR)

    SNR=vSNR(noise);  
      
    for channel = 1:n

        h = sqrt(1/2)*(randn(1,L)+1i*randn(1,L));

        H = zeros(M,N);

        for i = 1:M
            H(i,i:L+i-1) = h;
        end

        y = H*cPrefixData;
        yData = awgn(y,SNR,'measured');


        %%
        %------------------------------------------------------------------
        % Receiver
        %------------------------------------------------------------------

        %%
        % DFT -------------------------------------------------------------

        dftData = DFT*yData;

        %%
        % Equalization ----------------------------------------------------

        lambda = sqrt(M)*IDFT*[h.'; zeros(M-L,1)];
        delta = diag(lambda);

        data_e = (delta^-1)*dftData;

        %%
        % Detection -------------------------------------------------------

        data_re = real(data_e);
        data_im = imag(data_e);

        data_a = zeros(numbBits/2,1);

        if mOrder==4
            if (data_re < 0 && data_im < 0)
                data_a = -0.7071 - 0.7071i;
            elseif (data_re >= 0 && data_im >= 0)
                data_a = 0.7071 + 0.7071i;
            elseif (data_re < 0 && data_im >= 0)
                data_a = -0.7071 + 0.7071i;
            elseif (data_re >= 0 && data_im < 0)
                data_a = 0.7071 - 0.7071i;
            end
        else
            if (data_re < 0)
                data_a = -1;
            else
                data_a = 1;
            end
        end

        %%
        % Demodulation ----------------------------------------------------

        if (mOrder==4) 
            outputData = qam4_demod(data_a,mOrder);                        % 4-QAM demodulation
        else
            outputData = (data_a.' + 1)/2;                                 % BPSK demodulation
        end

        [ne,nBER] = biterr(outputData, data);
        totalBER(noise) = totalBER(noise) + nBER;
    
    end % Canal Loop

end % SNR Loop

totalBER = totalBER/n; 

%%
%--------------------------------------------------------------------------
% Analysis
%--------------------------------------------------------------------------

% Bit Error Rate ----------------------------------------------------------
figure(1)
semilogy(vSNR,totalBER, '-b')
title('Probabilidade de Erro x SNR')
xlabel('SNR (dB)')
ylabel('BER')
axis([-15 15 0.001 1])
grid on
