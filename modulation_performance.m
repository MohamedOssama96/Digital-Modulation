%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Experiment 1: Performance of Different Modultion Types
%==========================================================================
%% identifying Simulation parameters:
n = 1e6 ;                   % Number of bits/SNR = 1 million bits
SNR = 0:2:30 ;              % Signal to noise ratio range in dB
snr = 10.^(SNR/10) ;        % Signal to noise ratio
%initializing some empty matices to be filled with the sequences 
Pe_OOK = zeros(1 , 16) ;
Pe_PRK = zeros(1 , 16) ;
Pe_FSK = zeros(1 , 16) ;
Pe_OOK_tb = zeros(1 , 16) ;
Pe_PRK_tb = zeros(1 , 16) ;
Pe_FSK_tb = zeros(1 , 16) ;
Pe_QAM_tb = zeros(1 , 16) ;
%% TX 
bits = randi([0 , 1] , 1 , n) ;      %Generate random binary data vector(1 x n)
%OOK_Modulation(0 , 1)
signal_OOK = bits ; %No change in the bits will be required
%PRK modulation(-1 , 1)
signal_PRK = (bits.*2)-1 ;
%FSK modulation
%if bit to send = 0 send 1 , else send i
ones = find(bits) ;     %find all elements that are equal to one in the bits sequence
zeros = find(~bits) ;   %find all elements that are equal to one in the bits sequence
signal_FSK(ones) = i ;  %replace all the ones with i
signal_FSK(zeros) = 1 ; %replace all the zeros with 1

%% modulating using toolboxes
signal_OOK_tb = genqammod(bits , [0  1]) ;   %OOK modulation using toolbox 
signal_PRK_tb = pskmod(bits , 2) ;           %PRK modulation using toolbox 
signal_FSK_tb = genqammod(bits , [1  1i]) ;  %FSK modulation using toolbox 
M = 16 ;
QAM_TX = randi([0 M-1] , 1 , n) ;              %Generate random data vector(1 x n)
signal_QAM_tb = qammod(QAM_TX , 16) ;        %QAM modulation using toolbox

%% RX
%Applying noise to the bit stream with different values of SNR , receiving
%the bits stream , recovering the original bits by demodulating it 
%and calculating the BER and Pe for each type of modulated signal
for i = 1:16  
    %%applying the noise on OOK sequence
    noise = (sqrt(mean(abs(signal_OOK).^2)/(2*snr(i))))*(randn(1 , n)+(1j*randn(1 , n))) ;
    OOK_sequence = signal_OOK+noise ;
    %%applying the noise on PRK sequence
    noise = (sqrt(mean(abs(signal_PRK).^2)/(2*snr(i))))*(randn(1 , n)+(1j*randn(1 , n))) ;
    PRK_sequence = signal_PRK+noise ;
    %%applying the noise on FSK sequence
    noise = (sqrt(mean(abs(signal_FSK).^2)/(2*snr(i))))*(randn(1 , n)+(1j*randn(1 , n))) ;
    FSK_sequence = signal_FSK+noise ;
   
    %%aplying noise on toolboxes signals
    %%applying the noise on OOK_tb sequence
    noise = (sqrt(mean(abs(signal_OOK_tb).^2)/(2*snr(i))))*(randn(1 , n)+(1j*randn(1 , n))) ;
    OOK_sequence_tb = signal_OOK_tb+noise ;
    %%applying the noise on PRK_tb sequence
    noise = (sqrt(mean(abs(signal_PRK_tb).^2)/(2*snr(i))))*(randn(1 , n)+(1j*randn(1 , n))) ;
    PRK_sequence_tb = signal_PRK_tb+noise ;
    %%applying the noise on FSK_tb sequence
    noise = (sqrt(mean(abs(signal_FSK_tb).^2)/(2*snr(i))))*(randn(1 , n)+(1j*randn(1 , n))) ;
    FSK_sequence_tb = signal_FSK_tb+noise ;
    %%applying the noise on QAM_tb sequence
    noise = (sqrt(mean(abs(signal_QAM_tb).^2)/(2*snr(i))))*(randn(1 , n)+(1j*randn(1 , n))) ;
    QAM_sequence_tb = signal_QAM_tb+noise ;
    
    %% demodulation
    %recovering the original bit stream and removing the noise effect 
    %by comparing each bit with a threshold  = 
    %0.5 for OOK
    OOK_RX = (real(OOK_sequence) >= 0.5) ;
    %0 for PRK
    PRK_RX = (real(PRK_sequence) >= 0) ;
    %For Fsk we will check for each received point ;if its real component
    %is bigger than its imaginary one then the point was a zero ;else ,  it
    %would be a one
    FSK_RX = (real(FSK_sequence)<imag(FSK_sequence)) ;
    
    %%demodulation using toolboxes
    OOK_RX_tb = genqamdemod(OOK_sequence_tb , [0  1]) ;
    PRK_RX_tb = pskdemod(PRK_sequence_tb , 2) ;
    FSK_RX_tb = genqamdemod(FSK_sequence_tb , [1  1i]) ;
    QAM_RX_tb = qamdemod(QAM_sequence_tb , M) ;
    
    %% BER calculations 
    %Calculating the number of error bits in the sequence 
    % and the probability of error in each sequence
    %OOK
    [~ , Pe_OOK(i)] = biterr(bits , OOK_RX) ;
    %PRK 
    [~ , Pe_PRK(i)] = biterr(bits , PRK_RX) ;  
    % FSK
    [~ , Pe_FSK(i)] = biterr(bits , FSK_RX) ;
    %%for toolboxes
    [~ , Pe_OOK_tb(i)] = biterr(bits , OOK_RX_tb) ;
    [~ , Pe_PRK_tb(i)] = biterr(bits , PRK_RX_tb) ;
    [~ , Pe_FSK_tb(i)] = biterr(bits , FSK_RX_tb) ;
    [~ , Pe_QAM_tb(i)] = symerr(QAM_TX , QAM_RX_tb) ;
end
% %% Ploting the BER curve against SNR
%% Ploting all in one curve
figure
semilogy(SNR , Pe_OOK , 'r-*' , 'LineWidth' , 2)
hold on ;
semilogy(SNR , Pe_OOK_tb , 'm--o' , 'LineWidth' , 2)
hold on ;
semilogy(SNR , Pe_PRK , 'g-*' , 'LineWidth' , 2)
hold on ;
semilogy(SNR , Pe_PRK_tb , 'k--+' , 'LineWidth' , 2)
hold on ;
semilogy(SNR , Pe_FSK , 'b-*' , 'LineWidth' , 2)
hold on ;
semilogy(SNR , Pe_FSK_tb , 'y--s' , 'LineWidth' , 2)
hold on ;
semilogy(SNR , Pe_QAM_tb , 'c-p' , 'LineWidth' , 2)
title('BER vs. SNR ')
ylabel('BER')
xlabel('SNR')
legend('OOK' , 'OOK-toolbox' , 'PRK-toolbox' , 'PRK' , 'FSK' , 'FSK-toolbox' , 'QAM-toolbox')
grid on ;
