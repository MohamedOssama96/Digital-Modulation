%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Experiment 2: Channel coding
%Using 'LBC' ,  Repetition codes ,  Convolutional codes ,  Uncoded 
%==========================================================================
%% identifying Simulation parameters:
num = 1e6 ;               % Number of bits/SNR = 1 million bits
k = 4 ;                   % Block length
n = 7 ;                   % Chip length
SNR = 0:2:30 ;              % Signal to noise ratio range in dB
snr = 10.^(SNR/10) ;        % Signal to noise ratio
%initializing some empty matices to be filled with the sequences 
Encoded_data = zeros(1 , num*(n/k)) ;
decoded_data = zeros(1 , num) ;
Pe_LBC = zeros(1 , 16) ;
Pe_3 = zeros(1 , 16) ;
Pe_5 = zeros(1 , 16) ;
Pe_11 = zeros(1 , 16) ;
Pe_conv = zeros(1 , 16) ;
Pe_s = zeros(1 , 16) ;
%% (A) Linear Block Code
bits  =  randi([0 1] , 1 , num) ;     % Generate random binary data vector(1 x n)
%% Generator Matrix Implementation
pol  =  cyclpoly(n , k) ;                       % Generating polynomial matrix
ParityMatrix  =  cyclgen(n , pol) ;             % Generating parity matrix
GeneratorMatrix  =  gen2par(ParityMatrix) ;   % Generating matrix
trt  =  syndtable(ParityMatrix) ;             % Generating syndrome matrix
%% Data Encoding with LBC Technique
Encoded_signal  =  encode(bits , n , k , 'linear/binary' , GeneratorMatrix) ;%encoding the signal
signal_PRK = ((Encoded_signal.*2)-1)*(sqrt(k/n)) ; % BPSK Modulation and normalizing the energy 
l = length(signal_PRK) ;
%% RX 
for i = 1:16
    %applying noise on the sequence
    noise = (sqrt(1/(2*snr(i))))*(randn(1 , l)+(1j*randn(1 , l))) ;
    PRK_sequence = signal_PRK+noise ;
    % decision Making 
    PRK_RX = (real(PRK_sequence) >= 0) ;       % Comparing the data with the threshold value  
    Decoded_signal = decode(PRK_RX , n , k , 'linear/binary' , GeneratorMatrix , trt) ; % Data decoding 

    %Calculating the number of error bits and Pe in the sequence 
    [~ , Pe_LBC(i)] = biterr(bits , Decoded_signal) ; 
end

%% (B) Repetition Code
bits_repeated_3 = repelem(bits , 3) ;                        % repeating each bit 3 times
bits_repeated_5 = repelem(bits , 5) ;                        % repeating each bit 5 times
bits_repeated_11 = repelem(bits , 11) ;                      % repeating each bit 11 times
bits_repeated_3 = ((2*bits_repeated_3)-1)*(sqrt(1/3)) ;% modulating the bits and normalizing their energy 
bits_repeated_5 = ((2*bits_repeated_5)-1)*(sqrt(1/5)) ;% modulating the bits and normalizing their energy 
bits_repeated_11 = ((2*bits_repeated_11)-1)*(sqrt(1/11)) ;% modulating the bits and normalizing their energy

%% RX
%recovering the original bit stream and removing the noise effect 
%by comparing each bit with a certain threshold  =  1/2
for i = 1:16
    %% Applying the noise on the sequences
    noise = (sqrt(1/(2*snr(i))))*(randn(1 , 3*num)+(1j*randn(1 , 3*num))) ;
    Rx_sequence_3 = bits_repeated_3+noise ;
    noise = (sqrt(1/(2*snr(i))))*(randn(1 , 5*num)+(1j*randn(1 , 5*num))) ;
    Rx_sequence_5 = bits_repeated_5+noise ;
    noise = (sqrt(1/(2*snr(i))))*(randn(1 , 11*num)+(1j*randn(1 , 11*num))) ;
    Rx_sequence_11 = bits_repeated_11+noise ;
     
    %% recovering bits in case of n = 3
    Rx_bits_3 = (real(Rx_sequence_3) >= 0) ; %returning the sequence to 0s & 1s
    %%returning the sequence to it's normal size agian by checking every 3
    %%bits which bit is the dominant ;if the number of ones was the majority
    %%in the 3 bits then it was a 1 ;else ,  it was a zero
    Rx_bits_3 = reshape(Rx_bits_3 , 3 , num) ;
    s = sum(Rx_bits_3) ;
    Rx_3 = (s >= 2) ;
    % calculating the number of error bits and Pe in the sequence
    [~ , Pe_3(i)] = biterr(bits , Rx_3) ;
    
    %% recovering bits in case of n = 5
    Rx_bits_5 = (real(Rx_sequence_5) >= 0) ;%returning the sequence to 0s & 1s
    %%returning the sequence to it's normal size agian by checking every 5
    %%bits which bit is the dominant ;if the number of ones was the majority
    %%in the 5 bits then it was a 1 ;else ,  it was a zero
    Rx_bits_5 = reshape(Rx_bits_5 , 5 , num) ;
    s = sum(Rx_bits_5) ;
    Rx_5 = (s >= 3) ;
    % calculating the number of error bits and Pe in the sequence
    [~ , Pe_5(i)] = biterr(bits , Rx_5) ;
    
    %% recovering bits in case of n = 11
    Rx_bits_11 = (real(Rx_sequence_11) >= 0) ;%returning the sequence to 0s & 1s
    %%returning the sequence to it's normal size agian by checking every 11
    %%bits which bit is the dominant ;if the number of ones was the majority
    %%in the 11 bits then it was a 1 ;else ,  it was a zero
    Rx_bits_11 = reshape(Rx_bits_11 , 11 , num) ;
    s = sum(Rx_bits_11) ;
    Rx_11 = (s >= 6) ;
    % calculating the number of error bits and Pe in the sequence
    [~ , Pe_11(i)] = biterr(bits , Rx_11) ;
    
end

%% (C)Convolutional code
%encoding parameters:
%================================================= 
costlength = 9 ;             %number of shift registers
traceback = 5*costlength ;   %traceback
polynomial = [657 435] ;     %2 outputs then n = 2
% converting the polynomial matrix to a trellis object so it can be used to
% encode the signal
trellis = poly2trellis(costlength , polynomial) ; 

%Encoding the bit sequence:
%=================================================
encoded_bits = convenc(bits , trellis) ;
% modulation type: BPSK = PRK
%=================================================
Tx_sequence = ((encoded_bits.*2)-1)*sqrt(1/2) ;% modualting and normalizing the energy
l = length(Tx_sequence) ;
%% Rx
for i = 1:16
    %%Aplying the noise
    noise = (sqrt(1/(2*snr(i))))*(randn(1 , l)+(1j*randn(1 , l))) ;
    Rx_sequence = Tx_sequence+noise ;
    %demodulation:
    %==================================
    Rx_bits = (Rx_sequence >= 0) ; 
    %Decoding:
    %==================================
    decoded_bits = vitdec(Rx_bits , trellis , traceback , 'trunc' , 'hard') ;
    
    % calculating the number of error bits and Pe in the sequence
    [~ , Pe_conv(i)] = biterr(bits , decoded_bits) ;

end
%% (D)Simple detector (Uncoded)
Tx_sequence = ((bits.*2)-1) ; %bpsk modulation
for i = 1:16
    %applying the noise
    noise = (sqrt(1/(2*snr(i))))*(randn(1 , num)+(1j*randn(1 , num))) ;
    Rx_sequence = Tx_sequence+noise ;
    %recovering the bits
    Rx_bits = (Rx_sequence >= 0) ; 
    % calculating the number of error bits and Pe in the sequence
    [~ , Pe_s(i)] = biterr(bits , Rx_bits) ;
end
%% Ploting the BER curve against SNR for PRK
figure
semilogy(SNR , Pe_LBC , 'k-*' , 'LineWidth' , 2)
hold on ;
semilogy(SNR , Pe_3 , 'r-*' , 'LineWidth' , 2)
hold on ;
semilogy(SNR , Pe_5 , 'g-*' , 'LineWidth' , 2)
hold on ;
semilogy(SNR , Pe_11 , 'b-*' , 'LineWidth' , 2)
hold on ;
semilogy(SNR , Pe_conv , 'm-*' , 'LineWidth' , 2)
hold on ;
semilogy(SNR , Pe_s , 'c-*' , 'LineWidth' , 2)
title('BER vs. SNR')
ylabel('BER')
xlabel('SNR')
grid on ;
legend('LBC' , 'n = 3' , 'n = 5' , 'n = 11' , 'Conv' , 'Uncoded')
