% ======================================================================= %
% SSY125 Project
% ======================================================================= %
clc
clear

% ======================================================================= %
% Simulation Options
% ======================================================================= %
N = 1e4;  % 5 simulate N bits each transmission (one block)
maxNumErrs = 100; % get at least 100 bit errors (more is better)
maxNum = 1e5; % 6 OR stop if maxNum bits have been simulated
EbN0 = -1:8; % power efficiency range:

% ======================================================================= %
% Other Options
% ======================================================================= %
constellation = SymbolMapper.QPSK_GRAY; % Choice of constellation
convolutional_encoder = ConvEncoder.E1; % Choice of convolutional code
decoder_type = DecoderType.HARD; % Choice of HARD/SOFT decoding

decoder = ViterbiDecoder(convolutional_encoder.trellis, decoder_type, constellation);

% ======================================================================= %
% Simulation Chain
% ======================================================================= %
BER_coded = zeros(1, length(EbN0)); % pre-allocate a vector for BER results
BER_uncoded = zeros(1, length(EbN0));

for i = 1:length(EbN0) % use parfor ('help parfor') to parallelize
  totErr_u = 0;  % Number of uncoded errors observed
  totErr_c = 0;  % Number of coded errors observed
  num = 0; % Number of bits processed
  snr = EbN0(i);
  
  drawFirst = true;
  while((totErr_c < maxNumErrs) && (num < maxNum))
  % ===================================================================== %
  % Begin processing one block of information
  % ===================================================================== %
  % [SRC] generate N information bits 
  u = randi([0,1], N, 1);

  % [ENC] convolutional encoder
  c = convolutional_encoder.encode(u);

  % [MOD] symbol mapper  
  x_coded = constellation.map(c);
  x_uncoded = constellation.map(u);

  % [CHA] add Gaussian noise
  y_coded = AWGN_channel(x_coded, snr);
  y_uncoded = AWGN_channel(x_uncoded, snr);

  % Only draw on the first iteration
%   if drawFirst
%       figure("Name", "Symbols received for EbN0 = " + string(snr));
%       hold on;
%       grid on;
%       plt = plot_symbols(y, constellation, snr);
%   end

  cf_coded = decoder.decode(y_coded);
  cf_uncoded = constellation.unmap(y_uncoded);
  % ===================================================================== %
  % End processing one block of information
  % ===================================================================== %
  BitErrs_coded = sum(u~=cf_coded); % count the bit errors and evaluate the bit error rate
  BitErrs_uncoded = sum(u~=cf_uncoded);
  totErr_u = totErr_u + BitErrs_uncoded;
  totErr_c = totErr_c + BitErrs_coded;
  num = num + N; 

  disp(['+++ [' num2str(totErr_u) '] ' num2str(totErr_c) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr_c/num, '%10.1e') '. +++']);
  end 
  BER_coded(i) = totErr_c/num; 
  BER_uncoded(i) = totErr_u/num;
end
% ======================================================================= %
% End
% ======================================================================= %

% ======================================================================= %
% Plot results
% ======================================================================= %

figure()
hold on;

BER_theory = @(EbN0) 2*qfunc(sqrt(2*EbN0)) - (qfunc(sqrt(2*EbN0))).^2;
plot(EbN0,BER_theory(10.^(EbN0 / 10)),'--', 'Color','Black')
plot(EbN0, BER_uncoded, 'Color', 'Red')
plot(EbN0, BER_coded, 'Color', 'Blue')

title('Plot of coded and uncoded BER compared to the theoretical BER')
xlabel('E_b/N_0 [dB]')
ylabel('BER')
legend('Theoretical BER', 'Uncoded transmission', 'Coded transmission')
axis([EbN0(1) EbN0(end) 1e-4 1])
set(gca, 'YScale', 'log')


% ======================================================================= %
% Custom functions
% ======================================================================= %

function p = plot_symbols(y, constellation, snr)
    col = linspace(1,10,length(y));
    scatter(real(y), imag(y), col);
    cs = constellation.constellation();
    p = scatter(real(cs), imag(cs), 'red', 'filled');
end

function y = AWGN_channel(x, snr)
    % Simulates the awgn channel by adding complex gaussian noise
    % x : symbol stream
    % snr : sound to noise ratio in dB
    
    x = x(:); % Make sure that we have a column vector
    n_x = length(x);

    % Calculate the std from the snr (only add noise if finite)
    if isfinite(snr)
        power = (x' * x) / n_x; % Calculate power of signal
        std = sqrt(power.*10.^(-snr/10)); % calculate std from snr

        % Add complex gaussian noise to signal
        N = (std * 1/sqrt(2)) .* (randn(n_x, length(snr)) + 1i * randn(n_x, length(snr)));
        y = x + N;
    else 
        y = x;
    end
end






