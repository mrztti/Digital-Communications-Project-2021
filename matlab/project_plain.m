% ======================================================================= %
% SSY125 Project
% ======================================================================= %
clc
clear

% ======================================================================= %
% Simulation Options
% ======================================================================= %
N = 3e2;  % simulate N bits each transmission (one block)
maxNumErrs = 100; % get at least 100 bit errors (more is better)
maxNum = 1e3; % OR stop if maxNum bits have been simulated
EbN0 = [-10, 0, 10, 20, 30];%-1:8; % power efficiency range

% ======================================================================= %
% Other Options
% ======================================================================= %
constellation = SymbolMapper.QPSK_GRAY; % Choice of constellation


% ======================================================================= %
% Simulation Chain
% ======================================================================= %
BER = zeros(1, length(EbN0)); % pre-allocate a vector for BER results

parfor i = 1:length(EbN0) % use parfor ('help parfor') to parallelize
  totErr = 0;  % Number of errors observed
  num = 0; % Number of bits processed
  snr = EbN0(i);
  figure("Name", "Symbols received for EbN0 = " + string(snr));

  while((totErr < maxNumErrs) && (num < maxNum))
  % ===================================================================== %
  % Begin processing one block of information
  % ===================================================================== %
  % [SRC] generate N information bits 
  u = randi([0,1], N, 1);

  % [ENC] convolutional encoder
  c = u; % TODO

  % [MOD] symbol mapper  
  x = constellation.map(c);

  % [CHA] add Gaussian noise
  y = AWGN_channel(x, snr);

  % scatterplot: plot(y, 'b.')
  
  plt = plot_symbols(y, constellation, snr);
  % [HR] Hard Receiver
  % ...

  % [SR] Soft Receiver
  % ...
  % ===================================================================== %
  % End processing one block of information
  % ===================================================================== %
  BitErrs = 0; % count the bit errors and evaluate the bit error rate
  totErr = totErr + BitErrs;
  num = num + N; 

  disp(['+++ ' num2str(totErr) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr/num, '%10.1e') '. +++']);
  end 
  BER(i) = totErr/num; 
end
% ======================================================================= %
% End
% ======================================================================= %



% ======================================================================= %
% Custom functions
% ======================================================================= %

function p = plot_symbols(y, constellation, snr)
    hold on;
    grid on;
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






