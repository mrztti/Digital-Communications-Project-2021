% ======================================================================= %
% SSY125 Project
% ======================================================================= %
clc
clear

% ======================================================================= %
% Simulation Options
% ======================================================================= %
N = 3e5;  % 5 simulate N bits each transmission (one block)
maxNumErrs = 100; % get at least 100 bit errors (more is better)
maxNum = 3e6; % 6 OR stop if maxNum bits have been simulated
EbN0 = -1:8; % power efficiency range:

% ======================================================================= %
% Other Options
% ======================================================================= %
decoder_type = DecoderType.SOFT;

cons1 = SymbolMapper.BPSK;
enc1 = ConvEncoder.E3;
dec1 = ViterbiDecoder(enc1.trellis, decoder_type, cons1);

cons2 = SymbolMapper.QPSK_GRAY;
enc2 = ConvEncoder.E3;
dec2 = ViterbiDecoder(enc2.trellis, decoder_type, cons2);

cons3 = SymbolMapper.AMPM;
enc3 = ConvEncoder.E4;
dec3 = ViterbiDecoder(enc3.trellis, decoder_type, cons3);

% ======================================================================= %
% Simulation Chain
% ======================================================================= %
BER_coded1 = zeros(1, length(EbN0)); % pre-allocate a vector for BER results
BER_coded2 = zeros(1, length(EbN0));
BER_coded3 = zeros(1, length(EbN0));

BER_uncoded1 = zeros(1, length(EbN0)); % pre-allocate a vector for BER results
BER_uncoded2 = zeros(1, length(EbN0));
BER_uncoded3 = zeros(1, length(EbN0));

lb = LoadingBar(length(EbN0)*maxNum);

for i = 1:length(EbN0) % use parfor ('help parfor') to parallelize
  totErr1 = 0;  % Number of coded errors observed
  totErr2 = 0;  % Number of coded errors observed
  totErr3 = 0;  % Number of coded errors observed

  totErr1uc = 0;  % Number of coded errors observed
  totErr2uc = 0;  % Number of coded errors observed
  totErr3uc = 0;  % Number of coded errors observed

  num = 0; % Number of bits processed
  snr = EbN0(i);
  
  drawFirst = true;
  while((totErr1 + totErr2 + totErr3 < 3*maxNumErrs) && (num < maxNum))
  % ===================================================================== %
  % Begin processing one block of information
  % ===================================================================== %
  % [SRC] generate N information bits 
  u = randi([0,1], N, 1);

  % [ENC] convolutional encoder
  c1 = enc1.encode(u);
  c2 = enc2.encode(u);
  c3 = enc3.encode(u);

  % [MOD] symbol mapper  
  x1 = cons1.map(c1);
  x2 = cons2.map(c2);
  x3 = cons3.map(c3);
  x1uc = cons1.map(u);
  x2uc = cons2.map(u);
  x3uc = cons3.map(u);

  % [CHA] add Gaussian noise
  y1 = cons1.AWGN_channel(x1, snr, enc1);
  y2 = cons2.AWGN_channel(x2, snr, enc2);
  y3 = cons3.AWGN_channel(x3, snr, enc3);

  y1uc = cons1.AWGN_channel(x1uc, snr, ConvEncoder.NONE);
  y2uc = cons2.AWGN_channel(x2uc, snr, ConvEncoder.NONE);
  y3uc = cons3.AWGN_channel(x3uc, snr, ConvEncoder.NONE);

  % Only draw on the first iteration
%   if drawFirst
%       figure("Name", "Symbols received for EbN0 = " + string(snr));
%       hold on;
%       grid on;
%       plt = plot_symbols(y, constellation, snr);
%   end

  cf1 = dec1.decode(y1);
  cf2 = dec2.decode(y2);
  cf3 = dec3.decode(y3);

  cf1uc = cons1.unmap(y1uc);
  cf2uc = cons2.unmap(y2uc);
  cf3uc = cons3.unmap(y3uc);
  % ===================================================================== %
  % End processing one block of information
  % ===================================================================== %
  BitErrs1= sum(u~=cf1); % count the bit errors and evaluate the bit error rate
  BitErrs2 = sum(u~=cf2);
  BitErrs3 = sum(u~=cf3);
  BitErrs1uc = sum(u~=cf1uc); % count the bit errors and evaluate the bit error rate
  BitErrs2uc = sum(u~=cf2uc);
  BitErrs3uc = sum(u~=cf3uc);


  totErr1 = totErr1 + BitErrs1;
  totErr2 = totErr2 + BitErrs2;
  totErr3 = totErr3 + BitErrs3;

  totErr1uc = totErr1uc + BitErrs1uc;
  totErr2uc = totErr2uc + BitErrs2uc;
  totErr3uc = totErr3uc + BitErrs3uc;
  num = num + N; 
  lb = lb.step(N);
  end 

  BER_coded1(i) = totErr1/num;
  BER_coded2(i) = totErr2/num; 
  BER_coded3(i) = totErr3/num;

  BER_uncoded1(i) = totErr1uc/num;
  BER_uncoded2(i) = totErr2uc/num; 
  BER_uncoded3(i) = totErr3uc/num;
  lb = lb.set(i*maxNum);
end
% ======================================================================= %
% End
% ======================================================================= %

% ======================================================================= %
% Plot results
% ======================================================================= %

figure()
hold on;
plot(EbN0, BER_coded1, 'Color', 'Red')
plot(EbN0, BER_coded2, 'Color', 'Green')
plot(EbN0, BER_coded3, 'Color', 'Cyan')

plot(EbN0, BER_uncoded1, 'Color', '#AB0101', 'Marker','x')
plot(EbN0, BER_uncoded2, 'Color', '#01AB31', 'Marker','x')
plot(EbN0, BER_uncoded3, 'Color', '#1C01AB', 'Marker','x')


title('Plot of different systems with/without coding')
xlabel('E_b/N_0 [dB]')
ylabel('BER')
legend('System 1 - CODED', 'System 2 - CODED', 'System 3 - CODED', 'System 1 - Uncoded', 'System 2 - Uncoded', 'System 3 - Uncoded')
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








