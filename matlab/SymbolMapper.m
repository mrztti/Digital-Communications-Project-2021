% =========================================================================
% Symbol Mapper class
% =========================================================================
% Defines all given symbol constellations and their respective map/unmap
% functions.
% Usage:
%   sm = SymbolMapper.BPSK
%   x = sm.map(c)
%

classdef SymbolMapper
  properties
        mapping_function
        unmapping_function
        n
    end
    methods
        function sm = SymbolMapper(f, u, n)
            sm.mapping_function = f;
            sm.unmapping_function = u;
            sm.n = n;
        end

        function x = map(obj, c)
            x = obj.mapping_function(c);
        end

        function y = unmap(obj, symbols)
            y = obj.unmapping_function(symbols);
        end

        function c = constellation(obj)
            c = obj.mapping_function(binary_basis(obj.n));
        end
        function y = AWGN_channel(obj, x, snr, encoder)
            % Simulates the awgn channel by adding complex gaussian noise
            % x : symbol stream
            % snr : sound to noise ratio in dB
            
            x = x(:); % Make sure that we have a column vector
            n_x = length(x);
        
            % Calculate the std from the snr (only add noise if finite)
            if isfinite(snr)
                power = (x' * x) / (obj.n*n_x); % Calculate power of signal
                std = sqrt((power/encoder.rate).*10.^(-snr/10)); % calculate std from snr
        
                % Add complex gaussian noise to signal
                N = (std * 1/sqrt(2)) .* (randn(n_x, length(snr)) + 1i * randn(n_x, length(snr)));
                y = x + N;
            else 
                y = x;
            end
        end

    end
    enumeration
        BPSK(@BPSK_map, @BPSK_unmap, 1)
        QPSK_GRAY(@QPSKG_map, @QPSKG_unmap, 2)
        AMPM(@AMPM_map, @AMPM_unmap, 3)
    end
end


 function x = BPSK_map(c)
    % The BPSK mapping is equivalent to applying -2*c + 1
    x = 2 * c - 1;
 end

 function y = BPSK_unmap(symbols)
    y = (real(symbols) >= 0);
 end

function x = QPSKG_map(c)
    %   c : bit vector
    %
    %   We want the average energy per symbol to be 1 (Es = 1)
    %   => basis vectors are (j) sqrt(1/2)
    %

    % Securities (even length column)
    c = c(:); 
    if rem(length(c),2) ~= 0
            error('bit vector must be of even length');
    end
    
    a = sqrt(1/2);
    va = [a; 1i * a]; % Constant basis vectors
    c_gray = [c(1:2:end), c(2:2:end)]; % Split into length 2 chuncks
    x = - BPSK_map(c_gray) * va; % Re-use BPSK to obtain value along basis
end

function y = QPSKG_unmap(symbols)
    y = zeros(2*length(symbols), 1);
    s1 = (real(symbols) < 0);
    s2 = (imag(symbols) < 0);
    y(1:2:end) = s1;
    y(2:2:end) = s2;
end

function b = binary_basis(n)
    b = dec2bin(0:2^n-1)' - '0';
end

function x = AMPM_map(c)
    %   c : bit vector
    %
    %   We want the average energy per symbol to be 1 (Es = 1)
    %   => alpha = 4 / ( sqrt(2)*(4 + 2*sqrt(5)) ) = 0.3339 ~ 1/3
    %

    
    % Securities (mod 3 length column)
    c = c(:); 
    if rem(length(c),3) ~= 0
            error('bit vector must be a multiple of 3');
    end

    a = 4 / ( sqrt(2)*(4 + 2*sqrt(5)) ); % Constant defined by unit average energy
    % Split into length 3 chunks and apply 
    c_3b = - BPSK_map([c(1:3:end), c(2:3:end), c(3:3:end)]); 
    
    % Diagonal separation (1st bit)
    out = c_3b(:,1) * (-a + 1i*a);
    
    % Quadrant separation across imaginary axis (3rd bit)
    out = out + c_3b(:, 3) * 2*a;

    % Quadrant separation across real axis (2nd bit)
    x = out + c_3b(:,2).*(-2*a*1i * sign(real(out)));

end

function y = AMPM_unmap(symbols)
    s3 = (real(symbols) < 0);
    s2 = (imag(symbols).*sign(s3 - 0.5) < 0); % flipped rule for different quandrants

    a = 1 / sqrt(10);
    rs = sign((real(symbols) > 0) - 0.5);
    is = sign((imag(symbols) > 0) - 0.5);
    s1 = (real(symbols - 2 * a * rs - 1i * 2 * a * is) > 0);

    y = zeros(3*length(symbols), 1);
    y(1:3:end) = s1;
    y(2:3:end) = s2;
    y(3:3:end) = s3;
end



























