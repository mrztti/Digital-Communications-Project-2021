% =========================================================================
% Symbol Mapper class
% =========================================================================
% Defines all given symbol constellations
% Usage:
%   sm = SymbolMapper.BPSK
%   x = sm.map(c)
%

classdef SymbolMapper
  properties
        mapping_function
        n
    end
    methods
        function sm = SymbolMapper(f, n)
            sm.mapping_function = f;
            sm.n = n;
        end

        function x = map(obj, c)
            x = obj.mapping_function(c);
        end

        function c = constellation(obj)
            c = obj.mapping_function(binary_basis(obj.n));
        end

    end
    enumeration
        BPSK(@BPSK_map, 1),
        QPSK_GRAY(@QPSKG_map, 2),
        AMPM(@AMPM_map, 3)
    end
end


 function x = BPSK_map(c)
    % The BPSK mapping is equivalent to applying -2*c + 1
    x = -2 * c + 1;
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
    x = BPSK_map(c_gray) * va; % Re-use BPSK to obtain value along basis
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

    a = 4 / ( sqrt(2)*(4 + 2*sqrt(5)) );
    % Split into length 3 chunks and apply 
    c_3b = BPSK_map([c(1:3:end), c(2:3:end), c(3:3:end)]); 
    
    % Diagonal separation (1st bit)
    out = c_3b(:,1) * (-a + 1i*a);
    
    % Quadrant separation across imaginary axis (3rd bit)
    out = out + c_3b(:, 3) * 2*a;

    % Quadrant separation across real axis (2nd bit)
    x = out + c_3b(:,2).*(-2*a*1i * sign(real(out)));

end


