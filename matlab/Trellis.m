classdef Trellis
    properties
        numInputSymbols % Number of input symbols
        numOutputSymbols % Number of output symbols
        numStates % Number of states
        nextStates % Next state matrix
        outputs % Output matrix
        inputBasis
        outputBasis
    end
    methods
        function t = Trellis(tr)
            t.numInputSymbols = tr.numInputSymbols;
            t.numOutputSymbols = tr.numOutputSymbols;
            t.numStates = tr.numStates;
            t.nextStates = tr.nextStates;
            t.outputs = tr.outputs;
            t.inputBasis = dec2bin(0:2^log2(tr.numInputSymbols)-1)' - '0';
            t.outputBasis = dec2bin(0:2^log2(tr.numOutputSymbols)-1)' - '0';
        end
        function [new_state_index, output_index] = getNextState(obj, curr_index, input_index)
            new_state_index = oct2dec(obj.nextStates(curr_index, input_index)) + 1;
            output_index = obj.outputs(curr_index, input_index) + 1;
        end

        function y = getOutput(obj, x)
            n = log2(obj.numOutputSymbols);
            y = zeros(n,1);
            y_p = oct2poly(x-1);
            y((n-length(y_p)+1):end) = y_p;
        end

        function y = inputIdx2seq(obj, idx)
            split = arrayfun(@(x) obj.inputBasis(:, x), idx);
            y = multiplex(log2(obj.numInputSymbols), split);
        end

    end
    enumeration
        % constraint length = 3
        % first poly [1 0 1] = 5
        % second poly [1 1 1] = 7

        E1(generate_rate_half_trellis(3, [1,0,1;1,1,1])),

        % constraint length = 5
        % first poly [1 0 1 1 1] = 27 (octal!)
        % second poly [1 0 1 1 0] = 26

        E2(generate_rate_half_trellis(5, [1 0 1 1 1;1 0 1 1 0])),

        % constraint length = 5
        % first poly [1 0 0 1 1] = 23
        % second poly [1 1 0 1 1] = 33

        E3(generate_rate_half_trellis(5, [1 0 0 1 1;1 1 0 1 1]))

        %TODO E

    end
end

function tr = generate_rate_half_trellis(constraint_length, G)
    n = 2;
    state_size = constraint_length-1;
    state_number = 2^(state_size);
    all_states = (dec2bin(0:2^state_size-1)' - '0')';
    all_transitions = [0;1];
    state_matrix = zeros(state_number,2);
    output_matrix = zeros(state_number,2);

    for s_i = 1:state_number
        for inp_i = 1:2
            curr_state = all_states(s_i,:);
            p2 = poly2oct([all_transitions(inp_i,:), curr_state(:,1:end-1)]);
            state_matrix(s_i,inp_i) = p2;
            input = [all_transitions(inp_i, :), curr_state];
            output_matrix(s_i,inp_i) = poly2oct(mod(input *(G'), 2));
        end
    end

    tr.outputs = output_matrix;
    tr.nextStates = state_matrix;
    tr.numInputSymbols = 2;
    tr.numOutputSymbols = 2^n;
    tr.numStates = state_number;
end

function tr = generate_E4_trellis()
    n = 2^2;
    state_size = 3;
    state_number = 2^(state_size);
    all_states = (dec2bin(0:2^state_size-1)' - '0')';
    all_transitions = [0;1];
    state_matrix = zeros(state_number,2);
    output_matrix = zeros(state_number,2);

    for s_i = 1:state_number
        for inp_i = 1:2
            curr_state = all_states(s_i,:);
            p2 = poly2oct([all_transitions(inp_i,:), curr_state(:,1:end-1)]);
            state_matrix(s_i,inp_i) = p2;
            input = [all_transitions(inp_i, :), curr_state];
            output_matrix(s_i,inp_i) = poly2oct(input * G');
        end
    end

    tr.outputs = output_matrix;
    tr.nextStates = state_matrix;
    tr.numInputSymbols = 1;
    tr.numOutputSymbols = n;
    tr.numStates = state_number;
end

function n = poly2oct(poly)
    n = 0;
    v = fliplr(poly);
    for i = 1:size(poly,2)
        n = n + v(i)*2^(i-1);
    end
    n = str2num(dec2base(n, 8));
end

% Multiplex [y1;y2;...] into [y11;y21;...;y12;y22;...] 
function m = multiplex(rate, stream)
    N = size(stream, 1) * rate;
    m = zeros(1, N);
    for i = 1:rate
        m(i:rate:end) = stream(:, i);
    end    
end






