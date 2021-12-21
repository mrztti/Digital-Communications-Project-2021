% =========================================================================
% DecoderType class
% =========================================================================
% Defines the two types of Viterbi decoders (Hard and Soft)
%

classdef DecoderType
    properties
        decoder
    end
    methods
        function dt = DecoderType(decoder)
            dt.decoder = decoder;
        end

        function y = decode(obj, symbols, trellis, constellation)
            y = obj.decoder(symbols, trellis, constellation);
        end
    end
    enumeration
        HARD(@hard_decoder)
        SOFT(@soft_decoder)
    end
end


%==========================================================================
% HARD DECODER VITERBI
%==========================================================================
function out = hard_decoder(symbols, trellis, constellation)
    num_states = trellis.numStates;
    num_inputs = trellis.numInputSymbols;
    num_out = trellis.numOutputSymbols;
    n = log2(num_out);

    % Split sequence into separate sequences for each convolutional output
    sequences = de_multiplex(n, constellation.unmap(symbols))';
    L = size(sequences, 2);
    
    % Store all the current best distances
    cumulative_metrics = zeros(num_states, L+1) + inf;
    cumulative_metrics(1) = 0;

    % Store the index of the best previous node before this one
    best_path_index = zeros(num_states, L+1);

    % Store the input from the best previous node to this one
    best_path_input = zeros(num_states, L+1);
    
    
    %%% FEED FORWARD
    
    % Iterate until the end of the sequence
    for time = 1:L
    
        % Desired sequence
        seq = sequences(:, time);
        metrics = sum(abs(trellis.specificOutputs - seq),1);

        % Iterate over all states and inputs
        for s_i = 1:num_states
            for inp_i = 1:num_inputs
                next_state = trellis.convertedNextStates(s_i, inp_i);
                next_min = metrics(:, s_i, inp_i) + cumulative_metrics(s_i, time);
                if next_min < cumulative_metrics(next_state,time+1)
                    cumulative_metrics(next_state,time+1) = next_min;
                    best_path_index(next_state, time+1) = s_i;
                    best_path_input(next_state,time+1) = inp_i;
                end
            end
        end
    end

    %%% PROPAGATE BACK
    idx = zeros(L, 1);

    if isinf(cumulative_metrics(1, L))
        error("Symbol stream did not terminate in the all-zero state")
    end
    
    current_state = 1;

    % Go back in time
    for time = fliplr(1:L)
        
        % Follow the best path
        idx(time) = best_path_input(current_state, time+1);
        current_state = best_path_index(current_state, time+1);
    end
    % Convert inputs from idx to seq
    out = trellis.inputIdx2seq(idx)';
end


%==========================================================================
% SOFT DECODER VITERBI (TODO)
%==========================================================================
function out = soft_decoder(symbols, trellis, constellation)
    num_states = trellis.numStates;
    num_inputs = trellis.numInputSymbols;
    num_out = trellis.numOutputSymbols;
    n = log2(num_out);

    mapped_outputs = reshape(constellation.map(trellis.specificOutputs), num_states, num_inputs);
    L = length(symbols);
    
    % Store all the current best distances
    cumulative_metrics = zeros(num_states, L+1) + inf;
    cumulative_metrics(1) = 0;

    % Store the index of the best previous node before this one
    best_path_index = zeros(num_states, L+1);

    % Store the input from the best previous node to this one
    best_path_input = zeros(num_states, L+1);
    
    
    %%% FEED FORWARD
    
    % Iterate until the end of the sequence
    for time = 1:L
    
        % Desired sequence
        s = symbols(time,:);
        metrics = abs(mapped_outputs-s).^2;

        % Iterate over all states and inputs
        for s_i = 1:num_states
            for inp_i = 1:num_inputs
                next_state = trellis.convertedNextStates(s_i, inp_i);
                next_min = metrics(s_i, inp_i) + cumulative_metrics(s_i, time);
                if next_min < cumulative_metrics(next_state,time+1)
                    cumulative_metrics(next_state,time+1) = next_min;
                    best_path_index(next_state, time+1) = s_i;
                    best_path_input(next_state,time+1) = inp_i;
                end
            end
        end
    end

    %%% PROPAGATE BACK
    idx = zeros(L, 1);

    if isinf(cumulative_metrics(1, L))
        error("Symbol stream did not terminate in the all-zero state")
    end
    
    current_state = 1;

    % Go back in time
    for time = fliplr(1:L)
        
        % Follow the best path
        idx(time) = best_path_input(current_state, time+1);
        current_state = best_path_index(current_state, time+1);
    end
    % Convert inputs from idx to seq
    out = trellis.inputIdx2seq(idx)';
end

function m = de_multiplex(rate, stream)
    stream = stream(:);
    N = length(stream) / rate;
    m = zeros(N, rate);
    for i = 1:rate
        m(:, i) = stream(i:rate:end);
    end    
end

