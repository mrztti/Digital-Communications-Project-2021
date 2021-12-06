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
        seq_desired = sequences(:, time);

        % Iterate overall all possible states 
        for si = 1:num_states
            % Skip unreachable states
            if isinf(cumulative_metrics(si, time))
                continue
            end

            % Iterate over all possible inputs
            for inp = 1:num_inputs
                % compute the new state and the new
                [new_state, output] = trellis.getNextState(si, inp);
                seq_out = trellis.getOutput(output);
                dist = sum(abs(seq_out - seq_desired));
                past_metric = cumulative_metrics(si, time);

                % See if dist is smaller than any other for this state
                % If no continue
                % If yes, store all values
                if (dist + past_metric) < cumulative_metrics(new_state, time + 1)
                    cumulative_metrics(new_state, time + 1) = dist + past_metric;
                    best_path_index(new_state, time + 1) = si;
                    best_path_input(new_state, time + 1) = inp;
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
function dec = soft_decoder(symbols, trellis, constellation)
end

function m = de_multiplex(rate, stream)
    stream = stream(:);
    N = length(stream) / rate;
    m = zeros(N, rate);
    for i = 1:rate
        m(:, i) = stream(i:rate:end);
    end    
end
