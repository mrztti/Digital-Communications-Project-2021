classdef LoadingBar
    properties
        elapsed
        target
    end
    methods
        function lb = LoadingBar(steps)
            lb.elapsed = 0;
            lb.target = steps;
            lb_print(lb.elapsed/lb.target);
        end
        function obj = step(obj, size)
            obj.elapsed = obj.elapsed + size;
            lb_print(obj.elapsed/obj.target)
        end
        function obj = set(obj, elapsed)
            obj.elapsed = elapsed;
            lb_print(obj.elapsed/obj.target)
        end
    end
end
function lb_print(perc)
    clc;
    bar_string = make_bar(perc);
    fprintf(strcat('%.2f %% \t',bar_string, '\n'), perc*100);
end

function str = make_bar(perc)
    str = '[◻◻◻◻◻◻◻◻◻◻◻◻◻◻◻◻◻◻◻◻]';
    p = 2;
    while (p-1)/20 <= perc
        str(p) = '◼';
        p = p+1;
    end
end

