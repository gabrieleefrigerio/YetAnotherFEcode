function name = run_host()
%RUN_HOST Identifier of the machine a result was produced on.
%
%   name = RUN_HOST()
%
% Stored in every result file so that a cpu_time measured on one machine is
% never silently compared with one measured on another.
%
% The distinction matters because the two axes of the comparison do not behave
% the same way. ACCURACY is machine independent: same code, same tolerances,
% same numbers, so an error study can be split across whatever hardware is
% free. COST is not: a speedup assembled from a FOM timed on one host and a
% ROM timed on another is a ratio between two computers, not between two
% models, and nothing in the result file would reveal it.
%
% See also RUN_FOM, RUN_ROM_SWEEP.

    name = getenv('COMPUTERNAME');        % Windows
    if isempty(name)
        name = getenv('HOSTNAME');        % most Linux shells
    end
    if isempty(name)
        [status, out] = system('hostname');
        if status == 0, name = strtrim(out); end
    end
    if isempty(name)
        name = 'unknown';
    end
end
