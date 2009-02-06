function s = makeDistanceFcnString(fname, shlib)
%MAKEDISTANCEFCNSTRING Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 1
        error ('Too few input arguments.');
    end

    fnamesz = size(fname);
    if ~ ischar(fname) || fnamesz(1) ~= 1 || fnamesz(2) == 0
        error('Argument #1 must be a function name string.');
    end

    if nargin < 2
        shlib = ['distances.' mexext];
    end

    shlibsz = size(shlib);
    if ~ ischar(shlib) || shlibsz(1) ~= 1 || shlibsz(2) == 0
        error('Argument #2 must be a function name string.');
    end

    if shlib(1) ~= filesep
        pkg = what('CoverTree');
        shlib = [pkg(1).path filesep shlib];
    end

    s = [shlib ':' fname];
end
