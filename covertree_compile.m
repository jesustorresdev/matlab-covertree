function covertree_compile(target)
%COVERTREE_COMPILE Summary of this function goes here
%   Detailed explanation goes here

    MEXOPTS = [];

    if nargin > 0
        switch target
            case 'debug'
                MEXOPTS = '-g';
        end
    end

    pkg = what('CoverTree');
    origdir = cd([pkg(1).path]);
    
    try
        % Compile covertree_call() method of CoverTree class
        %mex covertree_call.cc -inline -Iinclude -ldl -lboost_serialization
        cmd = ['mex ' MEXOPTS ' covertree_call.cc ' ...
            '-inline -Iinclude -ldl -lboost_serialization-gcc41-1_36'];
        disp(cmd);
        system(cmd);

        % Compile auxiliary library of distances
        % % Supplement variables already defined (e.g. LDFLAGS) only works if
        % % mex is executed from shell.
        cmd = ['mex ' MEXOPTS ' distances.cc -lmwlapack -lmwblas ' ...
            'LDFLAGS=''$LDFLAGS \"-Wl,--version-script,distances.map\"'''];
        disp(cmd);
        system(cmd);
    catch ME
        cd(origdir);
        rethrow(ME)
    end

    cd(origdir);
end
