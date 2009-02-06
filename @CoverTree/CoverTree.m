classdef CoverTree < handle
%COVERTREE Cover Tree data structure.
%   See http://hunch.net/~jl/projects/cover_tree/cover_tree.html

   properties (SetAccess = private, GetAccess = private)
       CoverTreeHandle = []                       % Cover Tree root node

       % Callbacks. Called from covertree_call()
       DistanceFcn = @CoverTree.defaultDistance   % Distance between points
       PreSerializeFcn = []                       % Called before serialize
       PostDeserializeFcn = []                    % Called after deserialize
   end

   properties (Dependent)
       % Properties calculated from data structure
       DepthDistances
       HeightDistances
       BreadthDistances
   end

   methods
       function CT = CoverTree(varargin)
           if nargin == 0
               return
           end

           iarg = 1;

           if iscell(varargin{iarg})
               P = varargin{iarg};
               iarg = iarg + 1;
           else
               P =[];
           end

           CT.parse_options(varargin{iarg:end})

           if ~ isempty(P)
               CoverTree.covertree_call('batch_create', CT, P);
           end
       end

       function insert(CT, P)
           if nargin ~= 2
               error ('Too few input arguments.');
           end

           if isempty(CT.CoverTreeHandle)
               CoverTree.covertree_call('batch_create', CT, P);
           else
               CoverTree.covertree_call('insert', CT, P);
           end
       end

       function [A D] = kNN(CT, queryCT, k)
           if nargin ~= 3
               error ('Too few input arguments.');
           end

           CT.check_covertree();

           if ~ isa(queryCT, 'CoverTree')
               error ('Argument #1 must be a CoverTree object.');
           end

           if ~ (isscalar(k) && isnumeric(k) && isreal(k))
               error ('Argument #2 must be a positive integer scalar.');
           end

           % Convert k to the nearest positive integer
           new_k = max(round(k), 1);
           if new_k ~= k
               warning('CoverTree:Rounding', ...
                   'k rounded to the nearest positive integer: %d', new_k);
           end

           [A D] = CoverTree.covertree_call( ...
               'k_nearest_neighbor', CT, queryCT, new_k);
       end

       function [A D] = epsilonNN(CT, queryCT, epsilon)
           if nargin ~= 3
               error ('Too few input arguments.');
           end

           CT.check_covertree();

           if ~ isa(queryCT, 'CoverTree')
               error ('Argument #1 must be a CoverTree object.');
           end

           if ~ (isscalar(epsilon) && isnumeric(epsilon) && isreal(epsilon))
               error ('Argument #2 must be a positive real scalar.');
           end

           % Convert k to the nearest positive single
           new_epsilon = max(single(epsilon), eps(single(0)));
           if new_epsilon ~= epsilon
               warning('CoverTree:Rounding', ...
                   'epsilon rounded to the nearest positive single: %d', ...
                   new_epsilon);
           end

           [A D] = CoverTree.covertree_call('epsilon_nearest_neighbor', ...
               CT, queryCT, new_epsilon);
       end

       function [A D] = unequalNN(CT, queryCT)
           if nargin ~= 2
               error ('Too few input arguments.');
           end

           CT.check_covertree();

           if ~ isa(queryCT, 'CoverTree')
               error ('Argument #1 must be a CoverTree object.');
           end

           [A D] = CoverTree.covertree_call('unequal_nearest_neighbor', ...
               CT, queryCT);
       end

       function load(CT, filename, mode)
           if nargin < 2
                error ('Too few input arguments.');
           end

           if ~ isempty(CT.CoverTreeHandle)
               error('Cover Tree object must be empty to load.');
           end

           if ~ (ischar(filename) && isvector(filename))
               error('Argument #1 must be a string.');
           end

           if nargin < 3
               mode = 'binary';
           elseif ~ (ischar(mode) && isvector(mode))
               error('Argument #2 must be a string.');
           end

           CoverTree.covertree_call('load', CT, filename, mode);
       end

       function save(CT, filename, mode)
           if nargin < 2
                error ('Too few input arguments.');
           end

           CT.check_covertree();

           if ~ (ischar(filename) && isvector(filename))
               error('Argument #1 must be a string.');
           end

           if nargin < 3
               mode = 'binary';
           elseif ~ (ischar(mode) && isvector(mode))
               error('Argument #2 must be a string.');
           end

           CoverTree.covertree_call('save', CT, filename, mode);
       end

       function delete(CT)
           if ~ isempty(CT.CoverTreeHandle)
               CoverTree.covertree_call('delete', CT);
           end
           CT.CoverTreeHandle = [];
       end

       % Getters for distance properties.
       % Must be calculated from data structure

       function D = get.DepthDistances(CT)
           D = [];
           if ~isempty(CT.CoverTreeHandle)
               D = CoverTree.covertree_call('depth_dist', CT);
           end
       end

       function [H M] = get.HeightDistances(CT)
           H = []; M = [];
           if ~isempty(CT.CoverTreeHandle)
               [H M] = CoverTree.covertree_call('height_dist', CT);
           end
       end

       function B = get.BreadthDistances(CT)
           B = [];
           if isempty(CT.CoverTreeHandle)
               B = CoverTree.covertree_call('breadth_dist', CT);
           end

       end

       % Methods to support the typical behavior

       function disp(CT)
           disp(inputname(1));
       end

       function display(CT)
           disp(' ');
           disp(inputname(1));
           disp(' ');
       end
   end

   methods (Access = protected)
       function check_covertree(CT)
           if isempty(CT.CoverTreeHandle)
               error('Cover Tree object is empty.');
           end
       end

       function parse_options(CT, varargin)
           iarg = 1;

           while iarg <= length(varargin)
               option = varargin{iarg};
               switch option
                   case {'DistanceFcn', 'PreSerializeFcn', 'PostDeserializeFcn'}
                       iarg = iarg + 1;
                       value = varargin{iarg};
                       if ~ (isa(value, 'function_handle') || ischar(value))
                           error ('Invalid callback function for property %s.', option);
                       end
                       CT.(option) = value;
                   otherwise
                       error ('Invalid argument #%u.', iarg);
               end

               iarg = iarg + 1;
           end
       end
   end

   methods (Static)
       s = makeDistanceFcnString(fname, shlib);
   end

   methods (Static, Access = protected)
       [LHS1 LHS2] = covertree_call(fname, CT, RHS1, RHS2);

       function d = defaultDistance(P1, P2, ub)
           e = P1(:) - P2(:);
           d = sqrt(e' * e);
       end
   end
end
