# MATLAB Binding for Cover Tree

## Description

This is a MATLAB binding for the Cover Tree implementation made by John Langford
<jl@hunch.net> and released in:

http://hunch.net/~jl/projects/cover_tree/cover_tree.html

The paper where cover tree data structure is described, the original implementation
and a FAQ about that code could be retrieved from the previous link.

## Installation

### Requirements

 * A C++ compatible compiler for MATLAB (only tested with g++ in Linux)
 * A MATLAB configured to use that C++ compiler to build MEX.
 * Boost.Serialization (>= 1.36)

### Build & set up the environment

 1. Donwload the repository to a directory and go there from within MATLAB.
 2. Edit [covertree_compile.m](covertree_compile.m) and check that Boost.Serialization shared library filename at [line 23](covertree_compile.m#L23) is ok.
 3. Run [covertree_compile](covertree_compile.m).

After the compiler successfully ends, we need to take some steps before using the toobox:

 1. Because `covertree_call` MEX depends on Boost.Serialization shared library, we must add its path
    (p. ej. '/usr/lib/') to LD_LIBRARY_PATH environment variable:
    
        setenv('LD_LIBRARY_PATH', [getenv('LD_LIBRARY_PATH') '/usr/lib/']);

 2. Probably you what to add the repository directory to MATLAB search path with:
 
        addpath('/path/to/repository/directory');

## Usage

Cover tree data structures could be handled through CoverTree class instances or
calling directly the operations supported by the covertree_call MEX.

### Using the CoverTree class

CoverTree methods:

<dl>
  <dt>CoverTree()</dt><dd>Create a new CoverTree object, inserting the points
specified in the first argument and setting the object as indicated by the
remaining options.
  <ul>
    <li><strong>Input arguments:</strong>
    <ol>
      <li>An single point or a cell array of points <strong>(optional)</strong>.
      <li>The string "DistanceFcn" followed by a distance function to be used
          <strong>(optional)</strong>.</li>
      <li>The string "PreSerializeFcn" followed by the function which
          process the points before to serialize the tree <strong>(optional)</strong>.</li>
      <li>The string "PostDeserializeFcn" followed by the function which
          process the points after to serialize the tree <strong>(optional)</strong>.</li>
    </ol></li>
    <li><strong>Output arguments:</strong><br>A new CoverTree class object.</li>
  </ul></dd>
  <dt>insert()</dt><dd>Insert new points in the tree.
  <ul>
    <li><strong>Input arguments:</strong><br>
    A single point or cell array of points.</li>
  </ul></dd>
  <dt>kNN()</dt><dd>Find the k nearest neighbors to specified points.
  <ul>
    <li><strong>Input arguments:</strong>
    <ol>
      <li>A CoverTree object containing the points to query.</li>
      <li>The number of k neighbors to find.</li>
    </ol></li>
    <li><strong>Output arguments:</strong></li>
    <ol>
      <li>A cell array containing an element for each point queried.
Each of these elements is another cell array with k+1 elements. The first
one is the point queried and the remaining are the k nearest neighbors.</li>
      <li>A cell array containing an element for each point queried.
Each of these elements is a vector of dimension k+1 with the distances
between the nearest neighbors and the point queried, placed in the same order
as in the first output argument.</li>
    </ol>
  </ul></dd>
  <dt>epsilonNN()</dt><dd>Search the ε-approximate nearest neighbors to the
specified points.
  <ul>
    <li><strong>Input arguments:</strong>
    <ol>
      <li>A CoverTree object containing the points to query.</li>
      <li> The value of ε, which defines the radius around the queried point
           where you want to find the nearest neighbors.</li>
    </ol></li>
    <li><strong>Output arguments:</strong>
    <ol>
      <li>A cell array containing an element for each point queried.
Each of these elements is another cell array. The first one is the point
queried and the remaining are the ε-approximate neighbors.</li>
      <li>A cell array containing an element for each point queried.
Each of these elements is a vector with the distances
between the nearest neighbors and the point queried, placed in the same order
as in the first output argument.</li>
    </ol></li>
  </ul></dd>
  <dt>unequalNN()</dt><dd>Find the nearest neighbor but not equal to the points
specified.
  <ul>
    <li><strong>Input arguments:</strong><br>
    A CoverTree object containing the points to query.</li>
    <li><strong>Output arguments:</strong>
    <ol>
      <li>A cell array containing an element for each point queried.
Each of these elements is another cell array with 2 elements. The first
one is the point queried and the other is the nearest neighbor.</li>
      <li>A cell array containing an element for each point queried.
Each of these elements is a vector with the distances between the nearest
neighbors and the point queried, placed in the same order as in the first
output argument.</li>
    </ol></li>
  </ul></dd>
  <dt>load()</dt><dd>Load a Cover Tree from the specified file. It must be invoked
on an initialized CoverTree object without points.
  <ul>
    <li><strong>Input arguments:</strong>
    <ol>
      <li>Name of the file from which to load the data structure.</li>
      <li>The string "text" if we want to retrieve a structure stored in text
format or "binary" (by default) if we want to retrieve a structure stored in binary
format <strong>(optional)</strong>.</li>
    </ol></li>
  </ul>
  <dt>save()</dt><dd>Save the Tree Object Tree Cover to a file.
  <ul>
    <li><strong>Input arguments:</strong>
    <ol>
      <li>Name of the file from which to load the data structure.</li>
      <li>The string "text" to save the structure in text format or "binary"
(default) to save in binary format <strong>(optional)</strong>.</li>
    </ol></li>
  </ul></dd>
  <dt>delete</dt><dd>Destroy the CoverTree object.
  </dd>
</dl>

CoverTree properties (all are read-only):

<dl>
  <dt>BreadthDistances</dt><dd>It is a vector with an element per node that
contains the number of children of every node.</dd>
  <dt>DepthDistances</dt><dd>It is a vector with an element per node that
contains the implicit level distance from each explicit node to the root.</dd>
  <dt>HeightDistances</dt><dd>It is a vector with an element per node that
contains the maximum number of explicit levels below each node.</dd>
</dl>

### Using covertree_call MEX

Most CoverTree class methods are implemented using `covertree_call`. The entry
point of this MEX always expects as first arguments a string and a unsigned
64-bit integer. The string specified the command to do, the integer is a
handler which identifies the Cover Tree affected by that command and the
remaining arguments must be valid parameters for the command.

Commands like `batch_create`, which create new trees, return a Cover Tree handler
that we can use with other commands. For them, the second argument of
`covertree_call` must be set to 0.

<dl>
  <dt>batch_create</dt><dd>Initialize a new Cover Tree and inserts the points
specified.
  <ul>
    <li><strong>Input arguments:</strong><br>
    A single point or cell array of points.</li>
    <li><strong>Output arguments:</strong><br>
    A handle that identifies the Cover Tree created.</li>
  </ul></dd>
  <dt>insert</dt><dd>Insert new points in the Cover Tree.
  <ul>
    <li><strong>Input arguments:</strong><br> A single point or cell array of points.</li>
  </ul></dd>
  <dt>k_nearest_neighbor</dt><dd>Find the k nearest neighbors to specified points.
  <ul>
    <li><strong>Input arguments:</strong>
    <ol>
      <li>A Cover Tree handler containing the points to query.</li>
      <li>The number of k neighbors to find.</li>
    </ol></li>
    <li><strong>Output arguments:</strong>
    <ol>
      <li>A cell array containing an element for each point queried.
Each of these elements is another cell array with k+1 elements. The first
one is the point queried and the remaining are the k nearest neighbors.</li>
      <li>A cell array containing an element for each point queried.
Each of these elements is a vector of dimension k+1 with the distances
between the nearest neighbors and the point queried, placed in the same order
as in the first output argument.</li>
    </ol></li>
  </ul></dd>
  <dt>epsilon_nearest_neighbor</dt><dd>Search the ε-approximate nearest neighbors to the
specified points.
  <ul>
    <li><strong>Input arguments:</strong>
    <ol>
      <li>A Cover Tree handler containing the points to query.</li>
      <li> The value of ε, which defines the radius around the queried point
           where you want to find the nearest neighbors.</li>
    </ol></li>
    <li><strong>Output arguments:</strong>
    <ol>
      <li>A cell array containing an element for each point queried.
Each of these elements is another cell array. The first one is the point
queried and the remaining are the ε-approximate neighbors.</li>
      <li>A cell array containing an element for each point queried.
Each of these elements is a vector with the distances
between the nearest neighbors and the point queried, placed in the same order
as in the first output argument.</li>
    </ol></li>
  </ul></dd>
  <dt>unequal_nearest_neighbor</dt><dd>Find the nearest neighbor but not equal to the points
specified.
  <ul>
    <li><strong>Input arguments:</strong><br>
    A Cover Tree handler containing the points to query.</li>
    <li><strong>Output arguments:</strong>
    <ol>
      <li>A cell array containing an element for each point queried.
Each of these elements is another cell array with 2 elements. The first
one is the point queried and the other is the nearest neighbor.</li>
      <li>A cell array containing an element for each point queried.
Each of these elements is a vector with the distances between the nearest
neighbors and the point queried, placed in the same order as in the first
output argument.</li>
    </ol></li>
  </ul></dd>
  <dt>breadth_dist</dt><dd>Returns the number of children of each node of Cover Tree.
  <ul>
    <li><strong>Output arguments:</strong><br>
    It is a vector with an element per node that
contains the number of children of every node.</li>
  </ul></dd>
  <dt>depth_dist</dt><dd>Returns the number of implicit levels between each node and
the root of Cover Tree.
  <ul>
    <li><strong>Output arguments:</strong><br>
    It is a vector with an element per node that contains the implicit level
distance from each explicit node to the root.</li>
  </ul></dd>
  <dt>height_dist</dt><dd>Returns the maximum number of explicit levels below
each Cover Tree node.
  <ul>
    <li><strong>Output arguments:</strong><br>
    It is a vector with an element per node that
contains the maximum number of explicit levels below each node.</li>
  </ul></dd>
  <dt>load</dt><dd>Load a Cover Tree from the specified file.
  <ul>
    <li><strong>Input arguments:</strong>
    <ol>
      <li>Name of the file from which to load the data structure.</li>
      <li>The string "text" if we want to retrieve a structure stored in text
format or "binary" if we want to retrieve a structure stored in binary
format (default).</li>
    </ol></li>
  </ul></dd>
  <dt>save</dt><dd>Save the Cover Tree to a file.
  <ul>
    <li><strong>Input arguments:</strong>
    <ol>
      <li>Name of the file from which to load the data structure.</li>
      <li> The string "text" to save the structure in text format or "binary" (by default)
           to save in binary format.</li>
    </ol></li>
  </ul></dd>
  <dt>delete</dt><dd>Frees the memory and other resources reserved by the
specified Cover Tree.</dd>
</dl>

## License

The code is dual licensed under the GPL and LGPL. Use whichever you prefer.


-- Jesús Torres <jmtorres@ull.es>
