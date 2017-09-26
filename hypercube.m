function VERTICES = hypercube(X,vertex);
%HYPERCUBE Return the vertices of a hypercube.
%   HYPERCUBE Return the vertices of a hypercube relative to a given 
%   vertex (or origin).  
%
%   V = HYPERCUBE(X,v)  The lengths along each dimension are given in the
%   vector X.  The starting vertex (default is zeros(1,DIM)) is specified in v.  
%
%   If v is a matrix it must have DIM columns and it is assumed each row corresponds
%   to a separate vertex.  If X and v are matrices, they must be the same size.
%
%
%   Example:
%      verts = hypercube([1 2 3]);
%

%Written by Daniel P. Dougherty (dpdoughe@stat.ncsu.edu)
% 9:28PM 08/31/00
%Copyright 2000 Daniel P. Dougherty.  All Rights Reserved.


%Create logical indices to index into the X (which will have to be repmated
%along the way).

DIM = size(X,2);

if (nargin == 1)
	vertex = zeros(1,DIM);
end

if ((size(X,1) >1 ) & (size(vertex,1) > 1))
	if ~all(size(X) == size(vertex))
		error('Vertex and length matrices must be same size.');
	end
end



num_cubes = max(size(X,1),size(vertex,1));


%Pre-allocate for the output.
VERTICES = reshape(repmat(X,1,(2^DIM)*(num_cubes/size(X,1)))',DIM,(2^DIM)*num_cubes)';
vertex   = reshape(repmat(vertex,1,(2^DIM)*(num_cubes/size(vertex,1)))',DIM,(2^DIM)*num_cubes)';

%%%
x = {[0 1]};
x = repmat(x,1,DIM);

[x{:}] = ndgrid(x{:});
x = fliplr(x);
x = cat(DIM,x{:});
x = logical(reshape(x,2^DIM,DIM));
x = repmat(x,num_cubes,1);

%Done creating the logical indices.
%%%%
VERTICES(~x) = 0;
VERTICES = VERTICES + vertex;




