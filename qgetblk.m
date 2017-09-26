function [val,selected] = qgetblk(A,S,I,extent,border);
%QGETBLK Get block values in a N-D quadtree decomposition.
%   VALS = QGETBLK(X,S,[],EXTENT) returns the values of the blocks 
%   having the dimensions specified in the vector EXTENT in an N-D
%   array whose last dimension indexes the blocks.
%
%   [VALS,I] = QGETBLK(X,S,[],EXTENT) also returns a vector of indices
%   I such that S.locs(I(j),:) is the branch point whose values are
%   now in the jth layer of VALS. 
% 
%   [VALS,I] = QGETBLK(X,S,[],EXTENT,BORDER) returns the values of the
%   blocks including values surrounding the blocks that are
%   within an overlapping region around the blocks.  Thus, a block 
%   specified by EXTENT is increased by BORDER to a size of 
%   EXTENT+2.*BORDER.  QGETBLK pads the border with NAN's if necessary,
%   for blocks on the edges of X.  EXTENT and BORDER must both be 
%   vectors of length DIM where DIM is NDIMS(A);
%
%   QGETBLK(X,S,I,EXTENT,BORDER) is the same as above except that only 
%   the blocks with locations at S.locs(I,:) are considered.  Blocks
%   at those locations which do not match the given EXTENT are excluded.
%      
%   Example:
%      X = peaks(50);
%      S = qdecomp(X,1.5);
%      [B,I] = qgetblk(X,S,[],[1 2],[1 1])	
%

%Written by Daniel P. Dougherty (dpdoughe@stat.ncsu.edu)
%8:14PM 09/04/00
%Copyright 2000 Daniel P. Dougherty.  All Rights Reserved.


V = size(A);


if (nargin < 4)
	error('At least 4 inputs are required.');
elseif (nargin > 5)
	error('Too many input argumnents.');
elseif (nargin == 4)
	border = zeros(1,length(V));
elseif (nargin == 5)
	if ((any(border<0) | ~all(isreal(border))) | (any(mod(border,1) ~= 0)))
		error('Border must be non-negative integer.');
	end
end


if isempty(I)
	I = [1:size(S.locs,1)]';
else
	I = I(:);
end


if issparse(A)
	error('Sorry sparse matrices are not allowed.');
end
	
if (size(S.locs,2) == 1) %Column-ordered indexing is being used...

	if ~all(border == 0)
		%Need to add border into the calcuations.  To do this
		%simply create a new matrix padded with the appropriate
		%number of nan's,adjust the decomposition S to coincide
		%and then process as usual.

		B = repmat(nan,size(A)+2.*border);

		%Use some of the stuff we already have in this toolbox
		%making it easy to set the interior points in N-D array.	

		%Make a little decomposition struct on the fly!
		
		%Set the internal core to A.  This leaves a shell 
		%of nan's.
		
		T.size   = size(B);
		T.locs   = [border+1];
		T.extent = size(A);
		
		A = qsetblk(B,T,[],T.extent(1,:),[],A);
		extent = extent + 2.*border;
		
		
		%New Stuff	
		H = cell(1,length(S.size));
		[H{:}] = ind2sub(S.size,S.locs);
		
		%%%
		
		S.size = S.size + 2.*border;
		S.locs = sub2ind(S.size,H{:});
		
		%Make border ready for element-by-element addition
		%with S.extent
		
		border = repmat(border,size(S.locs,1),1);
		S.extent = S.extent + 2.*border;
		V = size(A); %Re-set size of A.
		
		%Now process as usual :) 
	end

	selected = I(find(all((S.extent(I,:) == repmat(extent,length(I),1))'))');%Get the row indices
	
									
	
	if all(extent == S.size) & (length(selected) == 1) %Then this is an easy one to do since it is the whole
			      %array.
		val = A;
		return;
	end

										%of S that correspond to
										%these blocks.
	if isempty(selected) %Nothing to do...

		val = [];
		return;
	else


		%NOTICE::  Any improvements on the efficiency of the
		%following code would be greatly appreciated!  
		
	 
		ind = S.locs(selected);%The locations of the branch points.
		
		for i = 1:length(extent) 
			g{i} = [1:extent(i)];
		end

		[g{:}] = ndgrid(g{:});


		I = reshape(sub2ind_spec(S.size,ones(1,length(extent)),extent(:)',g{:}),extent);


		[I,J] = meshgrid(ind-1,I);

		ind = I + J;

		val = A(ind);

		val = reshape(val,[extent length(selected)]);
		%Each layer (last dimension indexes layers)
		%contains data for a single block.



	end



else	%Row,column,.... indexing is being used.

	if ~all(border == 0)
		%Need to add border into the cacluations.  To do this
		%simply create a new matrix padded with the appropriate
		%number of nan's,adjust the decomposition S to coincide
		%and then process as usual.

		B = repmat(nan,size(A)+2.*border);

		%Use some of the stuff we already have in this toolbox
		%making it easy to set the interior points in N-D array.	

		%Make a little decomposition struct on the fly!
		T.locs   = [border+1];
		T.extent = size(A);
		T.size   = size(B);
		A = qsetblk(B,T,[],T.extent(1,:),[],A);
		extent = extent + 2.*border;

		S.size = S.size + 2.*border;
		%Make border ready for element-by-element addition.
		border = repmat(border,size(S.locs,1),1);
		S.extent = S.extent + 2.*border;
		V = size(A); %Re-set size of A.

		%Now process as usual :) 
	end



	selected = I(find(all((S.extent(I,:) == repmat(extent,size(S.extent,1),1))'))');%Get the row indices
										%of S that correspond to
	
										
	
	if all(extent == S.size) & (length(selected) == 1) %Then this is an easy one to do since it is the whole
			      %array.
			     
		val = A;
		return;
	end
	
										%these blocks.

	if isempty(selected) %Nothing to do...

		val = [];
		return;
	else


		%NOTICE::  Any improvements on the efficiency of the
		%following code would be greatly appreciated!  

		ind = S.locs(selected,:);%The locations of the branch points.

		ind = num2cell(ind,1);  %The first cell holds the coordinates
					%along the first dimension and so on...

		ind = sub2ind(S.size,ind{:}); %Convert the subscripts to equivalent
					       %column-orderd index.


		for i = 1:length(extent) 
			g{i} = [1:extent(i)];
		end

		[g{:}] = ndgrid(g{:});

		I = reshape(sub2ind_spec(S.size,ones(1,length(extent)),extent(:)',g{:}),extent);

		

		[I,J] = meshgrid(ind-1,I);

		ind = I + J;


		val = A(ind);

		val = reshape(val,[extent length(selected)]);
		%Each layer (last dimension indexes layers)
		%contains data for a single block.



	end


end

function ndx =  sub2ind_spec(siz,mn,mx,varargin)
%This is a special version of sub2ind that we can use when we know 
%the minimum and maximum values in each cell of the input and 
%all subscript arrays are the same size.
 
n = length(siz);
k = [1 cumprod(siz(1:end-1))];
ndx = 1;
for i = 1:n,
         ndx = ndx + (varargin{i}-1)*k(i);
end


