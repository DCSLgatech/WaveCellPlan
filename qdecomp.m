function S = qdecomp(varargin)

%QDECOMP Perform N-D quadtree decomposition.
%
%   S = QDECOMP(X) performs a quadtree decomposition of the
%      array X, and returns the tree structure in the
%      struct S.  The rows of S.locs contain the branch points
%      (or vertices) and the rows S.extent contain the dimensions
%      of the cell at that branch point.
%
%   S = QDECOMP(X,THRESHOLD) splits a block if the maximum value
%      of the block elements minus the minimum value of the block
%      elements is greater than THRESHOLD. 
%
%   S = QDECOMP(X,THRESHOLD,MINDIMS,MAXDIMS) splits any cell
%      with all of its dimensions larger that MAXDIMS.  Any block
%      with any of its dimensions less than MINDIMS is not split.
%      P1,P2,... are any additional arguments which should be 
%      included in the calls to FUN.  If MINDIMS or MAXDIMS are 
%      left empty, the default values of ones(length(V),1) and V 
%      are used where V is the size of X.  The algorithm will try
%      to meet both criteria simultaneously, but this will not 
%      always be possible since quadtree is based on bisection.  
%
%
%   S = QDECOMP(X,FUN,MINDIMS,MAXDIMS,P1,P2,...) uses the user-
%      supplied function FUN to determine whether to split a block. 
%      QDECOMP calls FUN with all the current blocks of a given size.
%      All the blocks of a given size are stored in an ND array in 
%      which the last dimension indexes the blocks. FUN should return 
%      a K-element vector whose values are 1 if the Kth block should 
%      be split, and 0 otherwise. 
%
%      FUN can be a string containing the name of a function, a string
%      containing an expression, or an inline function.  FUN should 
%      have the following input structure 
%
%         SPLIT = FUNC(BLOCKS,S,P1,P2,...).  
%
%      Here BLOCKS represent the ND array of blocks, S is the current 
%      tree decomposition and P1, P2, ... are any additional parameters
%      specified by the user. 
%
%      Example:
%         P  = peaks(50);
%         S = qdecomp(P,0.8)
%         P2 = plotqdecomp(S);
%

%Written by Daniel P. Dougherty (dpdoughe@stat.ncsu.edu)
%8:14PM 09/04/00
%Copyright 2000 Daniel P. Dougherty.  All Rights Reserved.



global ANIMATE_QDECOMP;

[A, func, params, minDim, maxDim,threshold, msg] = ParseInputs(varargin{:});

%A is the data.
%func is the function which returns a 1 if a block should be split and
%	a 0 otherwise.
%params include threshold,minDim,maxDim and user-supplied parameters
%	which should be passed to func.
%minDim is the smallest block allowed.
%maxDim is the larges block allowed.


if (~isempty(msg))
    error(msg);
end

if (nargout == 1) 
	ANIMATE_QDECOMP = logical(0);
elseif (nargout == 0)
	ANIMATE_QDECOMP = logical(1);
end


S = Recit(A,func,minDim,maxDim,threshold,params);

%By default, produce the decomposition in column-ordered indexing.
%This will reduce the memory used in storage (especially as the
%dimensionality gets large).

V = num2cell(S.locs,1);
S.locs = sub2ind(S.size,V{:}); 

%DONE!!!


function S = Recit(A,func,minDim,maxDim,threshold,params);
%A is a matrix which is possibly going to get broken up.  A may be
%a submatrix of a larger matrix.

global ANIMATE_QDECOMP;


%Start to create a decomposition struct.
S = struct('locs',zeros(1,length(size(A))),'extent',zeros(1,length(size(A))));

Sold = S; %Sold is a snapshot of the decomposition which we will 
	%compare with the new decomposition in order to determine 
	%which blocks are new blocks.  A new block is one which has 
	%_either_ a different extent _or_ a different loc.  Just comparing
	%the extent's or loc's will not reveal all new blocks!!

S.locs(1,:) = ones(1,length(size(A)));
S.extent(1,:) = size(A);
S.size = size(A);

index = 1;
newindex = index;

while 1; %Keep looping until we are done...
%%%	
	if ANIMATE_QDECOMP
		if (newindex == 1) %First time through
			plotqdecomp(S);
			drawnow;
		else	
			plotqdecomp(S,'anim');
			drawnow;
		end
	end
	
	Sold = S;
	
	this_extent = unique(S.extent(newindex,:),'rows'); %Get the new list of blocksizes
						%(extents) that we need to
						%deal with.
	
	%Get indices in  S that should be split.
	index = [];
	
	for i = 1:size(this_extent,1) 
	
		[blocks,idx] = qgetblk(A,S,[],this_extent(i,:)); %Blocks is an ND array containing
						      %The actual values of a block
						      %in each "layer".  The layer dimension
						      %is always the last dimension in the
						      %ND array.
						      %
						      %idx is a vector of the row index
						      %into S of the blocks selected.
		
		%See if this is 1D data...
		if ((length(S.size) == 2) & any(S.size == 1)) %Then this 
							%is 1D data which needs to 
							%be treated a little differently
							%than higher dimension data.
				if (max(this_extent(i,:))./2 < max(minDim))
					
					doSplit = logical(zeros(size(idx,1),1));
				else
					if (max(this_extent(i,:)) > max(maxDim))
						doSplit = logical(ones(size(idx,1),1));
					else

						if ~isempty(params)
							doSplit = feval(func,blocks,S,params); %A vector of zeros and
								%ones.  Each "layer" gets 
								%a one if it should be split.
								%Dosplit is the same length as idx!
						else
							doSplit = feval(func,blocks,S,threshold);
						end

						%Make sure the user's function is returning the 
						%correct kind of output.
						if ~all(islogical(doSplit))
								error('Function must return logical array.');
						end
				        end
						
				end	
			
				
		else %Higher dimension data.
			
			if (any(this_extent(i,:)./2 < minDim))
					
			 		doSplit = logical(zeros(size(idx,1),1));

			else
				if (any(this_extent(i,:) > maxDim))
					doSplit = logical(ones(size(idx,1),1));
				else

					if ~isempty(params)
						doSplit = feval(func,blocks,S,params); %A vector of zeros and
							%ones.  Each "layer" gets 
							%a one if it should be split.
							%Dosplit is the same length as idx!
					else
						doSplit = feval(func,blocks,S,threshold);
					end

					%Make sure the user's function is returning the 
					%correct kind of output.
					if ~all(islogical(doSplit))
							error('Function must return logical array.');
					end
				end
			end
		end
		
		
		idx = idx(doSplit); %Get only the indices into S of the
				% blocks that should be split.
		
		
		index = [index;idx]; %Gather-up the indices at each turn of the
				%loop (i.e. for each extent size).
	end
	
	if isempty(index) %We are done!
		break; %Break out of the while loop.
	end

	 
	if ((length(S.size) == 2) & any(S.size == 1)) %This is 1D data
	
		a = ceil(S.extent(index,2)./2); %Extents of new primary cells.
		new_locs = [S.locs(index,2) S.locs(index,2)+a];
		new_locs = [ones(length(new_locs(:)),1) new_locs(:)];
		new_extents = [a S.extent(index,2)-a];
		new_extents = [ones(length(new_extents(:)),1) new_extents(:)];
		 
	else	%Higher-dimensional data. 
	
		extent_old = S.extent(index,:); %Matrix of extents of the cells that are
						%to be split.
		extent_new = ceil(S.extent(index,:)./2); %Just the extent of the "primary"
							%cell". Now we need to find the
							%extents of the other cells. 

		extent_calc = extent_old - extent_new;
		vertex_calc = extent_new + S.locs(index,:);
		
		new_locs = hypercube(extent_new,S.locs(index,:));
		new_extents = hypercube(extent_calc,vertex_calc)-new_locs; 
	end
	

	newindex = [index(:);[size(S.locs,1)+1:size(S.locs,1)+size(new_locs,1)-size(index,1)]'];
	%Replace old branch points with new values.
	S.locs(index,:)   = new_locs(1:size(index,1),:);
	S.extent(index,:) = new_extents(1:size(index,1),:);
	%And append the remaining new values to the bottom of the list.
	S.locs(end+1:end+size(new_locs,1)-size(index,1),:) = new_locs(size(index,1)+1:end,:); 
	S.extent(end+1:end+size(new_extents,1)-size(index,1),:) = new_extents(size(index,1)+1:end,:); 
	 
	
end



%%%
%%% Subfunction RECDECOMP_Split
%%% (The default split-tester)

function doSplit = QDECOMP_Split(blocks,S,threshold)

if (threshold == 0)
	doSplit = logical(ones(1,size(blocks,length(S.size)+1)));
	
else

	if ((length(S.size) == 2) & any(S.size == 1))%This is 1D data.
		maxVal = max(blocks);
		minVal = min(blocks);
	else %Higher dimension data.
		maxVal = max(blocks);
		minVal = min(blocks);
		for i = 2:(length(S.size)) 
			maxVal = max(maxVal);   %Vector the same length as the number of blocks.
			minVal = min(minVal); %Vector the same length as the number of blocks.
		end
	end


	doSplit = (double(maxVal(:)) - double(minVal(:))) > threshold;
end		

 
%%%
%%% Subfunction ParseInputs
%%%
function [A, func, params, minDim, maxDim, threshold, msg] = ParseInputs(varargin)




if (nargin == 0)
    msg = 'Too few input arguments';
    return;
end
%%


%Defaults
func = 'QDECOMP_Split';
params = [];
threshold = 0;
msg = '';
minDim = ones(1,length(size(varargin{1})));  %Allows a block to be as small as one pixel.
maxDim = size(varargin{1});  %Allows a block to be as large as the original
				%block.



switch length(varargin)
case 1
	A = varargin{1};
case 2
	A = varargin{1};
	
	if (ischar(varargin{2}) | isa(varargin{2},'inline'))
		
		func = fcnchk(varargin{2},0);
		
	elseif (isnumeric(varargin{2}))
		threshold = varargin{2};
		if isempty(threshold)
			threshold = 0;
		end	
	else
		error('Second input must be a function or a number.');
	end

case 3
	A = varargin{1};
	
	if (ischar(varargin{2}) | isa(varargin{2},'inline'))
		
		func = fcnchk(varargin{2},0);
		
	elseif (isnumeric(varargin{2}))
		threshold = varargin{2};
		if isempty(threshold)
			threshold = 0;
		end		
	else
		error('Second input must be a function or a number.');
	end

	if ~isempty(varargin{3})
		if isnumeric(varargin{3})
			minDim = varargin{3};
		else
			error('Third input can be either empty or a vector.');
		end
	end	
		

case 4

	A = varargin{1};
	
	if (ischar(varargin{2}) | isa(varargin{2},'inline'))
		
		func = fcnchk(varargin{2},0);
		
	elseif (isnumeric(varargin{2}))
		threshold = varargin{2};
		if isempty(threshold)
			threshold = 0;
		end	
	else
		error('Second input must be a function or a number.');
	end

	if ~isempty(varargin{3})
		if isnumeric(varargin{3})
			minDim = varargin{3};
		else
			error('Third input can be either empty or a vector.');
		end
	end

	if ~isempty(varargin{4})
		if isnumeric(varargin{4})
			maxDim = varargin{4};
		else
			error('Fourth input can be either empty or a vector.');
		end
	end

otherwise
	A = varargin{1};
	
	if (ischar(varargin{2}) | isa(varargin{2},'inline'))
		
		func = varargin{2};
		
	elseif (isnumeric(varargin{2}))
		threshold = varargin{2};
		if isempty(threshold)
			threshold = 0;
		end	
	else
		error('Second input must be a function or a number.');
	end

	if ~isempty(varargin{3})
		if isnumeric(varargin{3})
			minDim = varargin{3};
		else
			error('Third input can be either empty or a vector.');
		end
	end

	if ~isempty(varargin{4})
		if isnumeric(varargin{4})
			maxDim = varargin{4};
		else
			error('Fourth input can be either empty or a vector.');
		end
	end
	
	params = varargin{4:end};
	
	fcnchk(func,length(params))
end

%Special for 1D data.

if ((length(size(A)) == 2) & any(size(A) == 1)) %This is 1D data.
	A = A(:)';
	maxDim = [1,max(maxDim)];
end


if any(maxDim < minDim)
	error('Max. block dimensions must be larger than min. block dimensions.');
end


if ~isnumeric(A)
	error('Input array must be numeric.');
end
