%{
Copyright (c) 2017 Raghvendra V. Cowlagi. All rights reserved.

Copyright notice: 
=================
No part of this work may be reproduced without the written permission of
the copyright holder, except for non-profit and educational purposes under
the provisions of Title 17, USC Section 107 of the United States Copyright
Act of 1976. Reproduction of this work for commercial use is a violation of
copyright.


Disclaimer:
===========
This software program is intended for educational and research purposes.
The authors and the institutions with which the authors are affiliated are not
liable for damages resulting the application of this program, or any
section thereof, which may be attributed to any errors that may exist in
this program.


Author information:
===================
1. Raghvendra V. Cowlagi, Ph.D,
Assistant Professor,
Aerospace Engineering Program,
Worcester Polytechnic Institute, Worcester, MA.
Email: rvcowlagi@wpi.edu
Website: http://www.wpi.edu/~rvcowlagi

2. Panagiotis Tsiotras
College of Engineering Dean's Professor,
School of Aerospace Engineering,
Georgia Institute of Technology, Atlanta, GA.
Email: tsiotras@gatech.edu
Website: http://soliton.ae.gatech.edu/labs/ptsiotra/

Description: Transition matrix from wavelet coefficients
%}

% *************************************************************************
% Notation and identifier list
% ----------------------------
% N		: level of decomposition
% G		: adjacency matrix, square matrix with no. of rows/columns equal to
%		  no. of nodes (blocks) in the graph
% V		: cell locations (top left vertex), dimensions, and elevations
% Cf	: *original* (unapproximated) coefficient array
% Sz	: original coefficient size array
% inOldnotNew	: set of coefficients in AOld but not in ANew
% inNewnotOld	: set of coefficients in ANew but not in AOld
% AOld			: set of old non-zero coefficients
% VOld			: set of old cells
%**************************************************************************
function [G, V] = getNewGraphV02(AOld, GOld, VOld, inOldnotNew, inNewnotOld, Cf, Sz)
%% Initialize
N		= log2(Sz(end, 1)/Sz(1, 1));										% Level of decomposition N
NApp	= Sz(1, 1)^2;														% Total number of approximation coefficients
ElevM	= [1 1 1 1; 1 -1 1 -1; 1 1 -1 -1; 1 -1 -1 1];

addList = [];																% (elev, j, k, l) of cells to be added
remList = [];																% (j, k, l) of cells to be removed

%% Apply rules 2-4 from CDC paper
ANzr = inNewnotOld;
for q = 1:size(ANzr,1)
	j = ANzr(q, 1);		k = ANzr(q,2);		l = ANzr(q,3);
	n = j + N + 1;
	
	kApp	= floor( k*(2^(-N-j)) );
	lApp	= floor( l*(2^(-N-j)) );
	
	cf_app	= Cf(kApp*Sz(2,1) + lApp + 1);
	
	ptr = NApp;
	for m = 1:(n-1)
		ptr = ptr + 3*Sz(m+1,1)^2;
	end	
	cf_hor	= Cf(ptr + k*Sz(n+1,1) + l + 1);
	cf_ver	= Cf(ptr + Sz((n+1), 1)^2 + k*Sz(n+1,1) + l + 1);
	cf_dgn	= Cf(ptr + 2*Sz((n+1), 1)^2 + k*Sz(n+1,1) + l + 1);
	
	near_nzr= -N-1;
	if j == -N
		temp	= (2^j)*ElevM*[cf_app cf_hor cf_ver cf_dgn]';
		crsElev	= cf_app*(2^(-N));
	else
		% Find next coarser non-zero coefficient in the "same area" for
		% calculating the elevation		
		for jHat = (j-1):-1:(-N)
			kHat= floor( k*(2^(jHat-j)) );
			lHat= floor( l*(2^(jHat-j)) );
			
			if ismember([jHat kHat lHat], AOld, 'rows')						% NOTE: ismember of AOld, not ANzr
				near_nzr= jHat;
				kCap	= floor( k*(2^(jHat-j+1)) );
				lCap	= floor( l*(2^(jHat-j+1)) );
				
				crsCell = [[kCap lCap]*(2^(-jHat-1)) 2^(-jHat-1)];
				[~,cRec]= ismember(crsCell, VOld(:, 1:3), 'rows');			% NOTE: ismember of VOld, not addList
				
% % 				if ~cRec
% 					[tmp1, tmp2] = ismember(crsCell(1:2), VOld(:,1:2), 'rows');
% % 					if tmp1
% 						[j k l]
% 						[jHat kHat lHat]
% 						[kCap lCap]
% 						crsCell
% 						VOld(tmp2, :)
% % 					end
% % 				end
				
				temp	= (2^j)*ElevM*[VOld(cRec,4)*(2^(-j)) ...
					cf_hor cf_ver cf_dgn]';
				crsElev	= VOld(cRec,4);
				break;
			end
			
			if ismember([jHat kHat lHat], ANzr, 'rows')						% Also need the "standard" thing
				near_nzr= jHat;
				kCap	= floor( k*(2^(jHat-j+1)) );
				lCap	= floor( l*(2^(jHat-j+1)) );
				
				[~,cRec]= ismember([jHat+1 kCap lCap], addList(:, 2:4), 'rows');
				temp	= (2^j)*ElevM*[addList(cRec,1)*(2^(-j)) ...
					cf_hor cf_ver cf_dgn]';
				crsElev	= addList(cRec,1);
				break;
			end			
		end
		if (near_nzr == -N-1)												% No coarser coefficient non-zero in this area
			fprintf('WTF??\n')
			temp	= (2^j)*ElevM*[cf_app*(2^(-j-N)) cf_hor cf_ver cf_dgn]';
			crsElev = cf_app*(2^(-N));
		end
	end

	addList	= cat(1, addList, [temp (j+1)*ones(4,1)...						% Rule 2
		 [2*k 2*k 2*k+1 2*k+1]' [2*l 2*l+1 2*l 2*l+1]']);
	
	for jHat = (near_nzr+1):(j-1)											% Rule 3
		kHat = floor( k*(2^(jHat-j)) );
		lHat = floor( l*(2^(jHat-j)) );		
		
		addList	= cat(1, addList, [crsElev*ones(4,1) (jHat+1)*ones(4,1) ...
			[2*kHat 2*kHat 2*kHat+1 2*kHat+1]' [2*lHat 2*lHat+1 2*lHat 2*lHat+1]']);	
	end	
	
	for jHat = (-N):j														% Rule 4
		kHat = floor( k*(2^(jHat-j)) );
		lHat = floor( l*(2^(jHat-j)) );
		remList	= cat(1, remList, [jHat kHat lHat]);
	end
end

%% Apply "inverses" of rules 2-4 from CDC paper
% Assumes that for each non-zero coefficient at level j, "corresponding"
% coefficient at level j-1 (next coarser) is also non-zero. This will
% happen in the "window" approximation.
ANzr = inOldnotNew;
for q = 1:size(ANzr,1)
	j = ANzr(q, 1);		k = ANzr(q,2);		l = ANzr(q,3);
	remList = cat(1, remList, [(j+1)*ones(4,1)...							% "Inverse" of Rule 2
		 [2*k 2*k 2*k+1 2*k+1]' [2*l 2*l+1 2*l 2*l+1]']);
	 	
	crsElev = 0;
	for mk = 0:1
		for ml = 0:1
			crsCell = [[2*k+mk 2*l+ml]*2^(-j-1) 2^(-j-1)];
			[~,cRec]= ismember(crsCell, VOld(:,1:3), 'rows');
			crsElev	= crsElev + VOld(cRec,4);
		end
	end
	addList = cat(1, addList, [crsElev/4 j k l]);							% "Inverse" of Rule 4
end

%% Manipulate records to get cells
addList = unique(addList, 'rows');											% Newly added
remList = unique(remList, 'rows');											% Newly (to be) removed

if numel(addList) 
	[~,remIndx]= ismember(remList, addList(:,2:4), 'rows');
	remIndx	= remIndx(remIndx > 0);

	addList(remIndx, :)= [];													% Update newly added list
	inNewnotOldV		= [];													% Cells, not coefficients
	for n = 1:size(addList, 1)
		cellSize	= 2^(-addList(n,2));
		inNewnotOldV= cat(1, inNewnotOldV, [addList(n,3:4)*cellSize cellSize addList(n,1)]);
	end
else
	inNewnotOldV = [];
end

remCells = [];
for n = 1:size(remList, 1)
	cellSize = 2^(-remList(n,1));
	remCells = cat(1, remCells, [remList(n,2:3)*cellSize cellSize]);
end
[~,remIndx]= ismember(remCells, VOld(:,1:3), 'rows');
remIndx	= remIndx(remIndx > 0);
inOldnotNewV= VOld(remIndx, :);

[G, V]	= findAdj03(inNewnotOldV, inOldnotNewV, remIndx, VOld, GOld);
