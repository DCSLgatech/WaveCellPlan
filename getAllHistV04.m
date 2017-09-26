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

Description: Compute all H-length histories of cells from a given cell
%}

function allHist = getAllHistV04(G, actNode, H, seqPos, rawSeq1, allHist)
% getAllHist returns a list of all H-length traversals from a given node in
% a given graph. The program makes recursive calls to itself, traversing
% one step forward in each call. Each history returned is of length H and
% includes the given starting node.
%
% Description of input variables
% ------------------------------
% G			: adjacency matrix of cells
% actNode	: active node from which all sequences start
% H			: fixed length of sequences to be found
% rawSeq1	: input for propagating histories generated thus far to
%			  further recursive calls
% seqList	: list of fully formed histories known at all levels of recursion

actNhbrs	= find(G(:, actNode));											% Step forward once

for count1 = 1:size(actNhbrs, 1)											% For all neighbors
	newNode	= actNhbrs(count1);
	if seqPos == H - 2
		temp = [rawSeq1 newNode];
		if any(rawSeq1 == newNode), continue; end
		if any(G(newNode, rawSeq1(1:(end-1)))), continue; end
		allHist = cat(1, allHist, temp);		
	else
		rawSeq2	= cat(2, rawSeq1, newNode);
		if any(rawSeq1 == newNode), continue; end
		if any(G(newNode, rawSeq1(1:(end-1)))), continue; end
% 		if size(unique(rawSeq2), 2) < size(rawSeq2, 2)						% If history large enough
% 			continue; 
% 		else																% Recursive call
			allHist	= getAllHistV04(G, newNode, H, seqPos + 1, rawSeq2, allHist);
% 		end
	end
end