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

Description: Remove Inf values from a vector; return the largest 
contiguous interval with non-Inf values
%}
function [u, indx] = remInf(v)

infty = (v == Inf) + (v == -Inf);
if ~nnz(infty), u = v; indx = 1:numel(v); return; end
if nnz(infty) == numel(v); u = []; indx = []; return; end

firstInf	= find(infty, 1, 'first');
nonInf		= (v(firstInf:end) ~= Inf).*(v(firstInf:end) ~= -Inf);
endFirstInf = find(nonInf, 1, 'first') - 1 + firstInf - 1;
if ~numel(endFirstInf), endFirstInf = numel(v); end

if (firstInf - 1) >= (numel(v) - endFirstInf)
	u	= v(1:(firstInf-1));
	indx= 1:(firstInf - 1);
else
	[u, indx]	= remInf(v((endFirstInf + 1):end));
	indx		= endFirstInf + indx;
end