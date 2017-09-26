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

Description: MR decomposition that forms a window around the current
%}
function [Cf, NzrData] = mrDecV05(Cforig, Sz, jmax, p, windw)
% Cforig: original wavelet coefficients
% jmin	: coarsest cell is of dimension 2^(-jmin)
% jmax	: finest cell is of dimension 2^(-jmax)
% p		: position of agent in absolute units
% windw	: window function

%% Initialization
N		= log2(Sz(end, 1)/Sz(1, 1));										% Level of decomposition N

ANzr	= [];
pCells	= [];

%% Add all approximation coefficients to the list
nzCf(1:Sz(1, 1)^2, :) = [ones(Sz(1,1)^2, 1) (1:Sz(1, 1)^2)' ...				% List of all non-zero coefficients (1, loc, value)
	(Cforig(1, 1:Sz(1, 1)^2))'];

%% Add detail coefficients within the window
for j = (-N):(jmax-1)														% From coarse to fine
	n = j + N + 1;
	posn	= floor((2^j)*p);
	pCells	= cat(1, pCells, [j posn']);
% 	fprintf('j = %i, n = %i, posn = (%i, %i)\n', j, n, posn(1), posn(2));
	for l = (posn(2) - windw(n)):(posn(2) + windw(n))
		for k = (posn(1) - windw(n)):(posn(1) + windw(n))
			ptr = Sz(1, 1)^2;												% End of app coeffs
			for m = 1:(n-1)
				ptr = ptr + 3*Sz(m+1,1)^2;
			end
			if (k>=0) && (k<Sz(n+1,1)) && (l>=0) && (l<Sz(n+1,1))			% if k and l within limits
				ANzr	= cat(1, ANzr, [j k l]);
				ptr		= ptr + k*Sz(n+1,1) + l + 1;
				nzCf	= cat(1, nzCf, [1 ptr Cforig(ptr); ...				% horizontal
					1 ptr+Sz(n+1,1)^2 Cforig(ptr+Sz(n+1,1)^2);...			% vertical
					1 ptr+2*Sz(n+1,1)^2 Cforig(ptr+2*Sz(n+1,1)^2)]);		% diagonal
			end
		end
	end
end
Cf				= sparse(nzCf(:,1), nzCf(:,2), nzCf(:,3), 1, Sz(end,1)^2);
NzrData.pCells	= [pCells; jmax floor((2^jmax)*p')];
NzrData.ANzr	= ANzr;
