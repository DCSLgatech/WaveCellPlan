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

Description: 
%}
function [NzrDataNew, inOldnotNew, inNewnotOld] = ...
	getNewNzrV02(NzrData, Sz, jmax, dlta, windw)

%% Initialization
N			= log2(Sz(end, 1)/Sz(1, 1));
rsFine		= 2^(-jmax-1);
% pCellsNew	= [((-N):jmax)' zeros(jmax+N+1, 2)];
ANzr		= NzrData.ANzr;
pCells		= NzrData.pCells;
pCellsNew	= pCells;

inOldnotNew = [];
inNewnotOld = [];
rsdu		= zeros(2,1);
incP		= zeros(2,1);

%%
for j = (-N):(jmax-1)
	n = j + N + 1;
	for m = 1:2
			
		rsdu(m,1)	= (2^(j-jmax))*(pCells(end,m+1) + rsFine) - pCells(n,m+1);
		incP(m,1)	= (2^(j-jmax))*dlta(m) + rsdu(m);
		
		if ((0 <= incP(m)) && (incP(m) < 1))
			pCellsNew(n, m+1) = pCells(n, m+1);
			continue;
		elseif (incP(m) >= 1)
			pCellsNew(n, m+1) = pCellsNew(n, m+1) + 1;
			if m == 1
				if (pCells(n, 2) - windw(n) >= 0)
					inOldnotNew = cat(1, inOldnotNew, ...
						[j*ones(2*windw(n)+1, 1) ...
						(pCells(n,2)-windw(n))*ones(2*windw(n)+1, 1) ...
						((pCells(n,3)-windw(n)):(pCells(n,3)+windw(n)))']);
				end
				if (pCellsNew(n, 2) + windw(n) < Sz(n+1,1))
					inNewnotOld = cat(1, inNewnotOld, ...
						[j*ones(2*windw(n)+1, 1) ...
						(pCellsNew(n,2)+windw(n))*ones(2*windw(n)+1, 1) ...
						((pCellsNew(n,3)-windw(n)):(pCellsNew(n,3)+windw(n)))']);
				end
			else
				if (pCells(n, 3) - windw(n) >= 0)
					inOldnotNew = cat(1, inOldnotNew, ...
						[j*ones(2*windw(n)+1, 1) ...
						((pCells(n,2)-windw(n)):(pCells(n,2)+windw(n)))' ...
						(pCells(n,3)-windw(n))*ones(2*windw(n)+1, 1)]);
				end
				if (pCellsNew(n, 3) + windw(n) < Sz(n+1,1))
					inNewnotOld = cat(1, inNewnotOld, ...
						[j*ones(2*windw(n)+1, 1) ...
						((pCellsNew(n,2)-windw(n)):(pCellsNew(n,2)+windw(n)))' ...
						(pCellsNew(n,3)+windw(n))*ones(2*windw(n)+1, 1)]);
				end
			end
		elseif (incP(m) < 0)
			pCellsNew(n, m+1) = pCellsNew(n, m+1) - 1;
			if m == 1
				if (pCells(n, 2) + windw(n) < Sz(n+1,1))
					inOldnotNew = cat(1, inOldnotNew, ...
						[j*ones(2*windw(n)+1, 1) ...
						(pCells(n,2)+windw(n))*ones(2*windw(n)+1, 1) ...
						((pCells(n,3)-windw(n)):(pCells(n,3)+windw(n)))']);
				end
				if (pCellsNew(n, 2) - windw(n) >= 0)
					inNewnotOld = cat(1, inNewnotOld, ...
						[j*ones(2*windw(n)+1, 1) ...
						(pCellsNew(n,2)-windw(n))*ones(2*windw(n)+1, 1) ...
						((pCellsNew(n,3)-windw(n)):(pCellsNew(n,3)+windw(n)))']);
				end
			else
				if (pCells(n, 3) + windw(n) < Sz(n+1,1))
					inOldnotNew = cat(1, inOldnotNew, ...
						[j*ones(2*windw(n)+1, 1) ...
						((pCells(n,2)-windw(n)):(pCells(n,2)+windw(n)))' ...
						(pCells(n,3)+windw(n))*ones(2*windw(n)+1, 1)]);
				end
				if (pCellsNew(n, 3) - windw(n) >= 0)
					inNewnotOld = cat(1, inNewnotOld, ...
						[j*ones(2*windw(n)+1, 1) ...
						((pCellsNew(n,2)-windw(n)):(pCellsNew(n,2)+windw(n)))' ...
						(pCellsNew(n,3)-windw(n))*ones(2*windw(n)+1, 1)]);
				end
			end
		end
	end
% 	disp([j incP'])
end

remOld = [];
for m = 1:size(inOldnotNew,1)
	n = inOldnotNew(m,1) + N + 1;
	if (inOldnotNew(m,2) < 0) || (inOldnotNew(m,2) >= Sz(n+1,1)) || ...
			(inOldnotNew(m,3) < 0) || (inOldnotNew(m,3) >= Sz(n+1,1))
		remOld = cat(1, remOld, m);
	end
	if ~ismember(inOldnotNew(m,:), ANzr, 'rows')
		remOld = cat(1, remOld, m);
	end
end
% AOldT is not the same as NzrData.ANzr, in general it contains fewer
% elements, i.e., elements which are themselves zero (without setting
% them so) are not in AOldT (they are in NzrData.ANzr)
inOldnotNew(remOld, :) = [];

if numel(inOldnotNew)
	ANzr = setdiff(ANzr, inOldnotNew, 'rows');
end

remNew = [];
for m = 1:size(inNewnotOld,1)
	n = inNewnotOld(m,1) + N + 1;
	if (inNewnotOld(m,2) < 0) || (inNewnotOld(m,2) >= Sz(n+1,1)) || ...
			(inNewnotOld(m,3) < 0) || (inNewnotOld(m,3) >= Sz(n+1,1))
		remNew = cat(1, remNew, m);
	end
	if ismember(inNewnotOld(m,:), ANzr, 'rows')
		remNew = cat(1, remNew, m);
	end
end
inNewnotOld(remNew, :) = [];

NzrDataNew.pCells	= pCellsNew;
NzrDataNew.ANzr		= cat(1, ANzr, inNewnotOld);