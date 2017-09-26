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
function drawCdV05(VCell, myAxes, myCol, myWd, txt, someCells, colFace, ...
	imgSize, map)

axes(myAxes); hold on;
if (nargin <= 5) || (numel(someCells) == 0)
	drawCells = 1:size(VCell, 1);
else
	drawCells = someCells;
end

if (nargin > 6) && (numel(colFace) == 1) && (colFace == 0)
	Z	= zeros(imgSize);
	for n = 1:size(VCell,1)
		Z((VCell(n,2) + 1):(VCell(n,2) + VCell(n,3)), ...
			(VCell(n,1) + 1):(VCell(n,1) + VCell(n,3))) = VCell(n,4);
	end
	image(Z, 'XData', 0.5, 'YData', 0.5); colormap(map);
end

for m = drawCells
	if (nargin > 6) && ((numel(colFace) ~= 1) || (colFace ~= 0))
		rectangle('Position', [VCell(m,1) VCell(m,2) VCell(m,3) ...
			VCell(m,3)], 'EdgeColor', myCol, 'LineWidth', myWd, ...
			'FaceColor', colFace);
	else
		rectangle('Position', [VCell(m,1) VCell(m,2) VCell(m,3) ...
			VCell(m,3)], 'EdgeColor', myCol, 'LineWidth', myWd);
	end
	if txt
		text(VCell(m,1) + 0.1*VCell(m,3), ...
			VCell(m,2) + 0.15*VCell(m,3), num2str(m), 'FontName', ...
			'Consolas', 'FontSize', txt + VCell(m,3), ...
			'FontWeight', 'bold', 'Color', 'k');
	end
end
drawnow;