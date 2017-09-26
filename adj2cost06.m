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
1. Raghvendra V. Cowlagi
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

Description: Transition cost matrix from adjacency matrix and intensity map
%}

% ************************** TRANSITION COST MATRIX ***********************
function G = adj2cost06(G, V, nodeG)

K1 = 1;
K2 = 0.1;
K3 = 0;
for n = 1:size(V, 1)
	nhbrs = find(G(n,:));
	for m = nhbrs
% 		G(n,m) = K1*((256 - V(m,4))^2)*(V(m,3))^2 + K2*(V(m,3))^2;				% G(n,m) = cost of going from n to m
% 		G(n,m) = K1*((256 - V(m,4)))*(V(m,3))^2 + K2*(V(m,3))^2;				% G(n,m) = cost of going from n to m
		if V(m, 4) < 1e-5
			G(n, m) = 1e8;
		else
% 			G(n,m) = K2*(V(m,3))^2 + K1*((256 - V(m,4)));
			G(n,m) = (K1*((256.1 - V(m,4)))*V(m,3) + K2*(V(m,3))^2)/10 ...
				+ K3*norm(V(m, 1:2) - V(nodeG, 1:2));
		end
% 		G(n,m) = K1*((256 - V(m,4)) + (256 - V(n,4)))/2 + K2*(V(m,3))^2;
% 		G(n,m) = K1*(V(m,4))*(V(m,3))^2 + K2*(V(m,3))^2;				% G(n,m) = cost of going from n to m
	end
end