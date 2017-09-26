Reference:
----------
R. V. Cowlagi and P. Tsiotras, "Multiresolution Motion Planning for Autonomous
Agents via Wavelet-Based Cell Decompositions." IEEE Transactions on Systems, Man, and 
Cybernetics- Part B: Cybernetics, vol. 42, no. 5, pp. 1455 - 1469, October 2012.
doi: 10.1109/TSMCB.2012.2192268

Instructions for sample code:
-----------------------------
1. Run the file "mr_demo_map.m" to execute the example shown in Fig. 16 - 17 in the paper.
This example demonstrates the execution of the proposed motion-planning algorithm in an 
environment with varying levels of "threat" to be avoided. This simulation will provide 
two views of the algorithm during execution: a "local" view indicating the execution
of the motion-planner in the fine resolution window (which considers a dynamical model of 
the vehicle), and a "global" view indicating the path in the multi-resolution map (where
vehicle dynamical constraints are ignored).

Note: This code was written in MATLAB R2009b. Running this code in newer versions is 
likely to be slow, with several warnings visible, primarily because of significant changes
in the MATLAB optimization toolbox.


2. Run the file "mr_demo_culdesac.m" to execute the example shown in Fig. 9 - 10 in the paper.
This example is an illustration of the behavior of the algorithm in recovery from a cul-de-sac.
Proposition 1 on p. 1459 in the article states the proposed multi-resolution path-planning
algorithm is complete, i.e., if a path to the goal exists, then this algorithm will find it. 
This example illustrates how a path to the goal is found even if a cul-de-sac is encountered
due to the multi-resolution approximation of the environment map.

3. Run the file  "mr_comparison.m" (with appropriate modifications given in code comments)
to regenerate the comparison data shown in Fig. 15.