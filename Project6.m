%% Author: Alec Sblendorio 
%% Due Date: December 6, 2020
%% Subject: EP 501 Numerical Methods for Scientists and Engineers 
%% Project #6
%%Problem 1: Electrostatic Potential in Dielectrics with Boundary
%%Conditions
%%Problem #2: Application of Ordinary Differential Equations
% Collaboration with Dennis Turano, and Bartosz Kwiecisnki, and Joseph
% Akash

%% Problem #1: Electrostatic Potential in Dielectrics with Boundary Conditions
%%Part A: Plot the dielectric function and note that it varies rapidly near the edges of the domain of interest.
%%Part B: Develop a system of finite difference equations for this system 
%%based on second order accurate centered differences (numerical differences
%%may also be used for dielectric function derivatives). 
%%Include this system in your homework submission.

%%Part C: Develop two finite difference equations for the boundary conditions
%%of this system. Use a first- order forward difference at the x = −a boundary.
%%Include these equations with your submission.

%%Part D: Solve your system of equations using the MATLAB “\” operator.
%%Part E: Since the dielectric function varies rapidly at the boundary, 
%%this is a problem where a second order (forward) difference may be useful
%%(see course notes for formula). Reformulate your matrix system to include
%%this for the x = −a boundary and solve the system numerically for 20 
%%grids points. Plot the result and compare it against the solution with a
%%first order forward difference.

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER END%%%%%%%%%%%%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER PLOTS%%%%%%%%%%%%%%%%%%');
%% Problem #2: Application of Ordinary Differential Equations
%%Part A:  Solve this system with the RK4 algorithm and plot the velocity
%%and position as a function of time. Compare this against the RK2 solution
%%in the repository and show they are roughly equal. Use 2πm 100 time steps
%%per particle oscillation period T = qB .

%%Part B: Show that the RK4 solution is better in the sense that it can
%%solve the problem accurately with fewer time steps.

%%Part C: Suppose the magnetic field is changed to vary linearly in the y-direction:
%%B(y)=B 1+(1/2)y in e_z
%%Use your RK4 solver to find the velocity. Plot the velocities and path of 
%%the particle in the x-y plane for a least five periods of oscillation. 
%%HINT: The particle should execute trochoidal motion.

disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER END%%%%%%%%%%%%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER PLOTS%%%%%%%%%%%%%%%%%%');


%% END PROJECT #6