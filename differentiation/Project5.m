%% Author: Alec Sblendorio 
%% Due Date: December 4, 2020
%% Subject: EP 501 Numerical Methods for Scientists and Engineers 
%% Project #5
%%Problem 1: Numerical Vector Derivatives (curl)
%%Problem #2: Numerical Vector Derivatives (gradient and laplacian)
%%Problem #3: Integration in Multiple Dimensions
%%Problem #4: Line Integration
% Collaboration with Sho Okayama, Dennis Turano, and Bartosz Kwiecisnki
%% Data Input 

%% Problem #1: Numerical Vector Derivatives (curl):
%%Part A: Plot the two components of the vector magnetic field defined by the 
%%piecewise function. Use an image plot (e.g. pcolor and shading flat in
% MATLAB for each magnetic field component (Bx,By) and have your plot 
%%show the region −3a ≤ x ≤ 3a, −3a ≤ y ≤ 3a. Make sure you add a colorbar 
%%and axis labels to your plot. You will need to define a range and resolution
%%in x and y, and create a meshgrid from that. Be sure to use a resolution 
%%fine enough to resolve important variations in this function.

%%Part B: Make a quiver plot of the magnetic field B; add labels, etc

%%Part C: Compute the numerical curl of B, i.e. ∇ × B. Use centered 
%%differences on the interior grid points and first-order derivatives on the edges.
%%Plot your result using imagesc, or pcolor.

%%Part D: Compute ∇ × B analytically (viz. by hand). Plot the alongside your
%%numerical approximation and demonstrate that they are suitably similar.

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');


 

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER END%%%%%%%%%%%%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER PLOTS%%%%%%%%%%%%%%%%%%');
%% Problem #2: Numerical Vector Derivatives (gradient and laplacian)
%%Part A:  Compute and plot the scalar field. and plot this function in the 
%%region−3a≤x≤3a,−3a≤y≤3a in the z=0 plane. Be sure to use a resolution 
%%fine enough to resolve variations in this function 
%%(aside from those associated with the singularity).

%%Part B: Write a function to numerically compute the Laplacian of a scalar field, i.e. ∇2Φ.
%%Plot your result with appropriate labels and colorbars.

%%Part C: Compute an analytical laplacian (viz. differentiate by hand),
%%plot the results alongside your numerical calculation, 
%%and demonstrate that your numerical laplacian is suitably accurate.

disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER END%%%%%%%%%%%%%%%%%%');

%% Problem #3: Integration in Multiple Dimensions.
%%Part A:  Numerically compute the electrostatic energy in the region
%%R ≡ −3a ≤ x ≤ 3a,−3a ≤ y ≤ 3a, −3a ≤ z ≤ 3a, defined by the integral.
%%using an iterated trapezoidal method (sweeps of single dimensional integrations) 
%%or multi-dimensional trapezoidal method program that you write.


disp('%%%%%%%%%%%%%%%%%%PROBLEM #3 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #3 ANSWER END%%%%%%%%%%%%%%%%%%');
%% Problem #4: Line Integration
%%Part A:  Compute and plot the parametric path
%%r(φ)≡x(φ)eˆ +y(φ)eˆ =r cosφeˆ +r sinφeˆ (0≤φ≤2π) xy0x0y
%%in the x,y plane on the same axis as your magnetic field components 
%%from problem 1 (create a new figure which plots the path on top of a 
%%pcolor plot of each component). Take r0 = 2a. You will need to define a 
%%grid in φ to do this.

%%Part B: Plot the two components of the magnetic field B(x(φ), y(φ)) at 
%%the x, y points along r and visually compare against your image plots of
%%the magnetic field and path to verify.

%%Part C: Numerically compute the tangent vector to the path r by performing the derivative:
%%Compare your numerical results against the analytical derivative (e.g. plot the two)
%%and refine your grid in φ (if necessary) such that you get visually 
%%acceptable results - i.e. such that the path appears circular.

%%Part D: Numerically compute the auxiliary magnetic field integrated around the path r, i.e.
%%where the differential path length is given by:

disp('%%%%%%%%%%%%%%%%%%PROBLEM #4 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #4 ANSWER END%%%%%%%%%%%%%%%%%%');

%% END PROJECT #5