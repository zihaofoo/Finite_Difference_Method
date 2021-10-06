function [Solution, FlowRate, I_xx] = ChannelFlow(N_xi, N_eta, bb, hh)
%==================================================================================
% ChannelFlow.m
%==================================================================================
% The following code is the solution to the channel flow problem outlined in PS#1
% in the 16.920J/2.097J/6.339J course. Since this code is intended as a skeleton
% code and should therefore be easy to read, we chose not to implement the code 
% using vectorization. This will cause slightly lower performance.
%
% Inputs:
% =======
%    N_xi : Number Of nodes in the xi direction.
%    N_eta: Number of nodes int he eta direction
%    bb : The length of the base
%    hh : the height of the channel
%
% Outputs:
% ========
%    Solution : A xi-eta matrix containing the solution
%    Flowrate : The flowrate throught eh channel
%    I_xx     : The moment of inertia of the channel
%
% Some additional information for the function:
%
% In this pseudo code, the following conventions are used:
%
%    1) The computational domain xi-eta has values of (0,0) corresponding to the 
%       bottom left corner of the physical and computational domain. Correspondingly,
%       the upper-right corner in the xi-eta domain has a (xi, eta) value of (1,1). 
%
%    2) We use a node map in the matrix 'Node'. Node(i,j) grabs node reference 
%       number. Using the Node matrix, we can easily stamp/stencil the related
%       values into the finite difference matrix. To determine how the Node matrix
%       looks, you can easily type the creation commands in the command prompt.
%
% Written by: D.J.Willis
%==================================================================================
close all;

%-----------------------------------------------------------------------------------
% Geometric Parameters.
%-----------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------
% Grid Details and parameters
%-----------------------------------------------------------------------------------
d_xi  = 1./(N_xi-1);                                 
d_eta = 1./(N_eta-1);                                
NumNodes = N_xi * N_eta;                             

%-----------------------------------------------------------------------------------
% Initializing the sparse matrix 'A'
%-----------------------------------------------------------------------------------
% Note: We use spalloc to "sparse-allocate" the A matrix. It is essential to allo-
% cate a sparse A-matrix due to memory restrictions. 
%-----------------------------------------------------------------------------------
A   = spalloc(NumNodes, NumNodes, 9*NumNodes);       

%-----------------------------------------------------------------------------------
% Initializing the RHS.
%-----------------------------------------------------------------------------------
RHS =                 

%-----------------------------------------------------------------------------------
% Creation of a node numbering scheme. The node numbering scheme is created here
% in order to simplify the overall solution process. The idea is as follows:
%
% We construct a matrix called Node, which has elements corresponding to the node
% numbers in the grid representation of the solution domain. This allows us to cycle
% through the resulting matrix, grab element (i,j) and easily find the (i +/- 1), and 
% (j +/- 1) node numbers. 
%-----------------------------------------------------------------------------------
Node = zeros(N_xi, N_eta);                          
Node(1 : NumNodes) = [1:NumNodes];                  


'Constructing The Jacobian'

for i = 1 : N_xi 
    for j = 1 : N_eta 
        xi(i,j)  =  
        eta(i,j) =    
        J(i,j)   =  
    end
end

'Constructing the "A" Matrix'
%-----------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------
%                            INNER REGION OF THE DOMAIN                   
%-----------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------
% We begin our construction of the matrix 'A' by considering the Inner region of the
% domain. This is essentially the fill-in for all nodes not touching the boundary.
% The boundary nodes are handled later.
%-----------------------------------------------------------------------------------
for i = 2 : N_xi - 1
    for j = 2 : N_eta - 1   
        ANode_i = Node(i,j);                        % Setting A_Matrix position for node i,j        
		%-----------------------------------------------------------------------------------
		%------------------------------ The Transformation ---------------------------------
		%-----------------------------------------------------------------------------------
		% The various components of the transformation  
		%-----------------------------------------------------------------------------------        
        
		a =  
		b =  
		c =  
		alpha =  
		beta  =  
		d =  
		e =  
        
		%-----------------------------------------------------------------------------------
		%---------------------------- FILLING UP THE A MATRIX ------------------------------
		%-----------------------------------------------------------------------------------
		% The filling of the matrix is done via the stamping of the computational molecule 
		% in the appropriate parts of the A-Matrix
		%-----------------------------------------------------------------------------------
		
		%-----------------------------------------------------------------------------------
		%---------------------  RHS Part Of Computational Molecule  ------------------------
		%-----------------------------------------------------------------------------------
		% Note, A(p,q), here is: p=the current node number on the grid, at which the 
		% differential equation is being approximated, and q refers to the neighboring
		% point to the current node. So, here we are using a stamping stencil based on the
		% ANode_i matrix. It may be worthwhile to view a reduced dimension version of 
		% ANode_i to fully grasp what is happening here.
		%-----------------------------------------------------------------------------------
        
        A(ANode_i, Node(i+1, j+1) ) =  
        A(ANode_i, Node(i+1, j  ) ) =  
	    A(ANode_i, Node(i+1, j-1) ) =  
      
        %-----------------------------------------------------------------------------------
        %--------------------  Middle Part Of Computational Molecule  ----------------------
        %-----------------------------------------------------------------------------------
        
        A(ANode_i, Node(i  , j+1) ) =  
	    A(ANode_i, Node(i  , j  ) ) =  
	    A(ANode_i, Node(i  , j-1) ) =  
                                    
        %-----------------------------------------------------------------------------------
        %---------------------  LHS  Part Of Computational Molecule  -----------------------      
        %-----------------------------------------------------------------------------------
        
        A(ANode_i, Node(i-1, j+1) ) =  
	    A(ANode_i, Node(i-1, j  ) ) =  
	    A(ANode_i, Node(i-1, j-1) ) =                  
        
	end
end

%-----------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------
%-------------------------- BOUNDARY CONDITIONS ------------------------------------
%-----------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------
%------------------------ Bottom of the Domain -------------------------------------
%-----------------------------------------------------------------------------------
j = 1;
for i=2:N_xi-1
    ANode_i = Node(i,j);
    A(ANode_i, Node(i,j))= 
    RHS(ANode_i)= 
end
%-----------------------------------------------------------------------------------
%------------------------ Top of the Domain ----------------------------------------
%-----------------------------------------------------------------------------------
j = N_eta;
for i=2:N_xi-1
    ANode_i = Node(i,j);
        
    % Fill in the boundary condition for the top of the domain
    
end
%-----------------------------------------------------------------------------------
%--------------------------- Left Side of the Domain -------------------------------
%-----------------------------------------------------------------------------------
i = 1;
for j=2:N_eta-1
    ANode_i = Node(i,j);
    
    % Fill in the boundary condition for the left of the domain
    
end
%-----------------------------------------------------------------------------------
%-------------------------- Right Side of the Domain -------------------------------
%-----------------------------------------------------------------------------------
i = N_xi;
for j=2:N_eta-1
    ANode_i = Node(i,j);

    % Fill in the boundary condition for the Right side of the domain
    
end

%-----------------------------------------------------------------------------------
%------------------------------  Domain Corners ------------------------------------
%-----------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------
%--------------------------------- BOTTOM LEFT -------------------------------------
%-----------------------------------------------------------------------------------

    ANode_i = Node(1,1);
    % Fill in the boundary condition for the bottom left corner of the domain

%-----------------------------------------------------------------------------------
%--------------------------------- Bottom Right ------------------------------------
%-----------------------------------------------------------------------------------
    ANode_i = Node(N_xi,1);
    
    % Fill in the boundary condition for the bottom right corner of the domain
    
%-----------------------------------------------------------------------------------
%----------------------------------- TOP Left --------------------------------------
%-----------------------------------------------------------------------------------

    ANode_i = Node(1,N_eta);  % Setting A_Matrix position for node i,j
    
    % Fill in the boundary condition for the top left corner of the domain
   
%-----------------------------------------------------------------------------------
%----------------------------------- TOP RIGHT  ------------------------------------
%-----------------------------------------------------------------------------------
    ANode_i = Node(N_xi,N_eta);
    
    % Fill in the boundary condition for the top right corner of the domain

    'Solving the system'
%-----------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLUTION OF Ax=b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------
Sol = A\RHS;



%-----------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------
%------------------------------- Post Processing the Solution ----------------------
%-----------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------
'Post-processing'
%-----------------------------------------------------------------------------------
%-------------------------- Computing the flow rate & I ----------------------------
%-----------------------------------------------------------------------------------
% With the velocity known, we can compute the flowrate and the moment of inertia. As
% a check of your flow rate integral try computing the area of the channel and compare 
% that with the analytical result. That is a good check.
%-----------------------------------------------------------------------------------


% Fill In Your Computations For the Flow Rate 


%-----------------------------------------------------------------------------------
%----------------------------- Moment of Inertia -----------------------------------
%-----------------------------------------------------------------------------------

% Fill In Your Computations For the Moment Of Inertia


%========================= Information for the plots ===============================
% Fill-in the transformation equations for x(xi,eta) and y(xi,eta)
for i = 1 : N_xi 
    for j = 1 : N_eta     
        xi  = (i-1)*d_xi;
        eta = (j-1)*d_eta;      
        x(i,j)=
        y(i,j)=
    end
end

%-----------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------
%                                SOLUTION OUTPUT
%-----------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------
% Uncomment the following lines to get the plot done
% %===================================================================================
% % Plot # 1
% %===================================================================================
% figure
% mesh(x,y,0*x,0*y);
% view([0,0,1]);
% axis equal
% 
% %===================================================================================
% % Plot # 2
% %===================================================================================
% figure
% [ccc,fff]=contour(x,y,abs(Solution),(0.02:0.02:max(max(abs(Solution))))');
% clabel(ccc,fff);
% hold on
% patch([0,bb/2,bb/2+aa ,0   ],[0,0,hh,hh],-ones(1,4),0,'facecolor',[0.8,0.8,0.8]);
% axis equal
