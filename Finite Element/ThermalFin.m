function [u, Troot] = ThermalFin(mesh, mu)
    %
    % -------------------------------------------------------------------------
    %   Computes the temperature distribution and root temperature for a fin
    %   using the Finite Element Method.
    % -------------------------------------------------------------------------
    %
    %   INPUT   mesh    grid label (coarse, medium, or fine)
    %           mu      thermal conductivities of sections and Biot number 1x5
    %
    %   OUTPUT  u       temperature disctribution in the fin
    %           Troot   root temperature
    %
    % -------------------------------------------------------------------------
    
    
    % Parameters setup: rearranges the values in mu
    kappa = ones(6,1);
    kappa(1:4) = mu(1:4);   % kappa(5) = 1
    kappa(6) = mu(5);       % Bi = mu(5)
    
    % Initialization
    A = sparse(mesh.nodes, mesh.nodes);
    F = sparse(mesh.nodes,1);
    
    % Domain Interior
    for i = 1:5         % interior regions
        for n = 1:length(mesh.theta{i})
            phi = mesh.theta{i}(n,:)';
            x_alp = mesh.coor(phi,:);           % Obtain global x,y coordinates 
            x_ones = ones(3,1);
            x_tilde = cat(2, x_ones, x_alp); 
            x_eye = eye(3);
            c_tilde = x_tilde \ x_eye;          % Solves for cx,alp , cy,alp
            
            Area_k = 0.5 * ((x_alp(1,1) * (x_alp(2,2) - x_alp(3,2))) ... 
                + (x_alp(2,1) * (x_alp(3,2) - x_alp(1,2))) + (x_alp(3,1) * (x_alp(1,2) - x_alp(2,2))));     % Calculates area of triangulation
            
            Alocal = zeros(3,3);
            
            for j1 = 1:3
                for j2 = 1:3
                    Alocal(j1,j2) = (c_tilde(j1,2) * c_tilde(j2,2)) + (c_tilde(j1,3) * c_tilde(j2,3));      % Local triangulation matrix
                end
            end
            
            Alocal = Alocal * Area_k * kappa(i);                % [ 3x3 matrix: 3 nodes per triangle ]
            A(phi,phi) = A(phi,phi) + Alocal;    
        end
    end
    
    
    % Boundaries (not root)
    i = 6;
        for n = 1:length(mesh.theta{i})
            phi = mesh.theta{i}(n,:)';
            x_alp = mesh.coor(phi,:);                            % Obtain global x,y coordinates 
            x_l2norm = norm(x_alp(1,:) - x_alp(2,:)); 
            Alocal = [2, 1; 1, 2];                              % [ 2x2 matrix: 2 nodes per edge ]
            Alocal = (x_l2norm * kappa(i) / 6.0) * Alocal;
            A(phi,phi) = A(phi,phi) + Alocal;
        end
    
    % Root Boundary
    i = 7;
        for n = 1:length(mesh.theta{i})
            phi = mesh.theta{i}(n,:)';
            x_alp = mesh.coor(phi,:);                            % Obtain global x,y coordinates 
            x_l2norm = norm(x_alp(1,:) - x_alp(2,:)); 
            Flocal = (x_l2norm / 2.0) * ones(2,1);              % [ 2x1 matrix: 2 nodes per edge ] 
            F(phi) = F(phi) + Flocal;
        end
    
    
    % Solve for the temperature distribution
    u = full(A\F);      % use full since A and F are sparse
    
    
    % Compute the temperature at the root
    Troot = 0; 
    
    i = 7;
    for n = 1:length(mesh.theta{i})
        phi = mesh.theta{i}(n,:)';
        x_alp = mesh.coor(phi,:);                            % Obtain global x,y coordinates 
        x_l2norm = norm(x_alp(1,:) - x_alp(2,:)); 
    end
    
end