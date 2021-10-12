N_xi = 3;
N_eta = 3;
NumNodes = N_xi * N_eta
Node = zeros(N_xi, N_eta);                          
Node(1 : NumNodes) = [1:NumNodes];                  

disp(Node)
