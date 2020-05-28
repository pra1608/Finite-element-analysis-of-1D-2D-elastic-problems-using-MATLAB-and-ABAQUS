%  File:    Q4_COARSE_MESH_DATA
%
global nnd nel nne nodof eldof n ngp
global geom connec dee nf Nodal_loads
%
% To change the size of the mesh, alter the next statements
%    
nnd = 28 ;                % Number of nodes:  				
nel = 18;                  % Number of elements: 				
nne = 4 ;                % Number of nodes per element:			
nodof =2;                % Number of degrees of freedom per node
ngp = 2;                  % number of Gauss points
eldof = nne*nodof;       % Number of degrees of freedom per element 	
%
%
% Nodes coordinates x and y  
geom = [0.0   -0.5; ...         %   x and y coordinates of node 1
        0.0   -1/6; ...         %   x and y coordinates of node 2
        0.0    1/6; ...         %   x and y coordinates of node 3
        0.0    0.5; ...
        1/3   -5/12; ...         %   x and y coordinates of node 4
        1/3   -1/9; ...         %   x and y coordinates of node 5
        1/3    7/36; ...         %   x and y coordinates of node 6
        1/3    1/2; ...   
        2/3   -1/3; ...         %   x and y coordinates of node 7
        2/3   -1/18; ...
        2/3    2/9; ...         %   x and y coordinates of node 8
        2/3    0.5; ...         %   x and y coordinates of node 9
         1    -1/4; ...         %   x and y coordinates of node 10
         1     0;    ...
         1     0.25; ...         %   x and y coordinates of node 11
         1     0.5; ...         %   x and y coordinates of node 12
        4/3   -1/6; ...         %   x and y coordinates of node 13
        4/3    1/18; ...
        4/3    5/18; ...         %   x and y coordinates of node 14
        4/3    0.5; ...         %   x and y coordinates of node 15
        5/3   -1/12; ...         %   x and y coordinates of node 16
        5/3    1/9; ...
        5/3    11/36; ...         %   x and y coordinates of node 17
        5/3    0.5; ...         %   x and y coordinates of node 18
         2      0; ...         %   x and y coordinates of node 19
         2     1/6; ...
         2     1/3; ...
         2     0.5];         %   x and y coordinates of node 20
                    %   x and y coordinates of node 21
%       
%      
%
disp ('Nodes X-Y coordinates')
geom
%
% Element connectivity
connec= [ 1   5    6    2 ;...   % Element 1  
          2   6    7    3 ;...   % Element 2  
          3   7    8    4 ;...   % Element 3 
          5   9    10   6 ;...   % Element 4
          6   10   11   7 ;...   % Element 5  
          7   11   12   8 ;...   % Element 6
          9   13   14   10 ;...   % Element 7
          10  14   15   11 ;...   % Element 8
          11  15   16   12 ;...   % Element 9
          13  17   18   14 ;...   % Element 10
          14  18   19   15 ;...    %11
          15  19   20   16 ;...   % Element 12
          17  21   22   18 ;...      % Element 13
          18  22   23   19 ;...              %14
          19  23   24   20 ;...              %15
          21  25   26   22 ;...              %16
          22  26   27   23 ;...              %17
          23  27   28   24];                 %18
     %   
%  
disp ('Elements connectivity')
connec
%
E = 3*10^7.;     % Elastic modulus in Pa
vu = 0.3;       % Poisson's ratio 
thick = 1;      % Beam thickness in m
%
% Form the elastic matrix for plane stress 
%
dee = formdsig(E,vu);
%
%
% Boundary conditions
%
nf = ones(nnd, nodof);    % Initialise the matrix nf to 1
nf(1,:) = [0   0];       % Node 19 is restrained in the x and y directions
nf(2,:) = [0   0];       % Node 20 is restrained in the x and y directions
nf(3,:) = [0   0];
nf(4,:) = [0   0];        % Node 21 is restrained in the x and y directions

%
% Counting of the free degrees of freedom
%
n=0;
for i=1:nnd
    for j=1:nodof
        if nf(i,j) ~= 0 
            n=n+1;
           nf(i,j)=n;
        end
    end
end
%
% loading
%
Nodal_loads= zeros(nnd, 2);   % Initialise the matrix of odal loads to 0
%
% Apply a concentrated at the node having x = 0, and y = 0.
%
for i = 1:nnd
    if (geom(i,1)==0&&geom(i,2)==0.5)||(geom(i,1)==2.&&geom(i,2)==0.5)
        Nodal_loads(i,:) = [0,-10/3];
    elseif geom(i,2)==0.5
        Nodal_loads(i,:) = [0,-20/3];
            
    end
    
end