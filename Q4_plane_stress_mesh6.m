% THIS PROGRAM USES AN 4-NODED QUADRILATERAL ELEMENT FOR THE LINEAR ELASTIC 
% STATIC ANALYSIS OF A TWO DIMENSIONAL PROBLEM
%
% Make these variables global so they can be shared by other functions
%
clc
clear all
global nnd nel nne nodof eldof n ngp
global geom connec dee nf Nodal_loads
global Length Width NXE NYE X_origin Y_origin dhx dhy 
global traction_loads

%
 format long g
fid =fopen('bar_6_results.txt','w');
disp('Results printed to file : bar_6_results.txt');
%
% To change the size of the mesh, alter the next statements
%    
Length = 4000.; % Length of the model
Width =400;    % Width
NXE = 50.;      % Number of rows in the x direction
NYE = 5;      % Number of rows in the y direction
dhx = Length/NXE; % Element size in the x direction
dhy = Width/NYE;  % Element size in the y direction
X_origin = 0. ;  % X origin of the global coordinate system
Y_origin = Width/2. ;   % Y origin of the global coordinate system
%
nne = 4;
nodof = 2;
eldof = nne*nodof;
ngp = 2;
%
Q4_mesh     % Generate the mesh 
%
E = 3*10^5;     % Elastic modulus in MPa
vu = 0.3;       % Poisson's ratio 
thick =1;      % Beam thickness in mm
%
% Form the elastic matrix for plane stress 
%
dee = formdsig(E,vu);
%
%
% Boundary conditions
%
nf = ones(nnd, nodof);    % Initialise the matrix nf to 1
%
% Restrain in all directions the nodes situated @ 
% (x = Length)
%
for i=1:nnd
    if geom(i,1) == 0
        nf(i,:) = [0 0];
    end
end
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
 % N
%

traction_loads = zeros(nnd,2);

p = 1;
for i = 1:nel
        [coord,g] = elem_q4(i);
        
       
        for k=1: nne
            
            if coord(k,2)== Width/2
                n_bc(p,1) = connec(i,k);
                p = p + 1;
            end
            
        end
        
end
n_bc  = reshape(n_bc,[2,size(n_bc,1)/2]);
l_bc = zeros(4,size(n_bc,2));
for i = 1:size(n_bc,2)
    l_bc(:,i) = [0; -1; 0; -1];
    
end
n_bc = [n_bc;l_bc];
    
  

% Compute nodal boundary force vector
for i = 1:size(n_bc,2)
    
        ft      = [0. 0. 0. 0.]';       % initialize nodal boundary force vector
        node2   = n_bc(1,i);        % first node
        node1   = n_bc(2,i);        % second node
        n_bce   = n_bc(3:6,i);      % traction value at node1
                        
        x1 = geom(node1,1); y1=geom(node1,2);    % coord of the first node
        x2 = geom(node2,1); y2=geom(node2,1);    % coord of the second node
    
        leng = sqrt((x2-x1)^2 + (y2-y1)^2);  % edge length
        J    = leng/2;                       % Jacobian 
        
        samp = gauss(ngp);                % get gauss points and weights
       
        for j=1:ngp                         % integrate in psi direction (1D integration)  
               
            psi = samp(j,1);               
            N   = 0.5*[1-psi    0      1+psi      0;       % 1D shape functions x-component
                        0     1-psi      0      1+psi];    % 1D shape functions y-component  
                                 
            T    = N * n_bce;
            
            ft  = ft + samp(j,2)*N' *T *J;
            
        end
        traction_loads(node1,1) = traction_loads(node1,1)+ft(1);
        traction_loads(node1,2) = traction_loads(node1,2)+ft(2);
        traction_loads(node2,1) = traction_loads(node2,1)+ft(3);
        traction_loads(node2,2) = traction_loads(node2,2)+ft(4);       % nodal bounbdary force vector

        % Assemble nodal boundary force vector
        
              
end    

    

%
%%%%%%%%%%%%%%%%%%%%%%%%%% End of input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% Assemble the global force vector
%
fg=zeros(n,1);
for i=1: nnd
    if nf(i,1) ~= 0   
        fg(nf(i,1))= Nodal_loads(i,1) + traction_loads(i,1);
    end
    if nf(i,2) ~= 0   
        fg(nf(i,2))= Nodal_loads(i,2) + traction_loads(i,2);
    end
end 
%
%  Form the matrix containing the abscissas and the weights of Gauss points
%
print_model
ngp = 2;
samp=gauss(ngp);
%
% Numerical integration and assembly of the global stiffness matrix
%
%  initialise the global stffness matrix to zero
kk = zeros(n, n);
%
for i=1:nel
    [coord,g] = elem_q4(i) ;       % coordinates of the nodes of element i,
                                   %  and its steering vector 
    ke=zeros(eldof,eldof) ;        % Initialise the element stiffness 
                                   %  matrix to zero
    for ig=1: ngp
        wi = samp(ig,2);
    for jg=1: ngp
        wj=samp(jg,2);
        [der,fun] = fmlin(samp, ig,jg);  % Derivative of shape functions 
                                         % in local coordinates
        jac=der*coord;                   % Compute Jacobian matrix
        d=det(jac);                      % Compute determinant of Jacobian
                                         %  matrix
        jac1=inv(jac);                   % Compute inverse of the Jacobian
        deriv=jac1*der;                  % Derivative of shape functions in
                                         % global coordinates
        bee=formbee(deriv,nne,eldof);    %  Form matrix [B]
        ke=ke + d*thick*wi*wj*bee'*dee*bee; % Integrate stiffness matrix

        
    end
    end
    kk=form_kk(kk,ke, g);               % assemble global stiffness matrix
end
%
%
%%%%%%%%%%%%%%%%%%%%%%%  End of assembly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
delta = kk\fg ;                         % solve for unknown displacements

disp('node       x_disp       y_disp ')     %
for i=1: nnd                                %
    if nf(i,1) == 0                         %
        x_disp =0.;                         %
    else
        x_disp = delta(nf(i,1));            %
    end
%    
    if nf(i,2) == 0                         %
        y_disp = 0.;                        %
    else
        y_disp = delta(nf(i,2));            %
    end
disp([i  x_disp  y_disp])                   % Display displacements of each node
DISP(i,:) = [x_disp  y_disp];
end

%
%
ngp=1;                              % Calculate stresses and strains at 
                                    %the centre of each element
samp=gauss(ngp);
%
for i=1:nel
    [coord,g] = elem_q4(i);      % coordinates of the nodes of element i, 
                                 % and its steering vector 
    eld=zeros(eldof,1);          % Initialise element displacement to zero
    for m=1:eldof                %
        if g(m)==0               %
            eld(m)=0.;           %
        else                     %
            eld(m)=delta(g(m));  % Retrieve element displacement from the 
                                 % global displacement vector
        end
    end
%
    for ig=1: ngp
        wi = samp(ig,2);
    for jg=1: ngp
        wj=samp(jg,2);
        [der,fun] = fmlin(samp, ig,jg); % Derivative of shape functions in
                                        % local coordinates
        jac=der*coord;                  % Compute Jacobian matrix
        jac1=inv(jac);                  % Compute inverse of the Jacobian
        deriv=jac1*der;                 % Derivative of shape functions in 
                                        % global coordinates
        bee=formbee(deriv,nne,eldof);   % Form matrix [B]
        eps=bee*eld                     % Compute strains
        sigma=dee*eps                   % Compute stresses
    end
    end
   SIGMA(i,:)=sigma ;           % Store stresses for all elements
end
%
% Average stresses at nodes
%
[ZX, ZY, ZT, Z1, Z2]=stresses_at_nodes_Q4(SIGMA);
%
%
% Plot stresses in the x_direction 
%
print_bar_result;
fclose(fid);
U2 = DISP(:,2);
cmin = min(U2);
cmax = max(U2);
caxis([cmin cmax]);
patch('Faces', connec, 'Vertices', geom, 'FaceVertexCData',U2, ...
      'Facecolor','interp','Marker','.');
colorbar;
figure(2)
p=1;
for i = 1:nnd
    if geom(i,2)== -Width/2
        coor(p,:)=geom(i,:);
        node(p,:) = (DISP(i,:));
        
        str(p,1)=ZX(i,1)
        p = p+1;
    end
end
subplot(2,1,1)
plot(coor(:,1),node(:,2),'r')
hold on
num = xlsread('disp.xlsx')
plot(num(:,1),num(:,2),'k')
legend('MATLAB','ABAQUS')

ylabel('Displacement U2')
subplot(2,1,2)
plot(coor(:,1),str(:,1),'r')
hold on
num = xlsread('str.xlsx')
plot(num(:,1),num(:,2),'k')
legend('MATLAB','ABAQUS')
xlabel('Distance along Length of Plate ')
ylabel('stress-xx')
