global Width  nbe connec nnd
global geom ngp traction_loads nel nne

traction_loads = zeros(nnd,2);

n_bc = zeros(6,nbe);
for j = 1:nbe
    for i = 1:nel
        [coord,g] = elem_q4(i);
        p = 1;
        for k=1: nne
            
            if coord(k,2)== Width
                n_bc(p,j) = connec(i,k);
                p = p+1;
            end  
        end       
    end           
end                               % Compute nodal boundary force vector
for i = 1:nbe
    
        ft      = [0. 0. 0. 0.]';       % initialize nodal boundary force vector
        node1   = n_bc(1,i);        % first node
        node2   = n_bc(2,i);        % second node
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
            
            ft  = ft + samp(j,2)*N' *T *J;      % nodal bounbdary force vector
        end

        % Assemble nodal boundary force vector
        
        traction_loads(node1,1) = ft(1);
        traction_loads(node1,2) = ft(2);
        traction_loads(node2,1) = ft(3);
        traction_loads(node2,2) = ft(4);
        
        
        
       
        
        
end    
