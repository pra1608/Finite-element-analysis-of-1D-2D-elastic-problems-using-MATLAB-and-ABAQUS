%
fprintf(fid, '-------------------------------------------------------- \n');
fprintf(fid, ' \n\n\n ******* PRINTING ANALYSIS RESULTS **************\n\n\n');
%
% Print global force vector
%
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid,'Global force vector  F \n');
fprintf(fid,'   %g\n',fg);
fprintf(fid,'\n');
%
%
% Print Displacement solution vector
%
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid,'Displacement solution vector:  delta \n');
fprintf(fid,' %8.5f\n',delta);
fprintf(fid,'\n');
%
% Print nodal displacements
%
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid, 'Nodal displacements \n');
fprintf(fid, 'Node                    disp_x                            disp_y\n');
for i=1:nnd
fprintf(fid,' %g,                     %8.5f                              %8.5f\n',i,DISP(i,1),DISP(i,2));
end
fprintf(fid,'\n');
%
% Print Members actions
fprintf(fid,'--------------------------------------------------------------------------------- \n');
fprintf(fid,' Print stresses at the Stresses at points \n');
fprintf(fid,'node          Stress_xx              stress_yy                  stress_xy                       stress_1                            stress_2\n');

for i = 1:nnd
    fprintf(fid,'%g            %08.3f                 %08.3f                    %08.3f                          %08.3f                            %08.3f\n',i,ZX(i,1),ZY(i,1), ZT(i,1),Z1(i,1),Z2(i,1) );
end
                                                                                                     
    
    