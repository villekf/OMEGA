function nbr = AssembleENBR(ny,nx)
% ---------------------------------
% This function is used to assemble the
% nearest neighborhood system for a grid.
% ---------------------------------
% Ville Kolehmainen, Nov. 2001
% ---------------------------------------------------

% Zeros are set "over the grid" pixels

dx = zeros((nx-1)*ny,2);
dy = zeros((ny-1)*nx,2);


n = 0;

for k = 1:nx
    for j = 1:ny-1
      
        n = n + 1; 
        
        dy(n,:) = (k-1)*ny + [j j+1];
   
    end
end

n = 0;

for k = 1:nx-1
    for j = 1:ny
      
        n = n + 1; 
        
        dx(n,:) = [(k-1)*ny+j k*ny+j]; 
   
    end
end

nbr = struct;
nbr.edx = dx;
nbr.edy = dy;




