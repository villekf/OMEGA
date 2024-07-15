function theta = Gauss_Newton(z,h_i,J_h,theta_i,k)


theta=theta_i+k*(J_h'*J_h)\J_h'*(z-h_i);
    
     
     