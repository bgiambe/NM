function [vect_or,d_or,varb,varb_eig, Zwork, Nsq22,Dwork] = n_modes(nn,Nsq2,Z,D,bctype,g)
% depth must be negative, bctypes:
%                                   1.   dp/dz = 0 at z=-H ,Free surface
%                                   2.   dp/dz = 0 at z=-H ,% Rigid surface (dp/dz = 0)
%                                   3.   p=0 at z = -H, % Rigid surface for W solution (W=0)

    % decimate the data
    m = floor(-Z(length(Z))/nn); 
    m(m==0)=1; % prevent m=0 

    Dwork = D(1:m:length(D));
    alpha = ones(length(Dwork),1)./Dwork;		% specific volume
    Zwork = Z(1:m:length(Z));			
    Nsq22 = Nsq2(1:m:length(Nsq2));	


    % Setting up the coefficient matrix

    n = length(alpha)-1;
    M = zeros(n);
    C = 1;

    for k = 2:(n-1)
    %Processing row k of a total of n-1
        	kp1 = k+1;
            km1 = k-1;
            aux = (Zwork(kp1)-Zwork(km1))/2; 
            A = 1/(Nsq22(k)*(aux^2));
            B = (1/((4*(aux^2))*Nsq22(k)))*((alpha(kp1)-alpha(km1))/alpha(k) - ...
                (Nsq22(kp1)-Nsq22(km1))/Nsq22(k));
            M(k,(k-1):(k+1)) = [A-B, -2*A, A+B];
    end
	
    %Applying boundary conditions
    % Bottom BC:
    aux2 = (Zwork(n+1)-Zwork(n-1))/2;		
    A = 1/(Nsq22(n)*(aux2^2));
    B = (1/((4*(aux2^2))*Nsq22(n)))*((alpha(n+1)-alpha(n-1))/alpha(n) - ...
                (Nsq22(n+1)-Nsq22(n-1))/Nsq22(n));

    if (bctype == 1 || bctype == 2)
        M(n,(n-1):n) = [A-B , B-A];	% dp/dz = 0 at z=-H 
    else
        M(n,(n-1):n) = [A-B , -2*A];	% W=0 at z = -H
    end

    % Top BC:
    aux = 1/(Zwork(2)-Zwork(1));			
    aux2 = (Zwork(3)-Zwork(1))/2; 
    A = 1/(Nsq22(2)*(aux2^2));
    B = (1/(4*(aux2^2)*Nsq22(2)))*((alpha(3)-alpha(1))/alpha(2) - ...
                (Nsq22(3)-Nsq22(1))/Nsq22(2));
    C = aux/((Nsq22(1)/g) - aux);

    if bctype == 1
        M(1,1:2) = [(-A*(2+C)+(C*B)) , A+B]; % Free surface
    elseif bctype == 2
		M(1,1:2) = [-A-B , A+B]; % Rigid surface (dp/dz = 0)
    else
		M(1,1:2) = [-2*A , A+B]; % Rigid surface for W solution (W=0)
    end

    %Computing the eigenvalues 
    [vect,d] = eig(M);  
    d = diag(d);
    
    
    %Order the eigenvalues and keep track of the
    %corresponding eigenvectors

    d = -d;			% Because of the negative sign on the RHS 
    [d_or, idb]= sort(d);

    vect_or=vect(:,idb);
    
    %Compute percentage of fit for each mode

	dotprod = zeros(n,1);
    for k=1:n
		dotprod(k) = (Dwork(1:n)')*vect(:,idb(k)); % This vector will also be ordered
    end
    
	soma = sum(abs(dotprod));
	varb = abs(dotprod)*100/soma; 

    %Compute energy content for each mode
   
    soma = sum(abs(d(idb).^2));
    varb_eig = abs(d(idb).^2)*100/soma; 
    varb_eig = sort(varb_eig,'descend');


    
    
end
    
