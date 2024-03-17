function [Nsq,D, g, Z,S,rhobar] = bruntvais(S,T, P, lat, p_ref,rr)
% P in dbar, S in PSU, T in °C

    g = gsw_grav(lat,p_ref);
    Z = -P;
    

    % density

    D = density(S,T,P./10); % (kg/m^3) %needs P in bar
   

    % Nsq
    k=length(D);
    
    rhobar = median(D(1:rr));
     
    cw = svel(S,T,P); % (m/sec)
    cw = (cw(2:k)+cw(1:(k-1)))/2;			% Average
    
	Nsq = (-g*(diff(D)./diff(Z))./rhobar);	% gradient
	Nsq = Nsq - (g*ones(size(cw))./cw).^2;	% Compress. correction for
							% in situ gradient
	Nsq = [Nsq(1);Nsq];

end