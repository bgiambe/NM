function [Nsq2,D,Z,Nsq22, Zwork,Dwork,vect_or,d,varb,varb_eig,hn,Lr,cn] = NM_fun(data, lat,fil,opt,nn,bctype)
%normal modes calculation for CTD cast (data points every 1 m)
%-----//INPUT
% data= table with VariableNames = {'Depth', 'Temperature', 'Salinity'}
% with units: m, °C, PSU respectively.
% lat = latitude in decimal degrees.
% ---filtering options:
% fil= 
%               'SG1' Savitsky-Golay smoothing filter opt=[nl,nr, order]
%                    NL and NR are the left and 
%                    right span in the convolution, respectively, and ORDER is the 
%                    order of the filter <=> statistical moment to be preserved in
%                    the filtered output.
%               'SG2' Savitzky-Golay (FIR) smoothing filter opt=[fl, order]
%                    FL frame length to the data (odd), ORDER is the 
%                    order of the filter
%               'IIR1' zero-pole lowpass filter of order 30, opt= [pass, cut]
%                     PASS is passband and CUT is stopband between [0,1]
%                     where 1 is the Nyquist wavenumber
%               'IIR2':zero-pole lowpass filter with adjustable order
%                      opt=[pass, cut, Rp, Rs] Rp dB passband ripple, 
%                      Rs dB  attenuation in the stopband,
%                      PASS is passband and CUT is stopband between [0,1]
%                      where 1 is the Nyquist wavenumber
%               'FIR1' all-zeros lowpass filter (moving average) of order 30,
%                      opt = [pass, cut] same as IIR
%               'FIR2':all-zeros lowpass filter (moving average) opt=[order, pass, cut]
%                      ORDER is the order of the filter, PASS is passband
%                      and CUT is stopband between [0,1] where 1 is the Nyquist wavenumber
%               'B':   Butterworth filter opt= [order, cut] ORDER is the order
%                      of the filter, CUT is stopband between [0,1] where 1 is the Nyquist wavenumber
% ---normal modes calculation options:
% nn = Decimate the data for every mod(ZMAX,nn) point
% bctypes = boundary conditions:
%                                   1.   dp/dz = 0 at z=-H ,Free surface
%                                   2.   dp/dz = 0 at z=-H ,% Rigid surface (dp/dz = 0)
%                                   3.   p=0 at z = -H, % Rigid surface for W solution (W=0)
%                                           --> in this case the first mode
%                                           is not the Barotropic, but the
%                                           first baroclinic mode
%                                          
%-----//OUTPUT
%Nsq2 = filtered N squared profile
%D = density profile
%Z = depth
%Nsq22 = decimated N squared profile
%Zwork = decimated depth
%Dwork = decimated density profile
%vect_or = ordered eigenvector (1 BT, 2 first BC,...)
%d = ordered eigenvalues
%varb=fit percentages
%varb_eig = ordered mode energies
%hn = equivalent depths
%Lr = internal radius
%cn = velocity


% set variables and remove missing entries

[T, tf]=rmmissing(data.Temperature);

if any(data.Depth<0)
    P=-data.Depth;
else
    P=data.Depth;
end
S=data.Salinity;

P(tf) = [];
S(tf) = [];

% nsq calculation
rr=length(P); 
p_ref = median(P); 

try
    [Nsq,D, g, Z,~,~] = bruntvais(S,T, P, lat, p_ref,rr);
catch
    [Nsq,D, g, Z,~,~] = bruntvais(S,T, P', lat, p_ref,rr);
end


% filtering 
Nsq2 = filtering(Nsq,fil,opt);

% modes calculation
[vect_or,d,varb,varb_eig, Zwork,Nsq22,Dwork] = n_modes(nn,Nsq2,Z,D,bctype,g);


% internal radius
f= gsw_f(lat); % radians/s % coriolis parameter
n = 1:1:length(Zwork);
hn = abs(1./(g.*real(d))) ;
N=real(sqrt(Nsq22));
Lr = zeros(length(N)-1,1);
for i=1:length(N)-1
    Lr(i) = N(i)*hn(i)/(n(i)*pi*f);  % (m)
end

%velocity
csq = abs(1./real(d));
cn = real(sqrt(csq)); %m/s


end