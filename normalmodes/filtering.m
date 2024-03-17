function x_filt = filtering(x,fil,opt)
%    fil= 
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

        if strcmp(fil, 'SG1')
        if length(opt)<3||length(opt)>3
            error(' You must supply 3 input arguments')
        end
        x_filt = savgol(x,opt(1),opt(2),opt(3));
    end
    
    if strcmp(fil, 'SG2')
        if length(opt)<2||length(opt)>2
            error(' You must supply 2 input arguments')
        end
        if mod(opt(1),2)==0
            error('Framelength must be odd')
        end
        x_filt = sgolayfilt(x,opt(2),opt(1));
    end
    
    if strcmp(fil, 'IIR1')
        if length(opt)<2||length(opt)>2
            error(' You must supply 2 input arguments')
        end
        if opt(1)<0||opt(1)>1 || opt(2)<0||opt(2)>1
            error('passband and stopband must be in [0,1]')
        end
        [N,Wn] = ellipord(opt(1),opt(2),.3,30);
        [B,A] = ellip(N,.3,30,Wn);
        x_filt = filter(B,A,x); 
        x_filt = [x_filt((N+1):length(x_filt)); x_filt(length(x_filt))*ones(N,1)];
    end
    
    if strcmp(fil, 'IIR2')
        if length(opt)<2||length(opt)>2
            error(' You must supply 2 input arguments')
        end
        if opt(1)<0||opt(1)>1 || opt(2)<0||opt(2)>1
            error('passband and stopband must be in [0,1]')
        end

        [N,Wn] = ellipord(opt(1),opt(2),opt(3),opt(4)); 
        [B,A] = ellip(N,Rp,Rs,Wn);
        x_filt = filter(B,A,x); 

    end
   
    if strcmp(fil, 'FIR1')
        if length(opt)<2||length(opt)>2
            error(' You must supply 2 input arguments')
        end
        if opt(1)<0||opt(1)>1 || opt(2)<0||opt(2)>1
            error('passband and stopband must be in [0,1]')
        end
        B = fir1(30,opt(2));
        nfilt = 31;
        itemp = 2*x(1)-x((nfilt+1):-1:2);
        [itemp,zi] = filter(B,1,itemp);
        [odata,zf]=filter(B,1,x,zi);
        itemp = zeros(2*nfilt,1);
        nd = max(size(x));
        itemp(:)=[2*x(nd)-x((nd-1):-1:(nd-2*nfilt))];
        x_filt = [odata;filter(B,1,itemp,zf)];
        x_filt(1:15) = [];x_filt=x_filt(1:nd);
    end
    
    if strcmp(fil, 'FIR2')
        if length(opt)<2||length(opt)>2
            error(' You must supply 2 input arguments')
        end
        if opt(2)<0||opt(2)>1 || opt(3)<0||opt(3)>1
            error('passband and stopband must be in [0,1]')
        end

        B = fir1(opt(1),[opt(2),opt(3)]);

        x_filt = filter(B,1,x);
     end


   if strcmp(fil, 'B')
        if length(opt)<2||length(opt)>2
            error(' You must supply 2 input arguments')
        end
        if opt(2)<0||opt(2)>1 
            error('stopband must be in [0,1]')
        end

        [b,a] = butter(opt(1),opt(2));
        temp = x(1)*ones(length(x),1);
        temp2=filter(b,a,[temp;x]);
        x_filt =temp2(length(temp)+1:end);
   end
end