function [rinit,vinit] = initialGuess(r0,v0,tend,N,mu)
%[xinit,vinit] = initialGuess(r0,v0,tstart,tend,N)
%generate an initial guess for the MPCI in an astrodynamics environment

tau = fliplr(cos((0:(N-1))*(pi/(N-1)))); %Gauss lobato nodes
T = (tau+1)*(tend)/2; %These are the timepoints of interest over the interval

rinit = zeros(length(T),3);
vinit = zeros(length(T),3);

for i = 1:length(T)
    [r,v] = propKep(r0,v0,T(i),mu,2,1e-13,500);
    if all(r==0)  %sometimes the kepler equation does not converge
        if i ~=1
            rinit(i,:) = rinit(i-1,:);
            vinit(i,:) = vinit(i-1,:);
        else
            err = ~(all(r==0));
            j = 1;
            while err
                [r,v] = kepler(r0,v0,T(i)+j,25);
                err = ~(all(r==0));
                j=j+1;
            end
            rinit(i,:) = r;
            vinit(i,:) = v;
        end
    else %no problem
        rinit(i,:) = r;
        vinit(i,:) = v;
    end
end


end