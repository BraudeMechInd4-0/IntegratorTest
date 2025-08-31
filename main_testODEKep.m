% Test different ODE solvers

%% general parameters
% This section is required to compute r and v from the data downloaded from
% NORAD. If you got the r0 and v0, you can just input them manually and
% delete this section

fprintf('loading...');
if ~exist('consts','var')
    [xpdotp, Re, J2, mu, whichconsts, consts, ~, dAT] = generate_parameters;
end
tumin = consts(1);

%% these values are from the CATCH paper
deltas = [4,8,16,32]; %number of sections per orbit
Ns = [8 16 32]; %number of points per section
%% create a vector of all satellites

%GPdata = readstruct("ParseGP/last-30-days.xml");
%first three - LEO, skynet - GEO, ASBM - HEO
satnames = ["IRIDIUM33deb","STARLINK1341","QIANFAN4","SKYNET4C", "ASBM2"];
fprintf("done!\n starting satellite loops\n")
%% start loop over satellites
for satname = satnames
    for N = Ns
        time78i = zeros(length(deltas),1);
        time45i = zeros(length(deltas),1);
        time113i = zeros(length(deltas),1);
        time8i = zeros(length(deltas),1);
        time8 = zeros(length(deltas),1);
        time4i = zeros(length(deltas),1);
        time4 = zeros(length(deltas),1);
        timeMCPI = zeros(length(deltas),1);
        deltaCounter = 0;
        for delta = deltas
            deltaCounter = deltaCounter+1; %putting this here so I don't forget at the end
            fprintf(satname+"_N"+N+"_delta"+delta+"...")
            GPdata = readstruct("Satellites/"+satname+".xml");
            satstructxml = GPdata.omm.body.segment.data;

            satrec = GPxml2rv(whichconsts,consts,satstructxml);
            [~,r0,v0] = sgp4(satrec,0,consts);
            clear GPdata

            %% Prepare time point list

            [a,e,inc,O,w,f0] = kep_elements(r0,v0,mu); %transfer the r and v to
            %  keplerian elements
            Sec = 2*pi*sqrt((a^3)/mu)/delta; %These is the orbital period of the
            % satellite, devided by delta. It's the time of one section

            %% pick the time horizon for the integration
            %tmax = 60*60*24*365; %one year
            %tmax = 60*60*24*182; %half a year
            %tmax = 60*60*24*30*3; %three months
            %tmax = 60*60*24*30*2; %two months
            tmax = 60*60*24*30; %one month
            %tmax = 60*60*24*14; %two weeks
            %tmax = 60*60*24*7; %one week
            %tmax = 60*60*24; %one day
            %tmax = 60*60*12; % 12 hours
            %tmax = 60*60*6; % 6 hours
            %tmax = 60*60*3; % 3 hours
            %tmax = 60*60; % 1 hour
            %tmax = 60*30; % 1/2 hour

            Q = fix(tmax/Sec)+1; %number of sections we have
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %In case there are several sections we have N points per-orbit
            %In case we have a fraction of an orbit we still want to have N points
            %nevertheless
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if Q > 1

                omega2 =Sec/2; %omega 2 is the same for all sections because all
                % sections are the same size
                tau_ = cos(pi*(0:(N-1))/(N-1)); %from 1 to -1
                tau_ = tau_(end-1:-1:1); %flip the vector because it's from 1 to -1.
                % Also exclude -1 (because we have sections, the end of one is the
                % start of the other, I take out the first and in the end will add
                % the first point which is -1 mapped back to t = 0

                omega1 = Sec/2:Sec:tmax+Sec/2; %a vector of omega1 values one for
                % each section
                tau = repmat(tau_,1,length(omega1)); %replicate tau to the amount of
                % sections we have
                omega1 = repelem(omega1,length(tau_)); %replicate the values of omega1
                % so that you can add it to all the points
                tspan = tau*omega2+omega1;
                tspan = [0 tspan]; %we removed the zero when we flipped this is OK
                tmax = tspan(end);
            else
                %if we have just one section it's way easier
                tau_ = cos(pi*(0:(N-1))/(N-1)); %from 1 to -1
                tspan = tmax/2 + tmax/2*tau_;
                tspan = fliplr(tspan);
            end

            %% Generate baseline
            fprintf("generating baseline...")

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %If we're using the basic TBP as a test we can use Kepler's TBP solution as a
            %baseline. This is not quite an analitical solution:In propKep we call
            %KepEq which is a numerical soltuion. However, This solution is within a
            %general tolerance for each points, and unlike step-by-step integrators the         
            %errors don't accumulate
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            rbaseline = zeros(length(tspan),6);
            rbaseline(1,1:3) = r0;
            rbaseline(1,4:6) = v0;

            for i = 2:length(tspan)
                [rt,vt] = propKep(r0,v0,tspan(i),mu,1,1e-15,500); %propagate to the timepoint that is tspan(i) with 1e-15 tolerance
                if all(rt ~= 0) %converged. The error is there as a debug, should not happen now
                    rbaseline(i,:) = [rt vt];
                else
                    error("huh...")
                end
            end
            fprintf("done!\n")
            %% RUN native ODE TEST and fixed step - such that no need to check the interim
            %set initial step size and max step to T*1e-2 - default
            options = odeset('RelTol',1e-9,'AbsTol', 1e-12,'Stats','on','MaxStep',Sec/N);
            time_ = tic();
            [~,r113i] = ode113(@(t,x)orbit_eq(t,x,mu),tspan,[r0 v0],options);
            time113i(deltaCounter) = toc(time_)
            err113i = rbaseline - r113i;
            clear r113i
            time_ = tic();
            [~,r45i] = ode45(@(t,x)orbit_eq(t,x,mu),tspan,[r0 v0],options);
            time45i(deltaCounter) = toc(time_)
            err45i = rbaseline - r45i;
            clear r45i
            time_ = tic();
            [~,r78i] = ode78(@(t,x)orbit_eq(t,x,mu),tspan,[r0 v0],options);
            time78i(deltaCounter) = toc(time_)
            err78i = rbaseline - r78i;
            clear r78i

            %% RK8
            time_ = tic();
            [~,r8i]=ode8(@(t,x)orbit_eq(t,x,mu),tspan,[r0 v0]);
            time8i(deltaCounter)= toc(time_)

            err8i = rbaseline - r8i;
            clear r8i
            %% RK4
            time_ = tic();
            [~,r4i]=ode4(@(t,x)orbit_eq(t,x,mu),tspan,[r0 v0]);
            time4i(deltaCounter)= toc(time_)


            err4i = rbaseline - r4i;
            clear r4i

            %% Run constant step on normal ODE and estimate at the points using averaging with RK-8
            h=tmax/(length(tspan)); %same amount of points
            time_ = tic();
            [t8,r8,k81]=ode8(@(t,x)orbit_eq(t,x,mu),[0 tmax],[r0 v0], h); %compute the values at the fixed step
            time8_1 = toc(time_)
            time_ = tic();
            t8estimate = zeros(size(tspan));
            r8estimate = zeros(length(tspan),6);
            t8estimate(1) = t8(1);
            r8estimate(1,:) = r8(1,:);
            ct8 = 2;
            ctspan = 2;
            while ct8 <= length(t8estimate) && ctspan <= length(tspan)
                if tspan(ctspan) <= t8(ct8) && tspan(ctspan) > t8(ct8-1)
                    t8estimate(ctspan) = tspan(ctspan);
                    theta = (tspan(ctspan)-t8(ct8-1))/h;
                    r8estimate(ctspan,:) = hermiteInterp(r8(ct8-1,:),r8(ct8,:),k81(ct8-1,:),k81(ct8,:),theta,h);
                    ctspan = ctspan+1;
                else
                    ct8 = ct8+1;
                end
            end
            time8_2 = toc(time_)
            time8(deltaCounter) = time8_1+time8_2

            err8 = rbaseline - r8estimate;
            clear r8estimate t8 t8estimate r8 k81 k8end

            %% Run constant step on normal ODE and estimate at the points using averaging with RK-4
            h=tmax/(length(tspan)); %same amount of points
            time_ = tic();
            [t4,r4,k41]=ode4(@(t,x)orbit_eq(t,x,mu),[0 tmax],[r0 v0], h); %compute the values at the fixed step
            time4_1 = toc(time_)
            time_ = tic();
            t4estimate = zeros(size(tspan));
            r4estimate = zeros(length(tspan),6);
            t4estimate(1) = t4(1);
            r4estimate(1,:) = r4(1,:);
            ct4 = 2;
            ctspan = 2;
            while ct4 <= length(t4estimate) && ctspan <= length(tspan)
                if tspan(ctspan) <= t4(ct4) && tspan(ctspan) > t4(ct4-1)
                    t4estimate(ctspan) = tspan(ctspan);
                    theta = (tspan(ctspan)-t4(ct4-1))/h;
                    r4estimate(ctspan,:) = hermiteInterp(r4(ct4-1,:),r4(ct4,:),k41(ct4-1,:),k41(ct4,:),theta,h);
                    ctspan = ctspan+1;
                else
                    ct4 = ct4+1;
                end
            end
            time4_2 = toc(time_)
            time4(deltaCounter) = time4_1+time4_2

            err4 = rbaseline - r4estimate;
            clear r4estimate t4 t4estimate r4 k14 k44

            %% Run boundery problem MPCI in sections according to the times defined
            %The MPCI integration method is defined per-section. So we need to iterate
            %and run it for each section.
            fprintf("MPCI...")

            time_=tic;

            tstart  = 0;
            if Q > 1
                tend = Sec;
            else
                tend = tmax;
            end
            rMPCISeg = zeros(length(err78i),6);

            i = 1;
            rstart = [r0 v0];
            rtmp = [];

            while i+size(rtmp,1)-1 <= length(tspan)
                [~,rtmp,err] = odeMPCI(@(t,x)orbit_eq(t,x,mu),[tstart tend],rstart,1e-12,1e-9,N);

                if ~err
                    rMPCISeg(i:i+size(rtmp,1)-1,:) = rtmp;
                    i = i+size(rtmp,1)-1;
                    %update the next section:The new tstart is now tend and the new tend is
                    %the next section.
                    tstart = tend;
                    tend = tend + Sec;
                    rstart = rtmp(end,:);
                else
                    disp(i)
                    error("This should not happen")
                end
            end

            timeMCPI(deltaCounter) = toc(time_)
            fprintf("done!\n")

            errMPCISeg = rbaseline - rMPCISeg;

            %% plot results norm

            normr113i = sqrt(sum(err113i(:,1:3).^2,2))';
            normr45i = sqrt(sum(err45i(:,1:3).^2,2))';
            normr78i = sqrt(sum(err78i(:,1:3).^2,2))';
            normr8 = sqrt(sum(err8(:,1:3).^2,2))';
            normr8i = sqrt(sum(err8i(:,1:3).^2,2))';
            normr4 = sqrt(sum(err4(:,1:3).^2,2))';
            normr4i = sqrt(sum(err4i(:,1:3).^2,2))';
            normrMPCI = sqrt(sum(errMPCISeg(:,1:3).^2,2))';

            %% save results
            % dt = datestr(now,'yyyymmddHHMM');
            % save("Results/"+satname+"_TBP_"+delta+"_"+dt+".mat")


            %% plot the results, first the comparison of each method to itself in the error norm
            plot_tikz_figure(tspan,[normr8;normr8i],"Results/TBP_B/"+satname+"_TBP_delta"+delta+"_N"+N+"_RK8")
            plot_tikz_figure(tspan,[normr4;normr4i],"Results/TBP_B/"+satname+"_TBP_delta"+delta+"_N"+N+"_RK4")
            
            plot_tikz_figure(tspan,[normr4i;normr8i;normr45i;normr78i;normr113i;normrMPCI],"Results/TBP_B/"+satname+"_TBP_delta"+delta+"_N"+N)


        end %delta
        %% generate runtime table
        varnames={'delta','RK4','RK8','ODE45','ODE78','ODE113','MPCI'};
        T = table(deltas',time4i,time8i,time45i,time78i,time113i,timeMCPI,'VariableNames',varnames);
        table2latex(T,char("Results/TBP_B/"+satname+"_TBP_N"+N+"_runtime.tex"));

        varnames={'delta','RK4 Defined Step','RK4 Hermite Interpolation','RK8 Defined Step','RK8 Hermite Interpolation'};
        T = table(deltas',time4i,time4,time8i,time8,'VariableNames',varnames);
        table2latex(T,char("Results/TBP_B/"+satname+"_TBP_N"+N+"_RK_runtime.tex"));
    end%N
end


%% Play sound
fs = 8192;
toneduration = 0.1;
spaceduration = 0.05;
tonefreq = 800;
nbeeps = 10;

t = linspace(0,toneduration,round(toneduration*fs));
y=0.8*sin(2*pi*tonefreq*t);
ys = zeros(1,round(spaceduration*fs));

Y=[repmat([y ys],[1 nbeeps-1]) y];
sound(Y,fs);


