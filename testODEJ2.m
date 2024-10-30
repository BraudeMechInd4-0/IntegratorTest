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
%% parameters for J2 and drag
Cd = 2.2;
As = [3.9,0.7,4,10,12]*1e-6;%km^2
m_s = [260,77,260,1250,2000];%kg

%% these values are from the CATCH paper
deltas = [4,8,16,32]; %number of sections per orbit
Ns = [8 16 32]; %number of points per section
%% create a vector of all satellites

%GPdata = readstruct("ParseGP/last-30-days.xml");
satnames = ["STARLINK1341","IRIDIUM33deb","QIANFAN4","SKYNET4C", "ASBM2"];
fprintf("done!\n starting satellite loops\n")
%% start loop over satellites
for i = 1:length(satnames)
    satname = satnames(i);
    A = As(i);
    m = m_s(i);
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
            deltaCounter = deltaCounter+1;
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
            tmaxs = 60*60*24*[7,14,30]; %one week, two weeks and a month
            tmax = tmaxs(end);

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
            %This is the baselone for a model including drag and J2 - very accurate integration that we'll compare to
            %The step size is determined to be such that will give about 10 times more
            %points-per-section
            %
            % @(t,x)orbit_eq_J2_drag is what one may define as f. It is the function to be
            % integrated
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            options = odeset('RelTol',1e-14,'AbsTol', 1e-15,'MaxStep',Sec/(N*10));
            fprintf("making baseline...")
            [~,rbaseline] = ode89(@(t,x)orbit_eq_J2_drag(t,x,mu,Cd,A,m,Re,J2),tspan,[r0 v0],options);


            fprintf("done!\n")
            %% RUN native ODE TEST and fixed step - such that no need to check the interim
            %set initial step size and max step to T*1e-2 - default
            options = odeset('RelTol',1e-9,'AbsTol', 1e-12,'Stats','on','MaxStep',Sec/N);
            time_ = tic();
            [~,r113i] = ode113(@(t,x)orbit_eq_J2_drag(t,x,mu,Cd,A,m,Re,J2),tspan,[r0 v0],options);
            time113i(deltaCounter) = toc(time_)
            err113i = rbaseline - r113i;
            clear r113i
            time_ = tic();
            [~,r45i] = ode45(@(t,x)orbit_eq_J2_drag(t,x,mu,Cd,A,m,Re,J2),tspan,[r0 v0],options);
            time45i(deltaCounter) = toc(time_)
            err45i = rbaseline - r45i;
            clear r45i
            time_ = tic();
            [~,r78i] = ode78(@(t,x)orbit_eq_J2_drag(t,x,mu,Cd,A,m,Re,J2),tspan,[r0 v0],options);
            time78i(deltaCounter) = toc(time_)
            err78i = rbaseline - r78i;
            clear r78i

            %% RK8
            time_ = tic();
            r8i = NaN(size(r78i));
            try %catch the error if the satellite decays
                [~,r8i]=ode8(@(t,x)orbit_eq_J2_drag(t,x,mu,Cd,A,m,Re,J2),tspan,[r0 v0]);
                time8i(deltaCounter)= toc(time_)
            catch
                warning("Error too big? r8i will e NaN")
                time8i(deltaCounter)= NaN
            end

            err8i = rbaseline - r8i;
            clear r8i
            %% RK4
            time_ = tic();
            r4i = NaN(size(r78i));
            try
                [~,r4i]=ode4(@(t,x)orbit_eq_J2_drag(t,x,mu,Cd,A,m,Re,J2),tspan,[r0 v0]);
                time4i(deltaCounter)= toc(time_)
            catch
                warning("Error too big? r8i will e NaN")
                time4i(deltaCounter) = NaN;
            end


            err4i = rbaseline - r4i;
            clear r4i

            %% Run constant step on normal ODE and estimate at the points using averaging with RK-8
            h=tmax/(length(tspan)); %same amount of points
            if h > Sec/N
                h = Sec/N;
            end
            time_ = tic();
            r8 = NaN(size(r78i));
            try
                [t8,r8,k81]=ode8(@(t,x)orbit_eq_J2_drag(t,x,mu,Cd,A,m,Re,J2),tspan,[r0 v0]);
                time8(deltaCounter)= toc(time_)
            catch
                warning("Error too big? r8i will e NaN")
                time8(deltaCounter)= NaN
            end
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

            %% Run constant step on normal ODE and estimate at the points using hermite interpolation with RK-4
            h=tmax/(length(tspan)); %same amount of points
            if h > Sec/N
                h = Sec/N;
            end
            time_ = tic();
            try
                [t4,r4,k41]=ode8(@(t,x)orbit_eq_J2_drag(t,x,mu,Cd,A,m,Re,J2),tspan,[r0 v0]);
                time4(deltaCounter)= toc(time_)
            catch
                warning("Error too big? r8i will e NaN")
                time4(deltaCounter)= NaN
            end
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

           
            %% plot results norm

            normr113i = sqrt(sum(err113i(:,1:3).^2,2))';           
            normr45i = sqrt(sum(err45i(:,1:3).^2,2))';
            normr78i = sqrt(sum(err78i(:,1:3).^2,2))';
            normr8 = sqrt(sum(err8(:,1:3).^2,2))';
            normr8i = sqrt(sum(err8i(:,1:3).^2,2))';
            normr4 = sqrt(sum(err4(:,1:3).^2,2))';
            normr4i = sqrt(sum(err4i(:,1:3).^2,2))';

            %% plot the results, first the comparison of each method to itself in the error norm
            plot_tikz_figure(tspan,[normr8;normr8i],"Results/J2B/"+satname+"_J2_delta"+delta+"_N"+N+"_RK8")
            plot_tikz_figure(tspan,[normr4;normr4i],"Results/J2B/"+satname+"_J2_delta"+delta+"_N"+N+"_RK4")

            %% run MPCI tests
            normrMPCIinit = zeros(length(rbaseline),length(tmaxs));
            normrMPCI = zeros(length(rbaseline),length(tmaxs));
            timeMCPIinit = zeros(1,length(tmaxs));
            timeMCPI = zeros(1,length(tmaxs));

            for k = 1:length(tmaxs)
                tmax = tmaxs(k);

                %% Run boundery problem MPCI in sections according to the times defined
                %The MPCI integration method is defined per-section. So we need to iterate
                %and run it for each section.
                fprintf("MPCI with initial guess...")

                time_=tic;

                tstart  = 0;
                if Q > 1
                    tend = Sec;
                else
                    tend = tmax;
                end
                rMPCISeg = zeros(size(rbaseline));

                j = 1;
                rstart = [r0 v0];
                rtmp = [];

                while j+size(rtmp,1)-1 <= length(tspan)
                    % The initial guess is only profitable when propagating in small time
                    % frames. The longer we propagate the further we drift from the
                    % originial TBP solution and then we get garbage in. It can be
                    % interesting to check when exactly we lose the TBP guess advantage.
                    %
                    % A solution is to have each section start where the previous ended and
                    % accumulate the "errors" on purpose. This makes for slightly shorter
                    % runtime? Need to check that as well.
                    %
                    %can we run this section 2 times: one time with no initial guess, and
                    %one time with the initial guess being the previous section's end?
                    [rinit,vinit] = initialGuess(rstart(1:3),rstart(4:6),Sec,N,mu);
                    [~,rtmp,err] = odeMPCI(@(t,x)orbit_eq_J2_drag(t,x,mu,Cd,A,m,Re,J2),[tstart tend],rstart,1e-12,1e-9,N,[rinit vinit]);

                    if ~err
                        rMPCISeg(j:j+size(rtmp,1)-1,:) = rtmp;
                        j = j+size(rtmp,1)-1;
                        %update the next section:The new tstart is now tend and the new tend is
                        %the next section.
                        tstart = tend;
                        tend = tend + Sec;
                        rstart = rtmp(end,:);
                    else
                        disp(j)
                        error("This should not happen")
                    end
                end

                %profile viewer
                %profsave
                timeMCPIinit(k) = toc(time_)
                fprintf("done!\n")

                errMPCISeginit = rbaseline - rMPCISeg;

                normrMPCIinit(:,k) = sqrt(sum(errMPCISeginit(:,1:3).^2,2));
                %normMPCI = sqrt(sum(errMPCISeg(:,4:6).^2,2))';


                %% Run boundery problem MPCI in sections according to the times defined
                %The MPCI integration method is defined per-section. So we need to iterate
                %and run it for each section.

                % The initial guess is only profitable when propagating in small time
                % frames. The longer we propagate the further we drift from the
                % originial TBP solution and then we get garbage in. It is
                % interesting to check when exactly we lose the TBP guess advantage.
                fprintf("MPCI without initial guess...")

                time_=tic;

                tstart  = 0;
                if Q > 1
                    tend = Sec;
                else
                    tend = tmax;
                end
                rMPCISeg = zeros(size(rbaseline));

                j = 1;
                rstart = [r0 v0];
                rtmp = [];

                % we run this section 2 times: one time with no initial guess, and
                %one time with the initial guess being the previous section's end
                %[rinit,vinit] = initialGuess(rstart(1:3),rstart(4:6),Sec,N,mu);

                while j+size(rtmp,1)-1 <= length(tspan)

                    % A solution is to have each section start where the previous ended and
                    % accumulate the "errors" on purpose. This makes for slightly shorter
                    % runtime? Need to check that as well.
                    %

                    [~,rtmp,err] = odeMPCI(@(t,x)orbit_eq_J2_drag(t,x,mu,Cd,A,m,Re,J2),[tstart tend],rstart,1e-12,1e-9,N);

                    if ~err
                        rMPCISeg(j:j+size(rtmp,1)-1,:) = rtmp;
                        j = j+size(rtmp,1)-1;
                        %update the next section:The new tstart is now tend and the new tend is
                        %the next section.
                        tstart = tend;
                        tend = tend + Sec;
                        rstart = rtmp(end,:);
                    else
                        disp(j)
                        error("This should not happen")
                    end
                end

                %profile viewer
                %profsave
                timeMCPI(k) = toc(time_)
                fprintf("done!\n")

                errMPCISeg = rbaseline - rMPCISeg;


                normrMPCI(:,k) = sqrt(sum(errMPCISeginit(:,1:3).^2,2))';
                %normMPCI = sqrt(sum(errMPCISeg(:,4:6).^2,2))';
            end %tmax
            %% generate runtime table
            %days = tmax/60/60/24;
            daylist = round(tmaxs/60/60/24)';
            varnames = {'Days','With Initial Guess', 'Without Initial guess'};
            T = table(daylist,timeMCPIinit',timeMCPI','VariableNames',varnames);
            %writetable(T,"Results/J2n/"+satname+"_J2_"+daytime+"_"+delta+"_MPCIruntime.csv");
            table2latex(T,char("Results/J2B/"+satname+"_J2_delta"+delta+"_N"+N+"_MPCI_runtime.tex"));

            %% save results
            % dt = datestr(now,'yyyymmddHHMM');
            % days = tmax/60/60/24;/home/eladden/Downloads/table2latex.m
            % save("Results/J2n/"+satname+"_J2_"+daytime+"_"+delta+"_"+dt+".mat")


            plot_tikz_figure(tspan,[normr4i;normr8i;normr45i;normr78i;normr113i;normrMPCIinit(:,end)';normrMPCI(:,end)'],"Results/J2B//"+satname+"_J2_delta"+delta+"_N"+N)


        end%delta
        varnames={'delta','RK4','RK8','ODE45','ODE78','ODE113'};
        T = table(deltas',time4i,time8i,time45i,time78i,time113i,'VariableNames',varnames);
        table2latex(T,char("Results/J2B/"+satname+"_J2_N"+N+"_runtime.tex"));

        varnames={'delta','RK4 Defined Step','RK4 Hermite Interpolation','RK8 Defined Step','RK8 Hermite Interpolation'};
        T = table(deltas',time4i,time4,time8i,time8,'VariableNames',varnames);
        table2latex(T,char("Results/J2B/"+satname+"_J2_N"+N+"_RK_runtime.tex"));
    end% N
end % satnames


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


