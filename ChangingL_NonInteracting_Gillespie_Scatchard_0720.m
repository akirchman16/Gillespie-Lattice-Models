clear all;
close all;

% This model will bind and unbind proteins at times according to the
% Gillespie Stochastic Simulation Algorithm (SSA). In order to make sure
% the model is accurate it is compared to theoretical results in a
% Scatchard plot from McGhee paper. The binding/unbinding reaction is
% modeled as: A+B<-->AB. In this model the exact populations of each
% species in the reaction is tracked, none are held constant. A variety of
% other plots can also be created.

N = 1000;    %length of DNA lattice
n = 10;  %length of each protein
k_on = 0.01;   %kinetic rate constant for binding
k_off = 0.01;  %kinetic rate constant for unbinding
Vol = 100;    %volume of system in which reaction occurs
maxPop = 1000;  %maximum population of free proteins for a run

Iterations = 500;  %number of events which will occur

K = k_on/k_off; %equilibrium constant
c_f = k_on/Vol; %constants for propensity functions
c_r = k_off;

L = zeros(1,ceil(N/n));   %memory allocation
v = zeros(1,ceil(N/n));
ScatchY = zeros(1,ceil(N/n));
BindingFrac = zeros(1,ceil(N/n));
UnbindingFrac = zeros(1,ceil(N/n));
AvgTime = zeros(1,ceil(N/n));
EqCover = zeros(1,ceil(N/n));
TotalTime = zeros(1,ceil(N/n));
Pop1 = zeros(1,ceil(N/n));

Loops = 1;

for p = ceil(N/n):maxPop
    Pop1(Loops) = p; %array to track number of runs
    
    DNA = zeros(1,N);   %empty DNA lattice
    BoundAtSpot = zeros(1,N);  %empty CurrentBound array

    xA = zeros(1,Iterations+1);   %memory allocation arrays
    xB = zeros(1,Iterations+1);
    xAB = zeros(1,Iterations+1);
    t = zeros(1,Iterations+1);
    a_f = zeros(1,Iterations);
    a_r = zeros(1,Iterations);
    a_0 = zeros(1,Iterations);
    dt = zeros(1,Iterations);
    Hist = zeros(2,Iterations); %top row is binding event, bottom row is unbinding events
    ProteinCount = zeros(1,Iterations+1);
    FracCover = zeros(1,Iterations+1);

    t(1) = 0;   %intial time is zero
    BindCounter = 0;    %initially counters are at zero
    UnbindCounter = 0;

    xA(1) = p; %initial populations for an empty lattice
    xB(1) = N-(n-1);
    xAB(1) = 0;

    ProteinCount(1) = xAB(1);   %how many proteins are currently bound; 0 for initially empty lattice
    FracCover(1) = (xAB(1)*n)/N; %initial fractional coverage; 0 for an initially empty lattice

    for i = 1:Iterations
        a_f(i) = c_f*xA(i)*xB(i);   %propensity functions
        a_r(i) = c_r*xAB(i);
        a_0(i) = a_f(i)+a_r(i);

        dt(i) = (1/a_0(i))*log(1/rand); %random time interval for Gillespie method

        if a_f(i) > rand*a_0(i) %a forward reaction occurs
            eventB = 0;
            while ~eventB   %repeats until a binding occurs
                SpotB = randi(N-(n-1));  %random location on lattice is chosen
                if DNA(SpotB:SpotB+(n-1)) == 0   %checks if location is free
                    DNA(SpotB:SpotB+(n-1)) = 1;   %binds protein to location
                    BoundAtSpot(SpotB) = 1;  %stores locatin in BoundAtSpot
                    BindCounter = BindCounter+1;    %updates bind counter
                    Hist(1,i) = SpotB;      %stores location of binding event

                    xA(i+1) = xA(i)-1;    %population updates for a binding event
                    FreeSpots = 0;          %process for calculating how many free spots remain on the lattice
                    for j = 1:N-(n-1)
                        if DNA(j:j+(n-1)) == 0
                            FreeSpots = FreeSpots+1;
                        end
                    end
                    xB(i+1) = FreeSpots;
                    xAB(i+1) = xAB(i)+1;

                    eventB = 1;
                end
            end
        else                       %otherwise an unbinding has to occur
            CurrentBound = find(BoundAtSpot == 1);
            pos = randi(length(CurrentBound));
            SpotU = CurrentBound(pos);  %random position for protein to unbind

            DNA(SpotU:SpotU+(n-1)) = 0; %unbinds protein
            BoundAtSpot(SpotU) = 0; %removes location from BoundAtSpot
            UnbindCounter = UnbindCounter+1;    %updates unbind counter
            Hist(2,i) = SpotU;  %stores location of unbinding event

            xA(i+1) = xA(i)+1;    %population updates for an unbinding event
            FreeSpots = 0;          %process for calculating how many free spots remain on the lattice
            for j = 1:N-(n-1)
                if DNA(j:j+(n-1)) == 0
                    FreeSpots = FreeSpots+1;
                end
            end
            xB(i+1) = FreeSpots;
            xAB(i+1) = xAB(i)-1;
        end
        t(i+1) = t(i)+dt(i);    %advances time
        ProteinCount(i+1) = sum(DNA)/n;   %all of these values should be integers
        FracCover(i+1) = sum(DNA)/N;
        
        L_loop(i) = xA(i+1)/Vol;    %free protein concentration after each loop
        v_loop(i) = (sum(DNA)/n)/N; %binding density after each loop
    end
    L(Loops) = L_loop(Iterations);     %equilibrium free protein concentration
    v(Loops) = v_loop(Iterations);      %equilibrium binding density
    ScatchY(Loops) = v(Loops)/L(Loops); %y-value for Scatchard plot
    
    BindingFrac(Loops) = BindCounter/Iterations;   %fraction of events which were binding for each loop
    UnbindingFrac(Loops) = UnbindCounter/Iterations;   %fraction of events which were unbinding for each loop
    AvgTime(Loops) = mean(dt);     %average time between events in each loop
    EqCover(Loops) = mean(FracCover(round(0.75*Iterations):Iterations));    %equilibrium fractional coverage for each loop
    TotalTime(Loops) = max(t);  %total time for each loop
    
    Loops = Loops+1;
end

x = 0:(1/n)/1000:1/n;    %x-values for theoretical Scatchard plot
TheorScatchY = K.*(1-(n.*x)).*(((1-(n.*x))./(1-((n-1).*x))).^(n-1));   %Eq. 10 in McGhee paper

figure();   
subplot(2,2,1); %plot fractional coverage over time (last run only)
scatter(t,FracCover,5,'r','filled');
xlabel('Time, t');
xlim([0 max(t)]);
ylabel('Fractional Coverage');
ylim([0 1]);
title('Fractional Coverage (Final Run)');

subplot(2,2,2);     %histogram of time intervals (only for last run)
histogram(dt,100);
xlabel('Time Interval, dt');
title('dt Histogram (Final Run)');

subplot(2,2,3);     %Scatchard plot benchmark (all runs)
scatter(v,ScatchY,5,'r','filled');
hold on;
plot(x,TheorScatchY,'black');
xlabel('v');
xlim([min(v)-0.1 max(v)+0.1]);
ylabel('v/L');
ylim([0 max(ScatchY)+0.1]);
title('Scatchard Plot (K=1)');

subplot(2,2,4);
scatter(t(1:Iterations),v_loop,5,'g','filled');   %plot of binding density over time (final run)
xlabel('Time, t');
xlim([0 max(t)]);
ylabel('Binding Density');
title('Binding Density (Final Run)');

figure();
subplot(2,2,1);
scatter(Pop1,BindingFrac,5,'r','filled');    %plots percentage of events which were binding events
hold on;
scatter(Pop1,UnbindingFrac,5,'g','filled');  %plots percentage of events which were unbinding events
xlabel('Initial Free Pop.');
xlim([ceil(N/n) maxPop]);
ylabel('Percentage');
ylim([0 1]);
title('Fraction of Events');
legend('Binding','Unbinding');

subplot(2,2,2);
scatter(Pop1,EqCover,5,'b','filled');   %plots equilibrium fractional coverage for each run
xlabel('Initial Free Pop.');
xlim([ceil(N/n) maxPop]);
ylabel('Equilibrium Cover');
ylim([0 1]);
title('Equilibrium Fractional Coverage');

subplot(2,2,3);
scatter(Pop1,AvgTime,5,'r','filled');    %plots average time step for each run
xlabel('Initial Free Pop.');
xlim([ceil(N/n) maxPop]);
ylabel('Avg Time Step');
title('Average Time Step');

subplot(2,2,4);
scatter(Pop1,TotalTime,5,'b','filled'); %plots total time for the number of iterations for each loop
xlabel('Initial Free Pop.');
xlim([ceil(N/n) maxPop]);
ylabel('Total Time');
title('Total Time');

% fit data (ScatchY vs v) with this equation: K*(1-(n*v))*(((1-(n*v))/(1-((n-1)*v)))^(n-1))