clear all;
close all;

% This model binds and unbinds proteins in time intervals according to the
% Gillespie algorithm (Stochastic Simulation Algorithm). This is a dynamic
% Monte Carlo model.Reactions are modeled as a reversible, bimolecular
% reaction: A+B<-->AB, where A is a free protein, B is a free location on a
% DNA lattice, and AB is a bound protein on the lattice. Results of the
% model are used to compare to a benchmark of a Scatchard plot with 
% theoretical values provided by McGhee and von Hippel (1974). Propensity
% functions are determined based on chemical kinetics. No cooperativity is
% included in this version of the code.

N = 1000;    %length of DNA lattice
n = 3;  %length of each protein
k_on = 1;   %kinetic rate constant for binding
k_off = 1;  %kinetic rate constant for unbinding

ProteinConc = [0.1:0.1:10];  %various concentrations of free protein (Molarity)
minIterations = 5000;  %number of events which will occur

K = k_on/k_off; %equilibrium constant

L = zeros(1,length(ProteinConc));   %memory allocation
v = zeros(1,length(ProteinConc));
ScatchY = zeros(1,length(ProteinConc));
BindingFrac = zeros(1,length(ProteinConc));
UnbindingFrac = zeros(1,length(ProteinConc));
AvgTime = zeros(1,length(ProteinConc));
EqCover = zeros(1,length(ProteinConc));
TotalTime = zeros(1,length(ProteinConc));

Loops = 1;

for m = ProteinConc
    A_Conc(Loops) = m;    %concentration of free protein for each loop
    
    DNA = zeros(1,N);   %empty DNA lattice
    BoundAtSpot = zeros(1,N);  %empty CurrentBound array

    L_A = zeros(1,minIterations+1);   %memory allocation arrays
    xB = zeros(1,minIterations+1);
    xAB = zeros(1,minIterations+1);
    t = zeros(1,minIterations+1);
    a_f = zeros(1,minIterations);
    a_r = zeros(1,minIterations);
    a_0 = zeros(1,minIterations);
    dt = zeros(1,minIterations);
    Hist = zeros(2,minIterations); %top row is binding event, bottom row is unbinding events
    ProteinCount = zeros(1,minIterations+1);
    FracCover = zeros(1,minIterations+1);

    t(1) = 0;   %intial time is zero
    BindCounter = 0;    %initially counters are at zero
    UnbindCounter = 0;

    L_A = m;    %sets concentration of free proteins
    xB(1) = N-(n-1);    %initial values for a free lattice
    xAB(1) = 0;

    ProteinCount(1) = xAB(1);   %how many proteins are currently bound; 0 for initially empty lattice
    FracCover(1) = (xAB(1)*n)/N; %initial fractional coverage; 0 for an initially empty lattice

    Events = 0;
    Equilibrium = 0;
    
    while ~Equilibrium
        
        a_f(Events+1) = k_on*(L_A)*(xB(Events+1));   %propensity functions (probability of each event happening)
        a_r(Events+1) = k_off*(xAB(Events+1));
        a_0(Events+1) = a_f(Events+1)+a_r(Events+1);     %sum of propensity functions used for determining dt

        dt(Events+1) = (1/a_0(Events+1))*log(1/rand); %random time interval for Gillespie method

        if a_f(Events+1) > rand*a_0(Events+1) %a forward reaction occurs
            eventB = 0;
            while ~eventB   %repeats until a binding occurs
                SpotB = randi(N-(n-1));  %random location on lattice is chosen
                if DNA(SpotB:SpotB+(n-1)) == 0   %checks if location is free
                    DNA(SpotB:SpotB+(n-1)) = 1;   %binds protein to location
                    BoundAtSpot(SpotB) = 1;  %stores locatin in BoundAtSpot
                    BindCounter = BindCounter+1;    %updates bind counter
                    Hist(1,Events+1) = SpotB;      %stores location of binding event

                    FreeSpots = sum(max(find(diff([1 DNA 1])==1)-find(diff([1 DNA 1])==-1)-n+1,0)); %function to find number of free spots now available
                    xB((Events+1)+1) = FreeSpots;    %updates populations
                    xAB((Events+1)+1) = xAB(Events+1)+1;

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
            Hist(2,Events+1) = SpotU;  %stores location of unbinding event

            FreeSpots = sum(max(find(diff([1 DNA 1])==1)-find(diff([1 DNA 1])==-1)-n+1,0)); %calculates number of free spaces
            xB((Events+1)+1) = FreeSpots;
            xAB((Events+1)+1) = xAB(Events+1)-1;
        end
        ProteinCount((Events+1)+1) = sum(DNA)/n;   %all of these values should be integers
        FracCover((Events+1)+1) = sum(DNA)/N;
        v_loop(Events+1) = (sum(DNA)/n)/N; %binding density after each loop
        
        if Events > minIterations
            EqTest = [FracCover(Events-500) FracCover(Events)]; %state of system 100 loops ago and now
            FracCoverChange = diff(EqTest);    %difference between the two states
            if abs((FracCoverChange)) <= 0.0001  %in equilibrium if fractional coverage hasnt changed by more than 0.1% in last 100 events
                Equilibrium = 1;
            else
                Equilibrium = 0;
            end
        end
        
        t((Events+1)+1) = t(Events+1)+dt(Events+1);    %advances time
        Events = Events+1;
    end
    
%     figure(1);
%     scatter(t,FracCover,1,'o','filled');    %plots fractional coverage over time for all concentrations
%     hold on;
%     xlabel('Time');
%     ylabel('Fractional Coverage');
%     ylim([0 1]);
%     title('Fractional Coverage vs. Time');
%     legend(['L = ' num2str(m)]);
    
    L(Loops) = L_A;     %equilibrium free protein concentration (moles of free protein/volume = M)
    v(Loops) = v_loop(Events);      %equilibrium binding density
    ScatchY(Loops) = v(Loops)/L(Loops); %y-value for Scatchard plot
    
    BindingFrac(Loops) = BindCounter/Events;   %fraction of events which were binding for each loop
    UnbindingFrac(Loops) = UnbindCounter/Events;   %fraction of events which were unbinding for each loop
    AvgTime(Loops) = mean(dt);     %average time between events in each loop
    EqCover(Loops) = mean(FracCover(Events-500:Events));    %equilibrium fractional coverage for each loop
    TotalTime(Loops) = max(t);  %total time for each loop
    
    Loops = Loops+1;
end

x = 0:(1/n)/1000:1/n;    %x-values for theoretical Scatchard plot
TheorScatchY = K.*(1-(n.*x)).*(((1-(n.*x))./(1-((n-1).*x))).^(n-1));   %Eq. 10 in McGhee paper

figure();
subplot(2,1,1);
scatter(v,ScatchY,5,'r','filled');  %Scatchard plot benchmark (all runs)
hold on;
plot(x,TheorScatchY,'black');
xlabel('v');
xlim([0 1/n]);
ylabel('v/L');
ylim([0 max(ScatchY)+0.1]);
title(['Scatchard Plot (K=' num2str(K) ', N=' num2str(N) ', n=' num2str(n) ')']);
legend('Model Data','MVH Theory');

subplot(2,1,2);
scatter(A_Conc,EqCover,5,'g','filled');   %plots equilibrium coverage w/ respect to p
xlabel('Moles of Free Protein');
xlim([min(ProteinConc) max(ProteinConc)]);
ylabel('Equilibrium Coverage');
ylim([0 1]);
title('Equilibrium Coverage');

% fit data (ScatchY vs v) with this equation: K*(1-(n*v))*(((1-(n*v))/(1-((n-1)*v)))^(n-1))