clear all;
close all;

% This is the model that will take a single DNA lattice and bind/unbind
% proteins over time. The proteins in this model are noncooperative and
% dynamic Monte Carlo methods are included to have random time intervals
% between events. It will compare results with the dynamic model to the
% non-dynamic model.

N = 1000;    %length of DNA lattice
n = 3;  %length of each protein
k_on = 1;   %kinetic rate constant for binding
k_off = 1;  %kinetic rate constant for unbinding
L = 0.5;  %concentration of free proteins (Molarity)

minIterations = 1000;  %minimum number of events which will occur

K = k_on/k_off; %equilibrium constant

DNA = zeros(1,N);   %empty DNA lattice
BoundAtSpot = zeros(1,N);  %empty CurrentBound array

xB = zeros(1,minIterations+1);  %memory allocation
xAB = zeros(1,minIterations+1);
t = zeros(1,minIterations+1);
a_f = zeros(1,minIterations);
a_r = zeros(1,minIterations);
a_0 = zeros(1,minIterations);
dt = zeros(1,minIterations);
Hist = zeros(2,minIterations); %top row is binding event, bottom row is unbinding events
ProteinCount = zeros(1,minIterations+1);
FracCover = zeros(1,minIterations+1);

t(1) = 0;   %initial values
BindCounter = 0;
UnbindCounter = 0;
xB(1) = N-(n-1);    %initial values for a free lattice
xAB(1) = 0;
ProteinCount(1) = xAB(1);   %how many proteins are currently bound; 0 for initially empty lattice
FracCover(1) = (xAB(1)*n)/N; %initial fractional coverage; 0 for an initially empty lattice
Equilibrium = 0;
Events = 0;

while ~Equilibrium
    a_f(Events+1) = k_on*(L)*(xB(Events+1));   %propensity functions (probability of each event happening)
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
    FracCover((Events+1)+1) = sum(DNA)/N;    %fractional coverage of the DNA lattice
    
    t((Events+1)+1) = t(Events+1)+dt(Events+1);    %advances time
    
    if Events > minIterations
        EqTest = [FracCover(Events-500) FracCover(Events)]; %state of system 500 loops ago and now
        FracCoverChange = diff(EqTest);    %difference between the two states
        if abs((FracCoverChange)) <= 0.0001  %in equilibrium if fractional coverage hasnt changed by more than 0.01% in last 100 events
            Equilibrium = 1;
        else
            Equilibrium = 0;
        end
    end
    
    Events = Events+1;
end

figure();
scatter(t,FracCover,2,'r','filled');    %fractional coverage vs. dynamic time
xlabel('Time, t (s)');
xlim([0 max(t)]);
ylabel('Fractional Coverage');
ylim([0 1]);
title(['Fractional Coverage (K = ' num2str(K) ', N = ' num2str(N) ', n = ' num2str(n) ')']);
legend([ num2str(Events) ' Events']);