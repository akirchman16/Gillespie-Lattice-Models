clear all;
close all;

% This is the model that will take a single DNA lattice and bind/unbind
% proteins over time. The proteins in this model are noncooperative and
% dynamic Monte Carlo methods are included to have random time intervals
% between events. Cooperativity is included by recording the probabilities
% of binding to each location on the lattice over time.

N = 100;    %length of DNA lattice
n = 3;  %length of each protein
k_on = 1;   %kinetic rate constant for binding
k_off = 1;  %kinetic rate constant for unbinding
w = 5;      %cooperativity parameter
L_A = 0.5;  %concentration of free proteins (Molarity)
Iterations = 100;  %number of events which will occur

K = k_on/k_off; %equilibrium constant
Locations = 2:1:N-(n-1)+1;  %list of all possible locations on the lattice

DNA = zeros(1,N+2);   %empty DNA lattice with dummy zeros
BoundAtSpot = zeros(1,N+2);  %empty CurrentBound array with dummy zeros

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
BindProb = zeros(Iterations,N); %probability of protein binding to each location if a binding event were to occur next
                                %rows are DNA lattices after each iteration
CoopEvents = zeros(3,Iterations);   %top row is isolated events, second row is singly contiguous, third row is doubly contiguous

t(1) = 0;
BindCounter = 0;
UnbindCounter = 0;

xB(1) = N-(n-1);    %initial values for a free lattice
xAB(1) = 0;

ProteinCount(1) = xAB(1);   %how many proteins are currently bound; 0 for initially empty lattice
FracCover(1) = (xAB(1)*n)/N; %initial fractional coverage; 0 for an initially empty lattice


for i = 1:Iterations
    a_f(i) = k_on*(L_A)*(xB(i));   %propensity functions (probability of each event happening)
    a_r(i) = k_off*(xAB(i));
    a_0(i) = a_f(i)+a_r(i);     %sum of propensity functions used for determining dt

    dt(i) = (1/a_0(i))*log(1/rand); %random time interval for Gillespie method
    
    for j = 2:N-(n-1)+1   %loop to calculate probability of binding at each spot if binding event occurs
        if (DNA(j) ~= 1) && (sum(DNA(j:j+(n-1))) == 0)
            if DNA(j-1) == 0 && DNA(j+n) == 0   %checks for isolated location
                BindProb(i,j-1) = 1/(xB(i));
            elseif (DNA(j-1) == 0 && DNA(j+n) == 1) || (DNA(j-1) == 1 && DNA(j+n) == 0) %checks for singly contiguous location
                BindProb(i,j-1) = w/(xB(i));
            elseif DNA(j-1) == 1 && DNA(j+n) == 1 %checks for doubly contiguous location
                BindProb(i,j-1) = (w^2)/(xB(i));
            end
        end
    end

    if a_f(i) > rand*a_0(i) %a forward reaction occurs
        eventB = 0;
        while ~eventB   %repeats until a binding occurs
            SpotB = randsample(Locations,1,true,BindProb(i,1:98));  %random location on lattice is chosen
            if DNA(SpotB:SpotB+(n-1)) == 0   %checks if location is free
               DNA(SpotB:SpotB+(n-1)) = 1;   %binds protein to location
               BoundAtSpot(SpotB) = 1;  %stores locatin in BoundAtSpot
               BindCounter = BindCounter+1;    %updates bind counter
               Hist(1,i) = SpotB;      %stores location of binding event

               FreeSpots = sum(max(find(diff([1 DNA 1])==1)-find(diff([1 DNA 1])==-1)-n+1,0)); %function to find number of free spots now available
               xB(i+1) = FreeSpots;    %updates populations
               xAB(i+1) = xAB(i)+1;

               eventB = 1;
            end
        end
        if DNA(SpotB-1) == 0 && DNA(SpotB+n) == 0 %isolated event occured
            CoopEvents(1,i) = SpotB;
        elseif (DNA(SpotB-1) == 0 && DNA(SpotB+n) == 1) || (DNA(SpotB-1) == 1 && DNA(SpotB+n) == 0) %singly contiguous event occured
            CoopEvents(2,i) = SpotB;
        elseif DNA(j-1) == 1 && DNA(j+n) == 1   %doubly contiguous binding event occured
            CoopEvents(3,i) = SpotB;
        end
    else                       %otherwise an unbinding has to occur
        CurrentBound = find(BoundAtSpot == 1);
        pos = randi(length(CurrentBound));
        SpotU = CurrentBound(pos);  %random position for protein to unbind

        DNA(SpotU:SpotU+(n-1)) = 0; %unbinds protein
        BoundAtSpot(SpotU) = 0; %removes location from BoundAtSpot
        UnbindCounter = UnbindCounter+1;    %updates unbind counter
        Hist(2,i) = SpotU;  %stores location of unbinding event

        FreeSpots = sum(max(find(diff([1 DNA 1])==1)-find(diff([1 DNA 1])==-1)-n+1,0)); %calculates number of free spaces
        xB(i+1) = FreeSpots;
        xAB(i+1) = xAB(i)-1;
    end
    t(i+1) = t(i)+dt(i);    %advances time

    ProteinCount(i+1) = sum(DNA)/n;   %all of these values should be integers
    FracCover(i+1) = sum(DNA)/N;    %fractional coverage of the DNA lattice
end

% figure();
% scatter(t,FracCover,2,'r','filled');    %fractional coverage vs. dynamic time
% xlabel('Time, t (s)');
% xlim([0 max(t)]);
% ylabel('Fractional Coverage');
% ylim([0 1]);
% title(['Fractional Coverage (K = ' num2str(K) ', N = ' num2str(N) ', n = ' num2str(n) ')']);
% legend([ num2str(Iterations) ' Events']);