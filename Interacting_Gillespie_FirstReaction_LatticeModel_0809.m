clear all;
close all;

% This is the model that will take a single DNA lattice and bind/unbind
% proteins over time. The proteins in this model are noncooperative and
% dynamic Monte Carlo methods are included to have random time intervals
% between events. Cooperativity is included. Results are compared to the
% benchmark Scatchard plot and equation provided by McGhee and von Hippel
% (1974).

N = 1000;    %length of DNA lattice
n = 1;  %length of each protein
k_on = 1;   %kinetic rate constant for binding
k_off = 1;  %kinetic rate constant for unbinding
w = 2;      %cooperativity parameter

Iterations = 2000;  %number of events which will occur
ProteinConc = [0.1:0.1:10];    %range of free protein concentrations used to form Scatchard plot

K = k_on/k_off; %equilibrium constant
Locations = 2:1:N-(n-1)+1;  %list of all possible locations on the lattice

Loops = 1;

L = zeros(1,length(ProteinConc));   %memory allocation
v = zeros(1,length(ProteinConc));
ScatchY = zeros(1,length(ProteinConc));
BindingFrac = zeros(1,length(ProteinConc));
UnbindingFrac = zeros(1,length(ProteinConc));
AvgTime = zeros(1,length(ProteinConc));
TotalTime = zeros(1,length(ProteinConc));
EqCover = zeros(1,length(ProteinConc));

for m = ProteinConc
    DNA = zeros(1,N+2);   %empty DNA lattice with dummy zeros
    BoundAtSpot = zeros(1,N+2);  %empty CurrentBound array with dummy zeros
    
    L_A(Loops) = m;  %sets free protein concentration for this loop

    xAB = zeros(1,Iterations+1);
    t = zeros(1,Iterations+1);
    a_If = zeros(1,Iterations);
    a_SCf = zeros(1,Iterations);
    a_DCf = zeros(1,Iterations);
    a_r = zeros(1,Iterations);
    a_0 = zeros(1,Iterations);
    dt = zeros(1,Iterations);
    FracCover = zeros(1,Iterations+1);
    BindProb = zeros(Iterations,N); %probability of protein binding to each location if a binding event were to occur next
                                    %rows are DNA lattices after each iteration
    v_loop = zeros(1,Iterations);
    BindHist = zeros(4,Iterations); %1: isolated binding, 2: SC binding, 3: DC binding, 4: unbinding
    xB_I = zeros(1,Iterations);
    xB_SC = zeros(1,Iterations);
    xB_DC = zeros(1,Iterations);
    
    t(1) = 0;

    xB_I(1) = N-(n-1);    %initial values for a free lattice
    xB_SC(1) = 0;
    xB_DC(1) = 0;
    xAB(1) = 0;
    
    BindCounter = 0;
    BindCounter_I = 0;
    BindCounter_SC = 0;
    BindCounter_DC = 0;
    UnbindCounter = 0;

    for i = 1:Iterations
        Isolated = 0;   %preparing arrays for types of available locations
        SinglyContiguous = 0;
        DoublyContiguous = 0;
        Counter_I = 0;  %resets counters to count types of locations
        Counter_SC = 0;
        Counter_DC = 0;
        for x = 2:N-(n-1)+1
            if DNA(x:x+(n-1)) == 0
                if DNA(x-1) == 0 && DNA(x+n) == 0   %records all isolated locations
                    Isolated(Counter_I+1) = x;
                    Counter_I = Counter_I+1;
                elseif (DNA(x-1) == 0 && DNA(x+n) == 1) || (DNA(x-1) == 1 && DNA(x+n) == 0) %records all singly contiguous locations
                    SinglyContiguous(Counter_SC+1) = x;
                    Counter_SC = Counter_SC+1;
                elseif DNA(x-1) == 1 && DNA(x+n) == 1   %records all doubly contiguous locations
                    DoublyContiguous(Counter_DC+1) = x;
                    Counter_DC = Counter_DC+1;
                end
            end
        end
        xB_I(i) = Counter_I;    %amounts of each location, used in propensity functions
        xB_SC(i) = Counter_SC;
        xB_DC(i) = Counter_DC;
        xAB(i) = sum(DNA)/n;
        
        a_If(i) = k_on*(L_A(Loops))*(xB_I(i));   %propensity functions (probability of each event happening)
        a_SCf(i) = k_on*(L_A(Loops))*(xB_SC(i))*w;
        a_DCf(i) = k_on*(L_A(Loops))*(xB_DC(i))*(w^2);
        a_r(i) = k_off*(xAB(i));
        a_0(i) = a_If(i)+a_SCf(i)+a_DCf(i)+a_r(i);     %sum of propensity functions used for determining dt

        dt_I(i) = (1/a_If(i))*log(1/rand); %random time intervals for first reaction methond
        dt_SC(i) = (1/a_SCf(i))*log(1/rand);
        dt_DC(i) = (1/a_DCf(i))*log(1/rand);
        dt_U(i) = (1/a_r(i))*log(1/rand);
        
        dt(i) = min([dt_I(i),dt_SC(i),dt_DC(i),dt_U(i)]);    %useful time interval (smallest)
        
        R_1 = rand;
        
        if dt(i) == dt_I(i) %tests for isolated binding event
            pos_I = randi(length(Isolated));
            SpotB_I = Isolated(pos_I); %chooses a random location for binding to occur
            DNA(SpotB_I:SpotB_I+(n-1)) = 1; %binds protein
            BindCounter = BindCounter+1;
            BindCounter_I = BindCounter_I+1;
            BindHist(1,i) = SpotB_I;
            BoundAtSpot(SpotB_I) = 1; %shows there's currently a protein bound at said location

            xAB(i+1) = xAB(i)+1;    %updates populations
        elseif dt(i) == dt_SC(i)    %tests for singly contiguous binding event
            pos_SC = randi(length(SinglyContiguous));
            SpotB_SC = SinglyContiguous(pos_SC);  %chooses a random location for binding
            DNA(SpotB_SC:SpotB_SC+(n-1)) = 1; %binds protein
            BindCounter = BindCounter+1;
            BindCounter_SC = BindCounter_SC+1;
            BindHist(2,i) = SpotB_SC;
            BoundAtSpot(SpotB_SC) = 1;    %shows there's currently a protein bound at location
        elseif dt(i) == dt_DC(i)   %tests for doubly contiguous binding event
            pos_DC = randi(length(DoublyContiguous));
            SpotB_DC = DoublyContiguous(pos_DC);  %chooses a random location for binding
            DNA(SpotB_DC:SpotB_DC+(n-1)) = 1; %binds protein
            BindCounter = BindCounter+1;
            BindCounter_DC = BindCounter_DC+1;
            BindHist(3,i) = SpotB_DC;
            BoundAtSpot(SpotB_DC) = 1;    %shows there's currently a protein bound at location
        elseif dt(i) == dt_U(i)    %otherwise an unbinding event occurs
            CurrentBound = find(BoundAtSpot == 1);
            pos = randi(length(CurrentBound));
            SpotU = CurrentBound(pos);  %selects an already bound protein to unbind

            DNA(SpotU:SpotU+(n-1)) = 0; %unbinds protein
            BoundAtSpot(SpotU) = 0;
            UnbindCounter = UnbindCounter+1;
        end
        v_loop(i) = (sum(DNA)/n)/N; %calculates binding density after each loop
        FracCover(i) = sum(DNA)/N;  %fractional coverage
    end
    L(Loops) = L_A(Loops);     %equilibrium free protein concentration (moles of free protein/volume = M)
    v(Loops) = v_loop(Iterations);      %equilibrium binding density
    ScatchY(Loops) = v(Loops)/L(Loops); %y-value for Scatchard plot
    
    BindingFrac(Loops) = BindCounter/Iterations;   %fraction of events which were binding for each loop
    UnbindingFrac(Loops) = UnbindCounter/Iterations;   %fraction of events which were unbinding for each loop
    AvgTime(Loops) = mean(dt);     %average time between events in each loop
    EqCover(Loops) = mean(FracCover(round(0.75*Iterations):Iterations));    %equilibrium fractional coverage for each loop
    TotalTime(Loops) = max(t);  %total time for each loop
    
    Loops = Loops+1;
end

X = 0:0.001:1/n;    %Theoretical Scatchard plot from McGhee paper
TheorR = sqrt(((1-((n+1).*X)).^2)+(4.*w.*X.*(1-(n.*X))));
TheorScatchY = K.*(1-(n.*X)).*((((((2.*w)-1).*(1-(n.*X)))+X-TheorR)./(2.*(w-1).*(1-(n.*X)))).^(n-1)).*(((1-((n+1).*X)+TheorR)./(2.*(1-(n.*X)))).^2);
TheorScatchY_w1 = K.*(1-(n.*X)).*(((1-(n.*X))./(1-((n-1).*X))).^(n-1));

figure();
% subplot(2,1,1);
scatter(v,ScatchY,5,'r','filled');  %Scatchard plot benchmark (all runs)
hold on;
plot(X,TheorScatchY,'black');
plot(X,TheorScatchY_w1,'--black');
xlabel('v');
xlim([0 1/n]);
ylabel('v/L');
ylim([0 max(ScatchY)+0.1]);
title(['Scatchard Plot (K = ' num2str(K) ', N = ' num2str(N) ', n = ' num2str(n) ', w = ' num2str(w) ')']);
legend('Model Data','MVH Cooperative','MVH NonCooperative');

% subplot(2,1,2);
% scatter(L_A,EqCover,5,'g','filled');   %plots equilibrium coverage w/ respect to p
% xlabel('Moles of Free Protein');
% xlim([min(ProteinConc) max(ProteinConc)]);
% ylabel('Equilibrium Coverage');
% ylim([0 1]);
% title('Equilibrium Coverage');