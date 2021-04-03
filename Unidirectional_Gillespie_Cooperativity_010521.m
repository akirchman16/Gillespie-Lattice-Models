clear all;
close all;

% This model will compare growth profiles of ssDNA as RAD51 binds and
% unbinds to the lattice until equilibrium is reached. It will compare
% growth profiles with different cooperativity values. The Gillespie
% algorithm will be used to model stochasticity. Cooperative binding will
% be unidirectional (towards larger numbers) and will only be expanding.

N = 8660;    %length of DNA lattice
n = 3;  %length of each protein
k_on = 1;   %kinetic rate constant for binding (units: 1/(nt*s)??? )
k_off = 1;  %kinetic rate constant for unbinding (units: 1/s??? )
L = 1;    %concentration of free proteins (units: Molarity???)

minIterations = 500;  %minimum number of evnts that will occur
K = k_on/k_off; %equilibrium constant (units: )

Cooperativity = [1,1,1,10,10,10,100,100,100]; %cooperativity values for the simulation
Loops = 0;

for w = Cooperativity
    Loops = Loops+1;

    DNA = zeros(1,N+2);   %empty DNA lattice with dummy zeros
    BoundAtSpot = zeros(1,N+2);  %empty CurrentBound array with dummy zeros

    xAB = zeros(1,minIterations+1); %memory allocation
    t = zeros(1,minIterations+1);
    a_If = zeros(1,minIterations);
    a_SCf = zeros(1,minIterations);
    a_DCf = zeros(1,minIterations);
    a_r = zeros(1,minIterations);
    a_0 = zeros(1,minIterations);
    dt = zeros(1,minIterations);
    FracCover = zeros(1,minIterations+1);
    BindHist = zeros(4,minIterations); %1: isolated binding, 2: SC binding, 3: DC binding, 4: unbinding
    xB_I = zeros(1,minIterations);
    xB_SC = zeros(1,minIterations);
    xB_DC = zeros(1,minIterations);
    Length_um = zeros(1,minIterations);
    Length_nm = zeros(1,minIterations);

    t(1) = 0;   %all initial values for variables
    xB_I(1) = N-(n-1);    %initial values for a free lattice
    xB_SC(1) = 0;
    xB_DC(1) = 0;
    xAB(1) = 0;
    FracCover(1) = sum(DNA)/N;
    x_Occupied(1) = 0;
    x_Unoccupied(1) = N;
    Length_nm(1) = 0.34*N;
    FracCoverStates = 0;
    FracCoverChange = 0;
    FilamentEnds_L = 0;
    FilamentEnds_R = 0;
    FilamentCount = 0;
    FilamentLengths = 0;
    Avg_FilamentLength = 0;
    BindCounter = 0;
    BindCounter_I = 0;
    BindCounter_SC = 0;
    BindCounter_DC = 0;
    UnbindCounter = 0;
    Equilibrium = 0;
    EqualTimes = 0;
    Events = 0;

    while (Equilibrium == 0) || (EqualTimes == 0)
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
                elseif (DNA(x-1) == 1 && DNA(x+n) == 0) %records all singly contiguous locations (expand only towards larger numbers)
                    SinglyContiguous(Counter_SC+1) = x;
                    Counter_SC = Counter_SC+1;
                elseif DNA(x-1) == 1 && DNA(x+n) == 1   %records all doubly contiguous locations
                    DoublyContiguous(Counter_DC+1) = x;
                    Counter_DC = Counter_DC+1;
                end
            end
        end
        xB_I(Events+1) = Counter_I;    %amounts of each location, used in propensity functions
        xB_SC(Events+1) = Counter_SC;
        xB_DC(Events+1) = Counter_DC;

        a_If(Events+1) = k_on*(L)*(xB_I(Events+1));   %propensity functions (probability of each event happening)
        a_SCf(Events+1) = k_on*(L)*(xB_SC(Events+1))*w;
        a_DCf(Events+1) = k_on*(L)*(xB_DC(Events+1))*(w^2);
        a_r(Events+1) = k_off*(xAB(Events+1));
        a_0(Events+1) = a_If(Events+1)+a_SCf(Events+1)+a_DCf(Events+1)+a_r(Events+1);     %sum of propensity functions used for determining dt

        dt(Events+1) = (1/a_0(Events+1))*log(1/rand); %random time interval for Gillespie method
        R_1(Events+1) = rand;

        if a_If(Events+1) > R_1(Events+1)*a_0(Events+1) %tests for isolated binding event
            pos_I = randi(length(Isolated));
            SpotB_I = Isolated(pos_I); %chooses a random location for binding to occur
            DNA(SpotB_I:SpotB_I+(n-1)) = 1; %binds protein
            BindCounter = BindCounter+1;
            BindCounter_I = BindCounter_I+1;
            BindHist(1,Events+1) = SpotB_I;
            BoundAtSpot(SpotB_I) = 1; %shows there's currently a protein bound at said location

            xAB((Events+1)+1) = xAB(Events+1)+1;    %updates populations
        elseif (a_If(Events+1)+a_SCf(Events+1)) > (R_1(Events+1)*a_0(Events+1))    %tests for singly contiguous binding event
            pos_SC = randi(length(SinglyContiguous));
            SpotB_SC = SinglyContiguous(pos_SC);  %chooses a random location for binding
            DNA(SpotB_SC:SpotB_SC+(n-1)) = 1; %binds protein
            BindCounter = BindCounter+1;
            BindCounter_SC = BindCounter_SC+1;
            BindHist(2,Events+1) = SpotB_SC;
            BoundAtSpot(SpotB_SC) = 1;    %shows there's currently a protein bound at location

            xAB((Events+1)+1) = xAB(Events+1)+1;    %updates populations
        elseif (a_If(Events+1)+a_SCf(Events+1)+a_DCf(Events+1)) > (R_1(Events+1)*a_0(Events+1))    %tests for doubly contiguous binding event
            pos_DC = randi(length(DoublyContiguous));
            SpotB_DC = DoublyContiguous(pos_DC);  %chooses a random location for binding
            DNA(SpotB_DC:SpotB_DC+(n-1)) = 1; %binds protein
            BindCounter = BindCounter+1;
            BindCounter_DC = BindCounter_DC+1;
            BindHist(3,Events+1) = SpotB_DC;
            BoundAtSpot(SpotB_DC) = 1;    %shows there's currently a protein bound at location

            xAB((Events+1)+1) = xAB(Events+1)+1;    %updates population
        else    %otherwise an unbinding event occurs
            CurrentBound = find(BoundAtSpot == 1);
            pos = randi(length(CurrentBound));
            SpotU = CurrentBound(pos);  %selects an already bound protein to unbind

            DNA(SpotU:SpotU+(n-1)) = 0; %unbinds protein
            BoundAtSpot(SpotU) = 0;
            UnbindCounter = UnbindCounter+1;
            BindHist(4,Events+1) = SpotU;

            xAB((Events+1)+1) = xAB(Events+1)-1;    %updates populations
        end

        FracCover((Events+1)+1) = sum(DNA)/N;  %fractional coverage

        if Events > minIterations && Equilibrium == 0
            FracCoverStates = (FracCover((Events-floor(0.25*Events)):1:Events));   %every 10 states for the last 20% of the model
            FracCoverChange = abs(diff(FracCoverStates));   %difference among each of those states
            if (mean(FracCoverChange) <= n/N) && (abs(FracCover(Events-floor(0.25*Events))-FracCover(Events)) <= 3*n/N)
                Equilibrium = 1;
            else
                Equilibrium = 0;
            end
        end

        t((Events+1)+1) = t(Events+1)+dt(Events+1);
        Events = Events+1;
        
        if Loops > 1       %makes sure each cooperativity value runs apporximately the same amount of time
            if max(t) >= max(TotalTime)
                EqualTimes = 1;
            else
                EqualTimes = 0;
            end
        elseif Loops == 1
            EqualTimes = 1;
        end

        x_Occupied(Events+1) = length(find(DNA ~= 0));
        x_Unoccupied(Events+1) = length(find(DNA(2:N+1) == 0));
        Length_nm(Events+1) = (0.34*x_Unoccupied(Events+1))+(0.51*x_Occupied(Events+1));   %length of molecule in nm

        Gillespie_Probabilities(1,Events) = a_If(Events)*dt(Events);      %probabilities of isolated binding
        Gillespie_Probabilities(2,Events) = a_SCf(Events)*dt(Events);     %probabilities of singly contiguous binding
        Gillespie_Probabilities(3,Events) = a_DCf(Events)*dt(Events);     %probabilities of doubly congicuous binding
        Gillespie_Probabilities(5,Events) = sum(Gillespie_Probabilities(1:3,Events));   %probabilities of binding in general
        Gillespie_Probabilities(7,Events) = a_r(Events)*dt(Events);       %probabilities of unbinding
        
        True_Probabilities(1,Events) = a_If(Events)/a_0(Events);
        True_Probabilities(2,Events) = a_SCf(Events)/a_0(Events);
        True_Probabilities(3,Events) = a_DCf(Events)/a_0(Events);
        True_Probabilities(5,Events) = (a_If(Events)+a_SCf(Events)+a_DCf(Events))/a_0(Events);
        True_Probabilities(7,Events) = a_r(Events)/a_0(Events);
        
        FilamentEnds_L = 1+find(diff(DNA)==1);  %locations of the left end of each filament
        FilamentEnds_R = 1+find(diff(DNA)==-1); %locations of the right end of each filament
        FilamentCount = length(FilamentEnds_R); %number of filaments on the lattice
        FilamentLengths(Events,1:FilamentCount) = FilamentEnds_R-FilamentEnds_L;    %length of filaments
        Avg_FilamentLength(Events) = mean(FilamentLengths(Events,1:FilamentCount)); %average filament length for each event
    end

    Length_um = Length_nm/1e+3; %length of molecule in um
    Equilibrium_Coverage(Loops) = mean(FracCoverStates);
    disp(['(', num2str(w),') - Equilibrium Saturation = ', num2str(Equilibrium_Coverage(Loops))]);

    TotalTime(Loops) = max(t);
    LoopLengths = sort(TotalTime,'descend');

    figure(1);
    
%     Plot of fractional coverage of the DNA
    
%     subplot(2,1,1);
    scatter(t,FracCover,2,'filled');
    hold on;
    xlabel('Time, t');
    xlim([0 max(LoopLengths)]);
    ylabel('Fractional Coverage');
    ylim([0 1]);
%     yline(Equilibrium_Coverage(Loops),'k',['\omega = ', num2str(w)]);
    title('Unidirectional Saturation of DNA Lattice');

%     plot of Length of Molecule over time (using values of 0.34nm and
%     0.51nm from van der Heijden paper)

%     subplot(2,1,2);
%     scatter(t,Length_um,2,'filled');
%     hold on;
%     xlim([0 max(TotalTime)]);
%     xlabel('Time, t (s)');
%     ylabel('Length (\mum)');
%     title('Length of DNA Molecule');
    
%     plot of average filament length over time
    
%     subplot(2,1,2);
%     scatter(t(2:length(t)),Avg_FilamentLength,2,'filled');
%     hold on;
%     xlabel('Time, t');
%     xlim([0 max(TotalTime)]);
%     ylabel('Average Filament Length (Sites)');
%     ylim([0 inf]);
%     title('Average Filament Length');
end

figure(1);
% subplot(2,1,1);
Legend = cell(length(Cooperativity),1);
for c = 1:length(Cooperativity)
    Legend{c} = ['\omega = ', num2str(Cooperativity(c))];
end
legend(Legend,'location','southeast');