clear all;
close all;

% This is the model that will take a single DNA lattice and bind/unbind
% proteins over time. The proteins in this model are noncooperative and
% dynamic Monte Carlo methods are included to have random time intervals
% between events. Cooperativity is included. Results are compared to the
% benchmark Scatchard plot and equation provided by McGhee and von Hippel
% (1974).

N = 1000;    %length of DNA lattice
n = 5;  %length of each protein
k_on = 0.1;   %kinetic rate constant for binding
k_off = 1;  %kinetic rate constant for unbinding
w = 0.1;      %cooperativity parameter

ProteinConc = [0.1:0.1:10];    %range of free protein concentrations used to form Scatchard plot
minIterations = 1.5*N;

K = k_on/k_off; %equilibrium constant

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

    xAB = zeros(1,minIterations+1);
    t = zeros(1,minIterations+1);
    a_If = zeros(1,minIterations);
    a_SCf = zeros(1,minIterations);
    a_DCf = zeros(1,minIterations);
    a_r = zeros(1,minIterations);
    a_0 = zeros(1,minIterations);
    dt = zeros(1,minIterations);
    FracCover = zeros(1,minIterations+1);
    v_loop = zeros(1,minIterations);
    BindHist = zeros(4,minIterations); %1: isolated binding, 2: SC binding, 3: DC binding, 4: unbinding
    xB_I = zeros(1,minIterations);
    xB_SC = zeros(1,minIterations);
    xB_DC = zeros(1,minIterations);
    
    t(1) = 0;

    xB_I(1) = N-(n-1);    %initial values for a free lattice
    xB_SC(1) = 0;
    xB_DC(1) = 0;
    xAB(1) = 0;
    
    FracCover(1) = sum(DNA)/N;
    
    BindCounter = 0;
    BindCounter_I = 0;
    BindCounter_SC = 0;
    BindCounter_DC = 0;
    UnbindCounter = 0;

    Equilibrium = 0;
    Events = 0;
    
    while Equilibrium == 0
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
        xB_I(Events+1) = Counter_I;    %amounts of each location, used in propensity functions
        xB_SC(Events+1) = Counter_SC;
        xB_DC(Events+1) = Counter_DC;
        
        a_If(Events+1) = k_on*(L_A(Loops))*(xB_I(Events+1));   %propensity functions (probability of each event happening)
        a_SCf(Events+1) = k_on*(L_A(Loops))*(xB_SC(Events+1))*w;
        a_DCf(Events+1) = k_on*(L_A(Loops))*(xB_DC(Events+1))*(w^2);
        a_r(Events+1) = k_off*(xAB(Events+1));
        a_0(Events+1) = a_If(Events+1)+a_SCf(Events+1)+a_DCf(Events+1)+a_r(Events+1);     %sum of propensity functions used for determining dt

        dt(Events+1) = (1/a_0(Events+1))*log(1/rand); %random time interval for Gillespie method
        R_1 = rand;
        
        if a_If(Events+1) > R_1*a_0(Events+1) %tests for isolated binding event
            pos_I = randi(length(Isolated));
            SpotB_I = Isolated(pos_I); %chooses a random location for binding to occur
            DNA(SpotB_I:SpotB_I+(n-1)) = 1; %binds protein
            BindCounter = BindCounter+1;
            BindCounter_I = BindCounter_I+1;
            BindHist(1,Events+1) = SpotB_I;
            BoundAtSpot(SpotB_I) = 1; %shows there's currently a protein bound at said location

            xAB((Events+1)+1) = xAB(Events+1)+1;    %updates populations
        elseif (a_If(Events+1)+a_SCf(Events+1)) > (R_1*a_0(Events+1))    %tests for singly contiguous binding event
            pos_SC = randi(length(SinglyContiguous));
            SpotB_SC = SinglyContiguous(pos_SC);  %chooses a random location for binding
            DNA(SpotB_SC:SpotB_SC+(n-1)) = 1; %binds protein
            BindCounter = BindCounter+1;
            BindCounter_SC = BindCounter_SC+1;
            BindHist(2,Events+1) = SpotB_SC;
            BoundAtSpot(SpotB_SC) = 1;    %shows there's currently a protein bound at location

            xAB((Events+1)+1) = xAB(Events+1)+1;    %updates populations
        elseif (a_If(Events+1)+a_SCf(Events+1)+a_DCf(Events+1)) > (R_1*a_0(Events+1))    %tests for doubly contiguous binding event
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
            
            xAB((Events+1)+1) = xAB(Events+1)-1;    %updates populations
        end
        
        v_loop(Events+1) = (sum(DNA)/n)/N; %calculates binding density after each loop
        FracCover((Events+1)+1) = sum(DNA)/N;  %fractional coverage
        
        if Events > minIterations
            EqTest = [FracCover(Events-500) FracCover(Events)]; %state of system 100 loops ago and now
            FracCoverChange = diff(EqTest);    %difference between the two states
            if abs((FracCoverChange)) <= 0.0001  %in equilibrium if fractional coverage hasnt changed by more than 0.1% in last 100 events
                Equilibrium = 1;
            else
                Equilibrium = 0;
            end
        end
        
        t((Events+1)+1) = t(Events+1)+dt(Events+1);
        Events = Events+1;
    end
    EventCount(Loops) = Events;
    
    L(Loops) = L_A(Loops);     %equilibrium free protein concentration (moles of free protein/volume = M)
    v(Loops) = v_loop(Events);      %equilibrium binding density
    ScatchY(Loops) = v(Loops)/L(Loops); %y-value for Scatchard plot
    
    BindingFrac(Loops) = BindCounter/Events;   %fraction of events which were binding for each loop
    UnbindingFrac(Loops) = UnbindCounter/Events;   %fraction of events which were unbinding for each loop
    AvgTime(Loops) = mean(dt);     %average time between events in each loop
    EqCover(Loops) = mean(FracCover(Events-500:Events));    %equilibrium fractional coverage for each loop
    TotalTime(Loops) = max(t);  %total time for each loop
    
%     if Loops <= 10
%         figure(1);      %figure to show if the first few iterations are reaching equilibrium
%         subplot(2,5,Loops);
%         scatter(t,FracCover,2,'g','filled');
%         hold on;
%         xlabel('Time, t');
%         ylabel('Fractional Coverage');
%         ylim([0 1]);
%         title('Fractional Coverage vs. Time');
%         legend([num2str(EventCount(Loops)) 'Events']);
%     end
    
    Loops = Loops+1;
end

X = 0:0.001:1/n;    %Theoretical Scatchard plot from McGhee paper
if w == 1
    TheorScatchY = K.*(1-(n.*X)).*(((1-(n.*X))./(1-((n-1).*X))).^(n-1));
else
    TheorR = sqrt(((1-((n+1).*X)).^2)+(4.*w.*X.*(1-(n.*X))));
    TheorScatchY = K.*(1-(n.*X)).*((((((2.*w)-1).*(1-(n.*X)))+X-TheorR)./(2.*(w-1).*(1-(n.*X)))).^(n-1)).*(((1-((n+1).*X)+TheorR)./(2.*(1-(n.*X)))).^2);
end

HL_X = 0:max(ProteinConc)/1000:max(ProteinConc);
Hill_Lang = ((HL_X).^(w))./((1/K)+((HL_X).^(w)));

figure();
% subplot(2,1,1);
scatter(v,ScatchY,5,'r','filled');  %Scatchard plot benchmark (all runs)
hold on;
plot(X,TheorScatchY,'b');
xlabel('v');
xlim([0 1/n]);
ylabel('v/L');
% ylim([0 max(ScatchY)+0.1]);
title(['Scatchard Plot']);
legend('Simulation','McGhee & von Hippel Model');

% subplot(2,1,2);       %Hill-Langmuir plot
% scatter(ProteinConc,EqCover,5,'g','filled');    %Fractional coverage vs. protein concentration
% hold on;
% plot(HL_X,Hill_Lang,'black');
% xlabel('Protein Concentration, L (M)');
% xlim([0 max(ProteinConc)]);
% ylabel('Fractional Coverage');
% ylim([0 1]);
% title('Hill-Langmuir Plot');
% legend('Model Data','Theoretical Values');

% Fit Scatchard data to following equation: K*(1-(n*v))*((((((2*w)-1)*(1-(n*v)))+v-(sqrt(((1-((n+1)*v))^2)+(4*w*v*(1-(n*v))))))/(2*(w-1)*(1-(n*v))))^(n-1))*(((1-((n+1)*v)+(sqrt(((1-((n+1)*v))^2)+(4*w*v*(1-(n*v))))))/(2*(1-(n*v))))^2)
%        search for values for K,n,v,w
% Fit Hill-Langmuir to following equation: ((H_X)^(w))/((1/K)+((H_X)^(w)))
%        search for values for K,w (w = Hill Coefficient)