clear all;
close all;

% This code will model noncooperative binding over time. Unbinding will be
% ignored in order to match Fig. 2B of the van der Heijden paper. This plot
% will be fractional coverage vs. time. No Gillespie timing will be
% involved.

N = 100;    %length of DNA lattice
n = 3;  %length of each protein
k_on = 0.1;   %kinetic rate constant for binding
k_off = 0.1;  %kinetic rate constant for unbinding
Vol = 1;    %volume of system in which reaction occurs

Iterations = 500;  %number of events which will occur
Pop = 100;  %maximum initial population of free proteins

K = k_on/k_off; %equilibrium constant
c_f = k_on/Vol; %constants for propensity functions
c_r = k_off;

% for n = 1:25
    DNA = zeros(1,N);   %empty DNA lattice
    CurrentBound = zeros(1,N);  %empty CurrentBound array

    xA(1) = Pop;  %initial populations for an empty lattice
    xB(1) = N-(n-1);
    xAB(1) = 0;

    dt = 0.01;

    t(1) = 0;
    BindCounter = 0;
    UnbindCounter = 0;

    a_f = zeros(1,Iterations);  %memory allocation to save time
    a_r = zeros(1,Iterations);
    a_0 = zeros(1,Iterations);

    for j = 1:Iterations
        a_f(j) = c_f*xA(j)*xB(j);   %propensity functions for each reaction

        Bind(j) = randi(N-(n-1));   %picks random location for an event to possibly occur
        if DNA(Bind(j):Bind(j)+(n-1)) == 0  %if location is free
            if a_f(j)*dt >= rand   %tests binding probability (see explanation in van der Heijden paper)
                DNA(Bind(j):Bind(j)+(n-1)) = 1; %binds protein
                CurrentBound(Bind(j)) = 1;  %stores location in CurrentBound
                BindCounter = BindCounter+1;    %counts binding events

                xA(j+1) = xA(j)-1;  %updates populations
                FreeSpots = 0;  %used to calculate how many free locations there are on the lattice
                for h = 1:N-(n-1)
                    if DNA(h:h+(n-1)) == 0
                        FreeSpots = FreeSpots+1;
                    end
                end
                xB(j+1) = FreeSpots;
                xAB(j+1) = xAB(j)+1;

                FracCover(j+1) = (xAB(j+1)*n)/N; %calculates fractional coverage after each loop
            else
                xA(j+1) = xA(j);    %if no reaction occured the populations don't change
                xB(j+1) = xB(j);
                xAB(j+1) = xAB(j);

                FracCover(j+1) = (xAB(j+1)*n)/N; %calculates fractional coverage after each loop
            end
        else
            xA(j+1) = xA(j);    %if no reaction occured the populations don't change
            xB(j+1) = xB(j);
            xAB(j+1) = xAB(j);

            FracCover(j+1) = (xAB(j+1)*n)/N; %calculates fractional coverage after each loop
        end
        t(j+1) = t(j)+dt;   %advances time in constant steps
    end
%     AppSize(n) = N/(sum(DNA)/n);    %apparent binding size
%     BindSize(n) = n;
% end

% x = 1:25/1000:25; %x-value for theoretical plot
% y = 1.295*x;    %y-value for theoretical plot

figure();
subplot(2,1,1);
scatter(t,xA,5,'r','filled');   %plots protein populations
hold on;
scatter(t,xAB,5,'g','filled');
xlabel('Time');
ylabel('Populations');
title('Protein Populations vs. Time');
legend('Free Proteins','Bound Proteins');

subplot(2,1,2);
scatter(t,FracCover,5,'r','filled');    %plots fractional coverage over time
hold on;
xlabel('Time');
xlim([0 max(t)]);
ylabel('Fractional Coverage');
ylim([0 1]);
title('Fractional Coverage vs. Time');

% figure();
% scatter(BindSize,AppSize,5,'r','filled');
% hold on;
% plot(x,y,'black');
% xlabel('Actual Binding Size');
% xlim([0 max(n)]);
% ylabel('Apparent Binding Size');
% title('Apparent Binding Size');
% legend('Data','Paper Fit');