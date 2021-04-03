clear all;
close all;

% This code will bind and unbind proteins to a DNA lattice over time. Each
% loop of the code will do one binding or unbinding event. Rather than
% previous codes that did multiple events over each loop, this code will
% only look at one protein at a time. An attempt to bind or unbind will
% occur at constant intervals of time. This will then be adjusted to have
% Gillespie timings. The reaction is modeled as an A+B<-->AB reaction and
% the proteins are non-interacting in this model.

N = 1000;    %length of DNA lattice
n = 3;  %length of each protein
k_on = 0.1;   %kinetic rate constant for binding
k_off = 0.1;  %kinetic rate constant for unbinding
Vol = 1;    %volume of system in which reaction occurs

Iterations = 1000;  %number of events which will occur
maxPop = 1000;  %maximum initial population of free proteins
dt = 1/maxPop;  %time between each binding event

K = k_on/k_off; %equilibrium constant
c_f = k_on/Vol; %constants for propensity functions
c_r = k_off;

for a = 1:maxPop
    DNA = zeros(1,N);   %empty DNA lattice
    CurrentBound = zeros(1,N);  %empty CurrentBound array
    
    xA(1) = a;  %initial populations for an empty lattice
    xB(1) = N-(n-1);
    xAB(1) = 0;
    
    t(1) = 0;
    
    a_f = zeros(1,Iterations);  %memory allocation to save time
    a_r = zeros(1,Iterations);
    
    for j = 1:Iterations
        a_f(j) = c_f*xA(j)*xB(j);   %propensity functions for each reaction
        a_r(j) = c_r*xAB(j);
        
        Bind(j) = randi(N-(n-1));   %picks random location for an event to possibly occur
        if DNA(Bind(j):Bind(j)+(n-1)) == 0  %if location is free
            if rand <= a_f(j)*dt   %tests unbinding probability
                DNA(Bind(j):Bind(j)+(n-1)) = 1; %binds protein
                CurrentBound(Bind(j)) = 1;  %stores location in CurrentBound
                
                xA(j+1) = xA(j)-1;  %updates populations
                FreeSpots = 0;  %used to calculate how many free locations there are on the lattice
                for h = 1:N-(n-1)
                    if DNA(h:h+(n-1)) == 0
                        FreeSpots = FreeSpots+1;
                    end
                end
                xB(j+1) = FreeSpots;
                xAB(j+1) = xAB(j)+1;
                vLoop(j+1) = (sum(DNA)/n)/N;  %updates binding density
            else
                xA(j+1) = xA(j);    %if no reaction occured the populations don't change
                xB(j+1) = xB(j);
                xAB(j+1) = xAB(j);
                vLoop(j+1) = (sum(DNA)/n)/N;  %updates binding density
            end
        elseif CurrentBound(Bind(j)) == 1   %tests if protein is bound at location
            if rand <= a_r(j)*dt   %tests unbinding probability
                DNA(Bind(j):Bind(j)+(n-1)) = 0; %unbinds protein
                CurrentBound(Bind(j)) = 0;  %removes location from current bound array
                
                xA(j+1) = xA(j)+1;  %updated populations
                FreeSpots = 0;
                for h = 1:N-(n-1)   %used to calculate free spots (xB) on lattice
                    if DNA(h:h+(n-1)) == 0
                        FreeSpots = FreeSpots+1;
                    end
                end
                xB(j+1) = FreeSpots;
                xAB(j+1) = xAB(j)-1;
                vLoop(j+1) = (sum(DNA)/n)/N;  %updates binding density
            else
                xA(j+1) = xA(j);    %no reaction occurs so no population change
                xB(j+1) = xB(j);
                xAB(j+1) = xAB(j);
                vLoop(j+1) = (sum(DNA)/n)/N;  %updates binding density
            end
        else
            xA(j+1) = xA(j);    %no reaction occurs so no population changes
            xB(j+1) = xB(j);
            xAB(j+1) = xAB(j);
            vLoop(j+1) = (sum(DNA)/n)/N;  %updates binding density
        end
        t(j+1) = t(j)+dt;   %advances time in constant steps
    end
    
    EqxA(a) = mean(xA(round(0.75*Iterations):Iterations));    %equilibrium populations
    EqxB(a) = mean(xB(round(0.75*Iterations):Iterations));
    EqxAB(a) = mean(xAB(round(0.75*Iterations):Iterations));

    v(a) = mean(vLoop(0.75*Iterations:Iterations));  %records equilibrium binding density
    L(a) = EqxA(a)/Vol; %calculates free protein concentration at equilibrium
    ScatchY(a) = v(a)/L(a); %calculates model Scatchard plot y-values
end

TheorScatchX = 0:(1/(1000*n)):(1/n);
TheorScatchY = K*(1-(n.*TheorScatchX)).*(((1-(n.*TheorScatchX))./(1-((n-1).*TheorScatchX))).^(n-1));

figure();
% subplot(3,1,1);
scatter(v,ScatchY,5,'r','filled');  %creates Scatchard plot to compare results
hold on;
plot(TheorScatchX,TheorScatchY,'black');
xlabel('v');
% xlim([0 1/n]);
ylabel('v/L');
% ylim([0 K+(0.25*K)]);
title('Scatchard Plot');

% subplot(3,1,2);
% scatter(t,xA,5,'r','filled');   %plots populations over time
% hold on;
% scatter(t,xAB,5,'g','filled');
% xlabel('Time, t');
% ylabel('Population, x');
% legend('Free Proteins','Bound Proteins');
% title('Populations vs. Time');
% 
% subplot(3,1,3);
% plot(t,vLoop,'black');    %plots binding density over time for last run
% xlabel('Time, t');
% ylabel('Binding Density, v');
% title('Binding Density vs. Time');