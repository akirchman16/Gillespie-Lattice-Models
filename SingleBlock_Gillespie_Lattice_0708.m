clear all;
close all;

% This code will model proteins binding to and unbinding from a DNA lattice
% with timings based on Gillespie algorithm. The proteins are
% non-interacting. For this simple model we will only be looking at a
% single location on the lattice which is large enough for the protein to
% bind to. In this code A represents a free protein, B represents a free
% location on the DNA lattice that is large enough to hold a protein, and
% AB represents a bound protein. It then plots the populations of species
% over time. Since there is only one location for binding/unbinding the
% reaction should simply alternate between a binding and unbinding event (j
% alternates between 1 and 2) and the populations should only vary by one.

N = 3;  %length of DNA lattice
n = 3;  %length of protein
k_on = 1;   %kinetic rate constant for binding
k_off = 1;  %kinetic rate constant for unbinding
omega = 1;  %volume of system
Iterations = 100;   %how many times the code runs

t(1) = 0;   %initializes time at 0
DNA = zeros(1,N);   %lattice to model DNA

xA(1) = 5;          %initial amounts of each species
xB(1) = N-(n-1);
xAB(1) = 0;

c_f = k_on/omega;   %constants needed for propensity functions
c_r = k_off;

BindCounter = 0;    %counters for binding and unbinding events
UnbindCounter = 0;

for a = 1:Iterations
    a_f(a) = c_f*xA(a)*xB(a);   %propensity functions for each reaction
    a_r(a) = c_r*xAB(a);
    a_0(a) = a_f(a)+a_r(a); %sum of propensity functions
    
    r_1(a) = rand;  %random numbers for Monte Carlo method
    r_2(a) = rand;
    
    tau(a) = (1/a_0(a))*log(1/r_1(a)); %random time step through Monte Carlo
    if a_f(a) > r_2(a)*a_0(a)   %which reaction will occur
        j(a) = 1;  %binding event
    else
        j(a) = 2;  %unbinding event
    end
    if j(a) == 1   %if a binding occurs
        DNA(1:n) = 1;   %binds protein
        BindCounter = BindCounter+1;    %increases counter
        xA(a+1) = xA(a)-1;  %updates populations
        xB(a+1) = xB(a)-1;
        xAB(a+1) = xAB(a)+1;
    elseif j(a) == 2   %if an unbinding occurs
        DNA(1:n) = 0;   %unbinds protein
        UnbindCounter = UnbindCounter+1;    %increases counter
        xA(a+1) = xA(a)+1;  %updates populations
        xB(a+1) = xB(a)+1;  
        xAB(a+1) = xAB(a)-1;
    end
    t(a+1) = t(a)+tau(a);   %advances time
end

figure();   %plots
scatter(t,xA,5,'r','filled');
hold on;
scatter(t,xB,5,'g','filled');
scatter(t,xAB,5,'b','filled');
xlabel('Time,t');
ylabel('Populations');
legend('xA','xB','xAB');