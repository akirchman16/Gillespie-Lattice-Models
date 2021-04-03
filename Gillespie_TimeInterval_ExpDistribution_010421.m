clear all;
close all;

% Set to create exponential distribution of time intervals that are
% possibly selected from the Gillespie Algorithm. This is all for the first
% event so the lattice would initially be completely empty.

Iterations = 5000;

N = 8660;       %input parameters
n = 3;
k_on = 1;
k_off = 1;
L = 1;

a_1 = k_on*L*(N-(n-1)); %propensity function for binding
a_2 = k_off*0;          %propensity function for unbinding
a_0 = a_1+a_2;          %sum of propensity functions

tau = zeros(1,Iterations);

for i = 1:Iterations
    r_1(i) = rand;  %random number from uniform distribution of unit interval
    
    tau(i) = (1/a_0)*log(1/r_1(i));  %time interval
end

disp(['1/a_0 = ', num2str(1/a_0)]);
disp(['Distribution Mean = ', num2str(mean(tau))]);

figure(1);
% histogram of time intervals for the first reaction
histfit(tau,100,'exponential');
hold on;
xlabel('Time Interval, \tau');
xlim([0 inf]);
ylabel('Occurences');
title('Time Intervals');