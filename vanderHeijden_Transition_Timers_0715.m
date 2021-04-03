clear all;
close all;

% This code should create some distribution based on method described in
% van der Heijden paper: transition occurs when transition probability is
% larger than a random value from unit uniform distribution.

dt = 0.01;  %time to each interaction
p_1 = 0.5;

for i = 1:1000
    event = 0;
    t = 0;
    while event == 0
        if p_1 >=rand
            event = 1;
            time(i) = t;
        else
            t = t+dt;
        end
    end
end

figure();
histogram(time,100);
hold on;
xlabel('Time to Transition');