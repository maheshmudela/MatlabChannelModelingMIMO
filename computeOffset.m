function  [theta_AOA_deg] = computeOffset(AOA_deg, AS)
% added for test.
%AOA_deg = randi([0,360],1,8);
%AS      = 2;

%for every AOA OR AOD, say N MS, compute phase of every M path for given n.
% here offset as cos(theat + phi), these phi, is based AS , angular spread of
% of any MS
theta = PasSubArry(AS);
%disp(theta);

M = length(theta);
N = length(AOA_deg);
for n=1:N
 for m=1:M

     % creating  two col vector for every m, with one value aoa + theta, other aoa -theta
     % theta_AOA =LOS angle +- phi;
     theta_AOA_deg(n,[2*m-1:2*m]) = AOA_deg(n) + [theta(m)  -theta(m)];


  end

end

%figure;
%t=1:1:N;
%plot(t,theta_AOA_deg(t,:));

figure;
plot(theta_AOA_deg);
title('Calculated AOA values across iterations (n)');
xlabel('Index n');
ylabel('Angle (degrees)');
legend('Col 1', 'Col 2', 'Col 3', 'Col 4', '...', 'Location', 'best');
grid on;


