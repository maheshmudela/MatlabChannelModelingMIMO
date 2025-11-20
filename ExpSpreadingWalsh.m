


pkg load communications;
pkg load signal;

clear;
close all;
lenStream = 100;
a = randi([1 10], 1, lenStream) + j*randi([1, 10], 1, lenStream); %[-2 -2 3 6 -8];

b = randi([1 10], 1, lenStream) + j*randi([1, 10], 1, lenStream); %[-2 -2 3 6 -8];

C1 = [ -1 -1 1 1 1 -1 1];
C2 = [ 1 -1 -1 1 1 1 -1];

d1 = zeros(length(a),length(C1));
d2 = zeros(length(a),length(C1));
d = zeros(length(a),length(C1));

for i=1:length(a)

  d1(i,:) = a(i).*C1(1:7);
  d2(i,:) = b(i).*C2(1:7);
  d(i,:)  = (d1(i,:) + d2(i,:));

  disp(i);
  disp(d1(i,:));
  disp(d2(i,:));
  disp(d(i,:));

end


txest = zeros(1, length(a));
db = zeros(length(a),length(C1));

for i=1:length(a)
 db = d(i,:).*C1;
 disp("code");
 disp(db);
 disp(d(i,:));
 disp(C1);
 disp("next");
 m = sum(d(i,:).*C1);
 txest(i) = m;
 %disp(txest(i));

end;

disp("send");
disp(a);
disp("receved");
disp(txest);
figure;
plot(abs(a),'r');
hold on;
title("Send symbol");
plot(abs(txest),'b');
hold off;
title("received");










