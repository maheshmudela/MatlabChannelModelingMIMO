
%ray based model, its dopller spectrum as
clear, clf
M = 8;
fc=5e9;  % 5ghz
fs=5e4;  % 50khz, symbol freq , no of symbol per secod.

speed_kmh=120;
Ts=1/fs;  %symbl period.

v_ms= speed_kmh/3.6;
wl_m = 3e8/fc;
lamda = 3e8/v_ms;
%pdbChannel parameter setting, of each doppl
pdf = [ 0 -1  -9 -10  -15  -20 ];
t_ns = [ 0 310 710 1090 1730  2510];
BsLosAngle=0;
MsLosAngle=0;
%base station angle spread
BsAsDeg=35;  %Laplcian Pass, narrow As mean more corrletation bewettn each element
BsAoDDeg=50*ones(length(pdf));
MsAsDeg = 2;

MsAoADeg =50*ones(length(pdf));
%speed of movile
V = 120;
thetaV = pi/4;
%At nth MS , AOA from each M path ,
%BsThetaLosDeg, BsAodDeg, AsBs , MsThetaLosDeg, MsAoADeg, AsMs , M)
[BsAodthetaDeg, MsAoathetaDeg, BsPhiRad] = genPhaseN_M(BsLosAngle, BsAoDDeg, BsAsDeg, MsLosAngle, MsAoADeg, MsAsDeg, M);
MsAoathetaRad = MsAoathetaDeg*pi/180;
%coNVERT TO RADIAN
PdfLinear = zeros(length(pdf));
PdfLinear = (10 .^pdf).*0.10;


%no of MS
N=length(pdf);
%for every mobile station, with M antneta, or ariving path, for N BS anetnta,NXM

h = zeros(N,M);
%t = (0:1/fs:1);
t = (1:1:M)
for n=1:length(pdf)

   tempP = exp(j*BsPhiRad(n,:));
   tempD = exp(j*2*(pi/lamda)*V*cos((MsAoathetaRad(n, :)' - thetaV)) .*t); %page 56, 2.32
   %tempD = exp(j * 2 * (pi/lamda) * V * cos((MsAoathetaRad(n, :)' - thetaV)) .* t);

   temp = tempP*tempD;
   disp(temp);

   h(n,:) = sqrt(PdfLinear(n)/M)*sum(temp);

end


 figure;
 t = 1:1:N;
 plot(t, h(t,:));
 xlabel('index n');
 ylabel('coeff ');
 hold on;


% Below is to make imapct of colletaion matrix, , corrleation of signal
% at tx, all antenna element, at tx, how its corrleated and same as rx side
% Later this is has to be analyzed.
%This uncorelated channel model ,
%spatial correlation
%  //rho = zeros(N,N);

% //for n1=1:N

%  //for n2=1:N
%    //rho(n1, :) = cross( h(n1,:), h(n2,:); %spatial ocrelation, of ray based channel mode)
    %or use eq 3.75
% //  end
% //end

 %generate toplitz

% finali MIMO channe coeffient , with corrletaion Rx ad Tx.
%








