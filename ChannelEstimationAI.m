%This still need to be updated wrt noise and also using bessel tabluated value for
% correlation in time domain, which show doppler spectrum, wrt fd max .
% Need to make use AI for kalman estimation for aaccurate doppler imapct, considering
% imapct CFO, freq drift and doppler due to speed.


%Channel estimation, for freq selctive fading, where freq domain response of channel
% is not flat , here  estimating doppler freq and the time , period,
% Tm = 1/fm, the pilot symbol interval has to < Tm, as well symbol period
% has to be Tsymb < Tm, so that impact of fading on channel is flat through output_precision
% bw assuming STO is corrcetd and no symbol time offset due sampling time offset.

% Channel estimation is based block type pilot, dmrs mapping type A

% X,Y , E{X.conj(Y), is correlation
% channelcoeff(tou) = a(tou)*exp(-j2pi*f*tou)*exp(j2pi*fd*t)
% theta =           -Bs.
%  v            -
%vv          - (theta)
%          Ms-----------
% doppler spectrum is freq domain of bessel function
% S(f) = fft(bessel(delt)
% u shape spectrum between -fd to fd, is kae model.
% Impulse response of channel , and PDF, impulse, show
% power distribution at diffrenet delay and its fourier trandform
% gives and its element dot product gives Autocorrleation in freq domain.

% or dopler, wrt doppler shft and its power pattren , and it is simulated
% through jake model..  and its time domain is bessel function and autocorelation
% in time domain.
% Freq domain power pattren that ipact how doppler channes power pattren, is auto corrleation,
% buts it time domain ifft , gives
% auto correlation in time domain, impact of doppler in time domain. Same appliaction
% power pattren in time domain auto corrleation, and we want to study how
% each freq elemnet is imapct is given fft of this , freq domain auto corrleation.



pkg load communications;
pkg load signal;

#pilotX = eye([0 1])
#pilotX
% Y=HX + N, this is mathemtica modelling.
% H = (Y-N)X-1
%H' = H.W
%
%H is ls estimate.?
Nfft = 64;
%H-H' = (Y-N)X-1-HW
snrdb = -20;
snr = 10^(snrdb*(0.1));

Nps = 2 ;% pilot spacing
Np = Nfft/Nps;
 %Nps = PilotInterval = Nfft/Nps;
k1 = 1:Nps:Nfft; %index location
k =  1:1:Nfft;
%kps = Nps*i +1:
%bits = Nfft*2
M= 4;
bpsk=2;
x = randi([0 1], (Nfft*2),1);

y = randi([0 1], (Nfft*2),1);

TxdataX= reshape(x(:,1), bpsk, [])';

dataSymbX =  bi2de(TxdataX, 'left-msb');

%// no need for padding
%grey code, as nebouring constetaion differ by 1 bits and lenth in bits not integer
dataX = qammod(dataSymbX,M);
X = fft(dataX,Nfft);

RxDataY=reshape(y(:,1),bpsk,[])';

dataSymbY =  bi2de(RxDataY, 'left-msb');
dataY = qammod(dataSymbY,M);
Y = fft(dataY,Nfft);

figure
plot(real(dataSymbX), imag(dataSymbX), '0', 'MarkerSize', 10, 'LineWidth', 2);
title('after QAM');

%X = randi([-1 1], Nfft,2) + j*randi([0 1],Nfft,2);

%H-H' = (Y-N)X-1-HW
%Y = randi([-1 1], Nfft,2) + j*randi([0 1],Nfft,2);

%Assuming Ls estimate is available.
%H = randi([0 1], Nfft,2) + j*randi([0 1],Nfft,2); %col matrix

H = Y(k1,1) ./X(k1,1);

h = ifft(H,Nfft);

%H-H' = (Y-N)X-1-HW
%X = randi([0 1], Nfft,2) + j*randi([0 1],Nfft,2);

%repmat(A, M, N): This function creates a large matrix by replicating (tiling) the input
% array A M times along the row dimension (vertically) and N times
% along the column dimension (horizontally).

% W = Rhh~ *(inverse(Rhh + (varinace_noise/varinace_signal) Unity.
%H_estimate = W*H~ ( H~= H+X-1z)

%Rhh  = here same freq , differnet time , imapct need to analyse
%that is ifft of doppler respose

%Rhh~ =freq domain corrleation, h and h~, freq domain reponse of power delay profile

 % h in CIR, channel impulse reposen
 k1 = 0:length(h)-1;

 %energy
 %h=Nx1, h'=complex conjugate transpose 1xN
 %NxN
 %spatial channel corrleation, matraix, channel covraince matrx
 %Rhh
 hh = h.*conj(h') % h', h hermitian, h transpose of conjugate. energy...

 %computing rms delay , which require to compute the first moment of pdf,wrt
 %delay, hh(0)*k1,
 % mean= sum(hh*k)/sum(hh);
 %second moment= sum(hh)*sqr(k)/sum(hh);
 sqrk1=square(k1);
 %rms delay = sqrt( secondMoment-sqr(first moment);
 FirstMoment  = sum(hh.*k1)/sum(hh);
 secondMoment = sum(hh.*sqrk1)/sum(hh);
 rmsDelay = sqrt(secondMoment-square(FirstMoment));

 df=1/Nfft; % resolution or ,
 rf=zeros(Nfft,Nfft);
 DnProfle= j*2*pi*rmsDelay*df;
 for i=1:Nfft
   for j=1:Nfft
     rf(i,j)=1/ (1 + (i-j)*DnProfle);

   end
 end





 % element wise complex multipication (x+iy) and x-jy)= variance
 % [ sqr(h1)*k1   sqr(h2)*k2   sqr(h3)*k3    ]

 % tmp = (h.*conj(h)) .*k1;
 %r = sum(tmp)/hh;



 % ' dot, perform conjugate transpose and .', non conjugate transpose
 % [h1  h2  h3] *[k1;k2;k3]= scaler, dot product/hh

 %r2 = k1*tmp.' /hh;

 %to claculae rms delay from channel estimateCfo
 % CIr = ifft(h),
 %P(tou)=,is how channel impulse resposne, wrt time varies
 % mag sqr is power delay profile.
 %H_cap = X-1.Y, P(tou) = [ sqr(h1) sqr(h2).] tou is ith tap i*tS
 %mean delay is sum(P(ith tou)*ith tou/sum(P(ith tou)
 %and rms delay is sqrt( sum( ( h1*h1 ...)*sqr( tou - mean)

 %tou_rms = sqrt(r2-r.^2);

 %Rhh~= as function of k1, = 1/(1 + j*2*pi*touRms*k1*df)
 %df=1/Nfft; % resolution or ,
 %DnProfle= j*2*pi*tou_rms*df;
 %REPEAT the array, A, r=1, and col Np
 % 1;2;3;  1;2;3;  so repeat the col by Np
 %K1 = repmat([0:Nfft-1].',1,Np);
 %K2 = repmat([0:Np-1],Nfft,1);
 %rf=zeros(Nfft,Nfft);
 %rf = 1 ./(1+DnProfle*Nps*(K1-K2));
 % creating toeplix
 % repeat [1,2..Np), row wise just once but repeat col wise Np
 K3 = repmat((0:Np-1).',1,Np);
 K4 = repmat([0:Np-1],Np,1);

 %sort of topltix matrix of t2-t1; t is
 K = K3-K4; %NpxNp;
 %1:Nps:Nfft is NP, so Np
 tmpDn = 1+ DnProfle*(K3-K4);
 rf2 = inv(tmpDn);
 Rhp = rf;
 % create MXN matrix with ones in main diagnol and zero else
 % eye(2,4), 1  0  0   0
 %           0  1  0   0
 Rpp = rf2 + (eye(Np,Np)/snr).';

 %InvRpp = inv(Rpp)*h.'

 newCord = linspace(1,Np,Nfft);
 [Xq,Yq] = meshgrid(newCord,newCord);
 Rpp1 = interp2(1:Np, 1:Np, Rpp, Xq,Yq,'Linear');

 h_mmse = transpose(Rhp*inv(Rpp1)*h);
 figure;
 plot(abs(fft(h_mmse,Nfft)));
 title('MMSE');




