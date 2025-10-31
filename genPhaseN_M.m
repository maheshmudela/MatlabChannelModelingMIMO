

%Phase NTxMR, for given n [ Nt, for each n  to M, Path to any Ms, as rx diversity
%for given path, AOD
%each AOD, from base station mean how ms, with M path is targeted, NXM
function [ BS_PHI_rad, BS_thetaOffsets, MS_thetaOffsets ] = genPhase(BsThetaLosDeg, BsAodDeg, AsBs , MsThetaLosDeg, MsAoADeg, AsMs , M)


%N=8
%BsThetaLosDeg = randi([0,360],1,N);
%BsAodDeg      = randi([0,360],1,N);
%AsBs      = 2;
%MsThetaLosDeg = randi([0,360],1,N);
%MsAoADeg      = randi([0,360],1,N);
%AsMs      = 2;
%M = 2;

  %if nargin==6
  %  M = 20;
  %end

  %for every n path to Mr, SAY n=1 is one MS, at AOD, theta as LOS,
  % here for every Ms, each MS, is with M recived antenta
  % each with given AS, angle of spread , 1 deg AS mean high correlation and beamforming
  % is higly confined, each M path with equal power , but each path will have random phasez

  N = length(BsAodDeg);
  %Phi is wrt each nth mobile station, but offset generated as per all M path to that MS
  BS_PHI_rad =  2*pi*rand(N,M); %cos(theta + phi), so possibe phase value.
  BS_thetaOffsets  = computeOffset((BsThetaLosDeg + BsAodDeg), AsBs);
  MS_thetaOffsets  = computeOffset((MsThetaLosDeg + MsAoADeg), AsMs );

  BS_ThetaDeg = zeros(N, M);
  MS_ThetaDeg = zeros(N, M);
  %RANDOM ORDER of phi

  N = size(MS_thetaOffsets,1);
  for n=1:N
    index = randperm(M);

    %BS_ThetaDeg(n,:) = BS_thetaOffsets(n,index);%*(180/pi);
    MS_ThetaDeg(n,:) = MS_thetaOffsets(n,index);% + BS_PHI_rad(n,:);%.*(180/pi);

  end

 figure;
 t = 1:1:N;
 plot(t, MS_ThetaDeg(t,:));
 xlabel('index n');
 ylabel('angle deg');
 hold on;
end




