%Moose appraoch for cfo estimation using traing with SRS, with comb like pilot signal
%Moose apparoch , convert tow cosecutive symbol preamble and take fft of both
% Y1 AND Y2 as fft, angle(sum((Y1 .* Y2))/2pi, its complex conjugate , that is inneer product
% angle(sum(Cp1 .* Cp2)), cp and cp2 as time doain complex symbol as Cp based appraoch
%Classen:  for dmrs based where.



pkg load communications;
pkg load signal;

clear, clf
cfo = 0.25; % normlized to samling freq;

%dmrs is
Bw=100e6;

fc = 2.5e9;
scs = 30e3;% 3
% x. conj(x')), is =(sqr(a) + sqr(b)), modulus of complex no
%x = a+jb
%x' = a-jb
% x. *conj(x') =
% 4 qpsk, power(2,2) 16= power(2, 4) qam,64 qam ( pow(2, 6), 256 = power(2,8)
M=4;
qpsk=log2(M); %// POSSIBLE consteation possible per 2 bits
Nfft = power(2, nextpow2(Bw/scs));

Tg = power(2,3)*16;
SymbLen = Nfft + Tg;

bitsLenth = Nfft*qpsk; %/ every Re HAS 2 BITS, as per qpsk

NoSymbol = 12; %grid, each colm is symbol
% encoded bit interlaved, scrambled data bits.
databits= randi([0 1],bitsLenth, 12); %/col as bits

%for every symbol
for i=1:12

  disp(['symbol ' num2str(i)] );

  dataMatrix = reshape(databits(:,i), qpsk, [])';

  dataSymb =  bi2de(dataMatrix, 'left-msb');

  %// no need for padding
  %grey code, as nebouring constetaion differ by 1 bits and lenth in bits not integer
  modsig = qammod(dataSymb,M);

  figure
  plot(real(modsig), imag(modsig), '0', 'MarkerSize', 10, 'LineWidth', 2);
  title('after QAM');

  %modsig = qammod(databits, M, 'InputType', 'bit');

  %disp(length(modsig));

  if length(modsig) ~= Nfft

       paddlen = Nfft - length(modsig);

       paddSig = [modsig;   zeros(paddlen,1)];

   else
       paddSig = modsig;

end

    qamTimeDomain = ifft(paddSig,Nfft);

    figure;

    plot(abs(qamTimeDomain));
    %hold on

    y_cp = zeros((Tg+Nfft), 1);

    y_cp(1:Tg,1) = qamTimeDomain(Nfft-Tg+1:Nfft,1);

    y_cp(Tg:Nfft+Tg-1,1) = qamTimeDomain(1:Nfft,1);

    figure;
    plot(mag2db(abs(y_cp(1:Tg))),'r');
    hold on;
    plot(angle(y_cp(1:Tg)),'r--');
    title('TX cp1');
    hold on;
    plot(mag2db(abs(y_cp(Nfft:Nfft + Tg))),'b');
    hold on;
    plot(angle(y_cp(Nfft:Nfft+Tg)),'b--');

    title('TX CP2');
    hold off;

    %y_cp', is conjugate transpose.
    y_mag = sqrt(y_cp .* y_cp');

    %figure;

%plot(y_mag, 'b', 'LineWidth', 2);
%hold off;
%plot(y_mag);

   power = (y_cp .* y_cp')/length(y_cp);

   powerDb = 10*log10(power);
   %1/fc= mean multiple Ts, sampling instance.
   t = (0:1:(length(y_cp)-1))/fc;

    %Updample.
    localOsc = exp(j*2*pi*fc*t);

    yUp  = y_cp .* localOsc.';

    snrdb = randi([ -20 20], 1,1);

    %gain = 10^(snrdb(i)/10);
    %yUp = gain*yUp;
    yawgn = awgn(yUp,snrdb);

    localOscilator=exp(-j*2*pi*fc*t);
    %downconversion.
    yawgnD = yawgn .* localOscilator.';
    nn = 1:1:length(yUp);
    cfoadd = exp(j*2*pi*cfo*nn);
    %// add cfo
    % elemnt wise cfo addition.
    yawgnCfo = yawgnD .* cfoadd.';

    %estimate using CP;
    Cp1 = zeros(1, Tg);
    Cp1 = yawgnCfo(1:Tg);

    Cp2 = zeros(1,Tg);
    Cp2 = yawgnCfo((Nfft:Nfft + Tg -1));

    figure;
    plot(mag2db(abs(Cp1)),'g');
    hold on;
    plot(angle(Cp1),'g--');
    title('Cp1');
    hold on;
    plot(mag2db(abs(Cp2)),'r');
    hold on;
    plot(angle(Cp2),'r--');
    title('Cp2');
    hold off;


    % this is normaized radian wrt fs, so actual fm is theta;
    % inner dot product = sum(elemtwise prdict)
    %AUTO CORRELATION.
    avgRadian = angle(sum(Cp1 .* Cp2));
    %fmn = fm/fs;
    %fmn = avgRadia8;
    %fm = avgRadian
    fm = (1/(2*pi))*avgRadian;
    disp(fm);
end



































% extract every



