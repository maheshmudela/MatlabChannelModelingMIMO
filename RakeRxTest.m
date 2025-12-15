
%CDMA2000, Under progress:
%Objective is to simuate signal channel reverse link with one traffic channel
% and pilot  and implmenet rake recver and channel esimation for L path SIMO
% Single user and based statation as L antenna Rx diversity at 1.2288e6
% RC 3, and further can be enahnced.


%https://helpfiles.keysight.com/csg/n7601/Content/Main/About_cdma2000.htm

%For reveserse link slot : 12 + 18
%All piot bits 111 = 1-2.1=-1 , i no , every bit is one symbl
% s(k,bit(i)) = 1-2*bit(i) + j 0, pilot has only I +j 0.
%Bpsk , qpsk, 9psk, , there is  s(k,bit(i,i+1,i+2) =  i(k) + j*q(k), where k = 0<k<DataBits/3..
%Reverse Traffic Channel with Radio Configurations 3 to 6
% 8 psk, reverse traffic use bpsk, qpsk and 8psk fllowed by orthogonal spreadingTechniq
% and piot it use orthogonal spreading witout any modulation

%But for other RC 1/2 ,reveser link use orthogonal modulation for traffic channel
% for this first 6 bits mapped to one of 64 orthogonal code. So no bpsk/qpsk
% or 8psk.. all or every 6 bits are grouped .

% The Reverse Packet Data Channel with Radio Configuration 7 shall use the W24 Walsh
% function for encoder packet sizes of 192, 408, 792, and 1560 bits. The Reverse Packet Data
% Channel shall use the W12 Walsh function for an encoder packet size of 3096 bits. The
% Reverse Packet Data Channel shall use the W12 and W24 Walsh functions for encoder
% packet sizes of 4632, 6168, 9240, 12312, 15384, and 18456 bits.


%For reverse pilot channel it use orthgonal spreading with WALS0(64);

%For reverse channel , enoder of 20ms encoded pasket with 192 bits
% at 9.6kbps , for hgher encoder size we can other walsh code.
% when encoder generate 9.6 kbps at 20ms enocoder size,
% it walsh1(2),wlash2(4)= H(2,:)
% 192/16 = per slot (1/2 rate) = 12 pair, each symbol qpsk modulated 2 bits per
% 12*2 = 24 chips and 1356 chips pr slot, mean 1152 chips pilot and
%384 next segment for traffic and power controll, last chip is punctured or overwriiten
% if power control bit is inserted at that chip location.
% what is spreading factor , lets say it using Wals1(2), then repeatation
% needed 12 symbol , 384/2= 192 qpsk but yu has 12 as actual so repeatation factor
% 192/12=16, so each symbol is repeated 16 times, to 192 symbol qpsk
% and finally walsh spread , and then finally complex long code scrammbled modulo 2 addition
% and finally I and Q split and


%Walsh function Wn
%N represents the n-t 1 h Walsh function (n = 0 to N – 1) of length N with the
%binary symbols of the Walsh function mapped to ±1 symbols using the mapping '0' to +1
%and '1' to –1. A binary-symbol Walsh function of length N can be serially constructed from
%the n-th row of an N × N Hadamard matrix with the zeroth row being Walsh function 0, the
%first row being Walsh function 1, etc. Within Walsh function n, Walsh chips shall be
%transmitted serially from the n-th row from left to right
% Now 1.2288e6, chips

%Pilot bits is continous stream of chips n forward link and time multiplex
% symbol on reverser link.
%But wcdma/umts has fixed on of pilot bits.

%Forward link pilot bits are transmitted continously, but orthogonally spread at
% at 1.2288e6 chips per second at osf 64, so
% bits per frame   1.2288e6/64 = 19.2kbs or per frame 384 bits

%In forward link where pilot bits are sent contously its time muliplex, 20ms frame, 16 slot, each 1.25 ms
% 384 bits per 20ms, 1.25ms has 24 bits in 1.25ms, but some bits are
% power cntrl bits , so 18 bits are pilot bits and rest bits are power
%so 18 bits obsvf = 18*64 = 1152 bits per 1.25  and per frame of 20ms
%ttal piot bits 16*18 = 288 bits.
% (288 + 96 )*64 , is in 20ms, so 1sec has 1.2288e6 bits


%Rake reciver BTS, processing reverse path, uplink single user , multiple
% Antena, L. 1xL , L receive path.

% Say Getting pilot ,
% Say Transmiited pilot bits and traffic bits
% Data bits and pilot bits in ds-spread spectrum , pilot bits and data bits
% are sent at the same but orthgonally spread with diffrent orthigonal codes
%The type of modulation depends on the encoder packet size,
%suppose encoder size increased as we need more bits to map to symbols
% so 4 psk, group of 4 bits on one symbl
% encder size/2, qpsk, 2 bits on one sybol, so that sybol on

%The Reverse Pilot Channel data shall be spread with W064 as specified in 2.1.3.1.13.


%for for Reverse link to avoid non lineary issue and battery power limitataion,
% it go for oqpsk, offsetQpsk, power saving modulation.
%so every 2 bits mapped to some phase but then limited/offset by, could
%have some mapping that can handle such power requirement or fadding and interfernce
%noise, doppler ?

%assuming mapped 2 bits per symbol so per secon no of symbol is same as data
% rate.. or we can have 4 bits/6 bits  so no of symbol per second
% depened dataRate and modulation scheme.


%For reverse link, data and pilot bits are sent in time muliplex way
% 1.2288e6 chips per second, per slot 1536 chips per slot, is shared by
% data symbol and pilot symbol.
% At 9.6 bist per second, 9.6 bits *1.25=12 bits per slot, 1/2 rate, 24 bits
% bits per slot or it means 1536 chips, now so every slot is not for data
% so seems may be 1152 for pilot and 384 for power control chips
% say 1
% 1356 chips pr slot, mean 1152 chips pilot and
%384 next segment for traffic and power controll, last chip is punctured or overwriiten

dataRate = 9.6e3 ;% kb per second (AMR EVS codec, support variable Rate )
%frame size 20ms, so per frame , data bits  in 1ms is 9.6 bits, 10ms ==>96bits
%so how many frames for 1 second, 100 frames
dataBits = 192; %bits per 20ms
chipRate = 1.2288e6;
Rate = 1/2; %1/2 is for reverse and forward channel and 1/3 for cch control
%Assuminf we perfroming 1/2 conv encoder with contstrained lenth 9, 8 stage
% shift register..just like modulo 2 GF Polynomial, but every time state changes
% and also for higher rates are supported with 1/3,turbo code (RSC) with
%g0= 753 octal and g1=561 octal
% 19200 bits per second, each 2 bits is mapped to qpsk I/Q symbol as per spec.
DataPerSec = dataRate/Rate;

%per slot 1536 chips, some bits with walsh2 spread and some bits with wasl4,
% wals4 4 means high protection resulsting chips/ symbol is fixed.

% then apply gain and any precoding impact??
%dataRate=9600 bits OR 2 bit group symbol ,
symbolPerSecond   =  randi([0 4], 1, dataRate) +j*randi([0 4], 1, dataRate);

% 9.6 symbol per ms, each slot is 1.25ms, so 12 symbols per slot, each slot holds
% Pilot chips per slot is 1152 and is spread by W0(64), so 18 bits per slot
% 18 bits per slot , 18*16=288 and per 20ms, 50 frames  , so 288*50 bits per second
% these are fixed.

pilotBits = 288*50;% per frame, 20ms, and  no of frame pre secon 50 frames
                   % per frame 16 slot of each 1.25ms, 6 bits per slot for power.

pilot      =  randi([1 1], 1, pilotBits );

pilotSymPerSec = 1-2*pilot(1:end);

osvfP = 64;
osvfD = 128;

Hp    = hadamard(osvfP);
Hd    = hadamard(osvfD);

%Channelisation code
W0   = Hp(1,:);
W1   = Hd(1,:);
orthogonalData=zeros(1, length(symbolPerSecond));
orthogonalPilot=zeros(1,length(symbolPerSecond));

%Generate LONG code, BASED ON GF POLY WITH??
%Each PN chip of the long code shall be generated by the modulo-2 inner product of a 42-bit
% mask and the 42-bit state vector of the sequence generator
%the mobile station shall use one of the following two long code
%masks unique to each channel: a public long code mask or a private long code mask.
LongCodeComplexP =[];



%The I long code for Spreading Rate 1 shall be the long code sequence specified in 2.1.3.1.16.
%The I long code for Spreading Rate 1 shall have a chip rate of 1.2288 MHz. The Q long code
%for Spreading Rate 1 shall be the I long code delayed by one chip

%The mask used for generating the I long code for Spreading Rate 1 (or equivalently, the first
%component sequence of the I long code for Spreading Rate 3) varies depending on the
%channel type on which the mobile station is transmitting. See Figure 2.1.3.1.16
UserLongMask = 0x8181;

%Generate 1.2288e6 chips per sedecond
ILongCodeUser = longCodeGnerator(UserLongMask);
% One chip delay, I chip delay mean 0QPSK? MODULATION at final stage.
%The final OQPSK signal(S_{OQPSK}(t)) is the sum of these two components: 
%\(S_{OQPSK}(t)=I(t)cos (omega _{c}t)+ Q (t-T_{b})sin (wc*t)\)(Note:
% \(T_{b}\) is the bit period, which equals \(T_{s}/2\), the symbol half-period delay.) 
QLongCodeUser = [zero(1,1) ILongCodeUser(1:end-1)];

%Generate 1.2288e6 chips per sedecond
%The I and Q PN sequences used for quadrature spreading shall be as specified in
% 2.1.3.1.17.1 and 2.1.3.1.17.2. These sequences are periodic with a period of 215 chips for
% Spreading Rate 1
% Spreading Rate 1 The PN sequences shall be based upon the following characteristic polynomials:
% PI(x) = x15 + x13 + x9 + x8 + x7 + x5 + 1
% (for the in-phase (I) sequence)
% and
% PQ(x) = x15 + x12 + x11+ x10 + x6 + x5 + x4 + x3 + 1
% (for the quadrature-phase (Q) sequence).
IchannelShortCode = IchannelShortCodeGenerator();
QchannelShortCode = QchannelShortCodeGenerator();

ILongCodeComplexD =[];
QLongCodeComplexD =[];


%This whole will work per frame of encoded audio with 20ms? Call for every frame
pilotSymPerFrame20ms=[]
orthogonalDataPerFrame20=[]
PilotSymbolPerSlot=18; %bpsk modulated
DataSymbolPerSlot=12;
pbp =1; % 1 symbol is punctured, for on/off of or near or far.
powerSymbol=one(1,1,pbp);
RepeatationFactor = 30;

%This will process slot by slot every 20ms frame encoded audio
for slotIndex=1:16

    %Crc/Convolution enoding 1/2, interleaving done  is already assumed to be done

    %Prepare for pilot symbol , i +j 0, all symbol bits mapped to I
    pilotSymbolPerSlot=[pilotSymPerFrame20ms((slotIndex-1)*PilotSymbolPerSlot+1:slotIndex*PilotSymbolPerSlot)];

    %Osvf spread , channelization , to make 1152 chips per slot
    pilotSymbolPerSlotSpread =remat(pilotSymbolPerSlot,1,Hp);

    %Scrambing the pilotSymbolPerSlotSprea,1152
    pilotSymbolPerSlotScrambled = LongCodeComplexP.*pilotSymbolPerSlotSpread;

    %As per data 9.6kbps, with 9.6 symbol per ms, each slot is 1.25ms, so 12 symbols per slot, qpsk modulated
    % In currnet slot , current symbol are spread by W1(2)
    DataSymbolPerSlot = [orthogonalDataPerFrame20((slotIndex-1)*DataSymbolPerSlot+1:slotIndex*DataSymbolPerSlot)];

    % 24 symbol after spread but 384-24=360 space is free, so repeatatation of symbol
    % by 30 time, to match the chip rate 1.2288e6 per second.
    DataSymbolPerSlotRepeat = repmat(DataSymbolPerSlot, 1, RepeatationFactor);

    %Walsh spread channelization code by W1(2);
    %This Walsh function repetition factor is the number of Walsh function sequence repetitions per
    %interleaver output symbol.
    DataSymbolPerSlotSpread = DataSymbolPerSlotRepeat.*[ Hd(1:end)];

    %Scrambling with LONG CODE, Direct sequence spreading using the long code
    %The in-phase spreading sequence shall be formed by a modulo-2 addition of the I-channel
    %PN sequence and the I long code sequence.
    ILongCodeComplexD = modulo(ILongCodeUser + IchannelShortCode,2);

    %The quadrature-phase spreading sequence shall
    %be formed by the modulo-2 addition of the following three terms: the 2
    %W1 Walsh function,
    %the modulo-2 addition of the I-channel PN sequence and the I long code sequence, and the
    %decimated-by-2 output of the modulo-2 addition of the Q-channel PN sequence and the Q
    %long code sequence. The decimator shall provide an output that is constant for the two
    %chips corresponding to the two symbols of the 2 W1 Walsh function, and the value of the
    %decimator output for the 2 W1 Walsh function period shall be equal to the first of the two
    %symbols into the decimator in that period. The 2 W1 Walsh function time alignment shall be
    %such that the first Walsh chip begins at the first chip of a frame.
    QLongCodeComplexD = QLongCodeUser .*QchannelShortCode;

    QLongCodeComplexD_Dec2= [QLongCodeComplexD(1:2:end)]

    %Walsh spread by 2;
    QLongCodeComplexD_Spread = QLongCodeComplexD_Dec2.*H1(2);

    %SHORT CODE QUDARTURE SPREADING WITH WALSH2.
    LongCodeComplexD = QLongCodeComplexD_Spread.*ILongCodeComplexD;

    DataSymbolPerSlotScrambled = DataSymbolPerSlotSpread*LongCodeComplexD;

     %If there are mutplie SCH OR MORE TRAFFIC channel then , but now its one
    % DataSymbolPerSlotScrambled = DataSymbolPerSlotSpreadChanel2*LongCodeComplexD;

    % Fina data DataSymbolPerSlotScrambled =  DataSymbolPerSlotSpreadChanel1+DataSymbolPerSlotSpreadChanel2
    % so on, all data channel with pilot and power control bits are scambled with same scarmbing code.

    %Time multiplex   1152+384 =1156 slot with pilot and data bits
    FinalSlotData = [pilotSymbolPerSlotScrambled DataSymbolPerSlotScrambled];
    powerSymbol = [];

    %Puncturing the data bits for adding power control bits
    %Rc 3, 800 bps that 1 bits
    FinalSlotDataPunctured=[FinalSlotData[1:end-pbp powerSymbol[1:pbp]];

    %Split I/q
     %1156 InPhase 1156 QuadPhase
    [Idata  Qdata] = splitIQ(FinalSlotDataPunctured);

    IdataGain = Idata *gain; %1.2288e6 = fs, and lpf, has t0 around < 614 khz
    QdataGain = Qdata *gain; % 16bit I , 16 bit Q, DAC

    %Fow now lpf/baseband filtering is not impmeneted into this simulation.
    %I and Q are used to modulate the carrier  0QPSK

end

% When we hav channel h(k), , we can have
% 1. we can break the channel svd decompistion, to have procoding
% and channel matrix =  H = UZV, Z, we can have rank of this matrix
% based on rank, we can most L path to be used for rake reciver
% or also causing predistortion at sent so that such channel
% is made flat for all w /Ts,

% This 16 slot 20ms frame is channel

% Add channel artificat/fading/awgn and fading


% Rake Receiver at base station

%
%The received input is downconverted to based band using inpahse and quadrature oscillator

% In phase and  Quadrature

IRxData =  [IdataGain(1:end)];

QRxData =  [QdataGain(1:end)];

%  match filtering for frame sync
% FOR ALL slot in 20ms
IRxDataS = [IRxData[1:(2*Ts/Tc))]
for ci=0: 2*Osvf;
     Eng = E{IRxData .*PnSeqI};

     if Max < Eng
        Max = Eng
        offset = ci;
     end
      IRxDataS = [zeros(1,ci) IRxDataS(ci:end]


end
%Make sure Rxdata buffer is atleast holdng 40ms data, and read every
% 20ms from the buffer at offset ci, and from appliaction, this 40ms buffer is always popluated
% at every 20ms.
RxCorrected = [IRxData[ci:end)]



%long code for user lave mask? mobile term id?

%Later I and Q are scrambled  with short m-sequence, the short code
% diffrenetiate teh signal emiited by diffrenet base station,
%  All base station uses same code but with differnet timing offset and
% recever attemts to syncronize to this short code  strongnest corrleation
% idenetifies the nearest base stataion


%Gpss receiver seraches for presence of gpss satellite signal
% by correlting the receved signal with diffrent gold sequence
% since period, is 20ms, same signal seq is repeated in period
% of 1 ms.

%Since the current slot cntain orthogonal data from all user
%           1152   384 chips          1556
%   slot = [pilot dataUser1] (+) [user2 data ]  (+)
% so same slot contain all the user but orthogonally , This means better orthgonallygnallity
% or long orthogonal code give more possibliy of more users!

% So first tracking recviver for coode sync, frame sync


%This is recever specifucally as base station , that
% has data from all users ms.
%Capture 20ms frame and process in ping pong way

%prev 20ms and current 20ms , 40ms buffer
%
bufferIQ = [IQ20msBuffer1  IQ20msBuffer2]

[basebandI basebandQ]= HardwareDownconversion(bufferIQ);


% Generate I AND q short code PN seq PN seq of length 15, bases on currnet basestation
%

[pnSeqI  pnSeqQ] = generatePN(stationId);

[basebandI basebandQ] = HardwareDownconversion(bufferIQ);

% Later stages it will be analysed  for data . error
[descsrambledI]  = basebandI .*pnSeqI;
[descsrambledIQ] = basebandI .*pnSeqQ;


%Need to anayze more on this channel estmation based Pilot
% hbar = [hbar0 ...];
%  l=0 [ a0  a1  a(L-1     ]  [h0] = [hbar0]
   l=1                        [h1] =

% estimated pilot = hbar . *yRx ; after sommtehung of h
%
% Here y = [], as recved descrambled  demultiplex pilot, from  l=0 path
% and same form other path l=2



%Fir LPF for smotthening the decsison of channel coeffent of lth path
% let say fir coeff , Need to check ,
% Adapative fir filter, the no of coeffient are based on max delay supported
% a= R
%h = SmoothingChannel(h(l,:);

%LMS UPDATE OF CHANNEL, in Rake reciver . FIR as L tapped delay , time varinat
% channel and its performance.
% It somethinging like passing or recveing the symbols with persion Ts, all consectuive symbols
% passed from such tapped delay lne channel.. as

% symbol rate or Ts, is limited by time variant channel ,
% For all pilot symbols in slot.
% this channel estimation.. in this loop based on rake rever apparoch
%




for step=1:M % M now can be seen as no of symbols in slot?

  %L say Ts , Ts/L=   or 2pi/L, so, N or L is some sort
  % sort L PATHS OR , What i can say, more L mean better SNR for any modulation
  % N or L can be choice of modulation /SNR condiation, better or SNR, LESS
  % FADING slow fading, no fading we can have L as say 1 or 2, so
  %  less perfomrance power needed, so it means L has to be function of
  % SNR as recved and analyzed based on feedback path.
  % For every next symbol this FIR filter gain and updated to accuracy
  % finally for any time < Td, time dipersion
  % from  awgn estimate of channel, we get dpller spectrum, or time vairant
  % correlation pattren.. and found td/fm,
  % Ts <td, ideally .. IF not, then Ts/Td, where there is channges
  % n of L or N =Ts/td, this ddteremine resoltion, phase
  % channeged as info  or phase chnages as info < 2*pi/N, so
  % what psk it can support is given by <= 2*Pi/N, witout any phase error
  % same for fdm, ??
  %Ts/Td, td time dispersion, FIR order is Ts/Td,?
  for l:1 L-1

      %Descrambling wrt station code,
      %hk = conjugate(bk) . *ypk  =  *bk .*(ak*bk + wk) = *bk.*bk*ak + *bk.*wk where |bk*bk| = 1
      %  = ak + noisek, ak is channel coeffient?? which is cmplex??
      % Genrate pilot symbol, walsh sperad and then scmable with short code
      % and then same
      % first cut hk, is found to be for any path.., generally pilot
      % signal are speard with only short code .

      % all delay in colmwise , l col by adding zero, as shift
      descsrambledI(l,1) = sum([ zeros(1,l) basebandI] .*pnSeqI); %Tc/Ts=N
      descsrambledQ(l,1) = sum([ zeros(1,l) basebandQ] .*pnSeqQ);


      rateFactorI = prvSymbolI(l,1) *errorI ;
      rateFactorQ = prvSymbolQ(l,1) *errorQ ;

      %lms fir filter that smmthen the h, lth path channel estimate
       hsmothenI(l) = lpf(rateFactorI, hi(l,1)) ;
       hsmothenQ(l) = lpf(rateFactorQ, hq(l,1)) ;

       estimatedPilotSymbol(l) = descsrambledI(L,1)*hsmothenI(l,1);
       estimatedPilotSymboQ(l) = descsrambledI(L,1)*hsmothenI(l,1);

       TempErrorI(l) = ActualPilotSymbol - estimatedPilotSymbol(l)
       TempErrorI(l) = ActualPilotSymbol - estimatedPilotSymbol(l)

     end

    %fINAL ESTIMATED i SYMBOL FROM ALL PATH

    errorI = sum(TempErrorI)
    errorQ = sum(TempErrorI)

end

for user=1: K
   %Now based on abve channel estimate,

   %Dequalize the data in same slot  and as well other orthignal trafiic users

end









  end:






end


for userIdx=0: MaxUser

     scrambler = scrambler[userIdx,:]






end



























