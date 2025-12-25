

pkg load communications;
pkg load signal;

clc; clear all; close all;

%CDMA2000, Under progress:3GPP2 C.S0002-D

%Objective is to simuate signal channel reverse link with one traffic channel
% and pilot  and implmenet rake recver and channel esimation for L path SIMO
% Single user and based statation as L antenna Rx diversity at 1.2288e6
% RC 3, and further can be enahnced.

% Capture 1 second data, break into 50 frames of 20ms, and each frame with 16 slot of 1.25ms


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
frameCount=0;
dataRate = 9.6e3 ;% kb per second (AMR EVS codec, support variable Rate )
%frame size 20ms, so per frame , data bits  in 1ms is 9.6 bits, 10ms ==>96bits
%so how many frames for 1 second, 100 frames
dataBits = 192; %bits per 20ms
chipRate = 1.2288e6;
osf = chipRate/dataRate;

Rate = 1/2; %1/2 is for reverse and forward channel and 1/3 for cch control
%Assuminf we perfroming 1/2 conv encoder with contstrained lenth 9, 8 stage
% shift register..just like modulo 2 GF Polynomial, but every time state changes
% and also for higher rates are supported with 1/3,turbo code (RSC) with
%g0= 753 octal and g1=561 octal
% 19200 bits per second, each 2 bits is mapped to qpsk I/Q symbol as per spec.
DataPerSec = dataRate/Rate;
%Every 2 bit , every 4 or every 6 or 9 bits mapped to complex constellation
%ppoint , with I AND Q

%Pilot channel in revser , bpsk with orthogonal spreading , same for traffic channel

%transmitting on the Reverse Packet Data Channel, the mobile station uses BPSK
%modulation, QPSK modulation, or 8-PSK modulation followed by orthogonal spreading as
%specified in 2.1.3.1.11.


%per slot 1536 chips, some bits with walsh2 spread and some bits with wasl4,
% wals4 4 means high protection resulsting chips/ symbol is fixed.

% then apply gain and any precoding impact??
%dataRate=9600 bits OR 2 bit group symbol ,

%This data bits are time multiplex with pilot bits in Pich?.
DataRateP = 4800; %


DataBitsP  =  randi([0 1], 1, DataRateP);

%19.2
dschBits  =  randi([0 1], 1, dataRate);

%9.6*1.25=12, 12*16=192 bits per 20ms, x50, frames , mean 9.6kbs.
DataBitsPerSecondP   = 1 - 2*DataBitsP;


% 16 psk, every 4 bits mapped to iq
dschBitsSymbol = reshape(dschBits, 4, []);

%
%Each 4 bits set generate one index  that points to array of some complex no.
%dschBitsSymbolMapped = MappingFun(dschBitsSymbol);

%For Now mapped to consrellation point, for now its just random
DataBitsSymbolMappedP = [ DataBitsPerSecondP] ;% randi([-15 16], 1, length(DataBitsPerSecond)) + j*randi([-15 16], 1, length(DataBitsPerSecond));


%For Now mapped to consrellation point, for now its just random
dsch2BitsSymbolMapped = randi([-15 16], 1, length(dschBitsSymbol)) + j*randi([-15 16], 1, length(dschBitsSymbol));


% 9.6 symbol per ms, each slot is 1.25ms, so 12 symbols per slot, each slot holds
% Pilot chips per slot is 1152 and is spread by W0(64), so 18 bits per slot
% 18 bits per slot , 18*16=288 and per 20ms, 50 frames  , so 288*50 bits per second
% these are fixed.
%Pilot  per secondd?
pilotBits = (288*50);% per frame, 20ms, and  no of frame pre secon 50 frames
                     % per frame 16 slot of each 1.25ms, 6 bits per slot for power.

pilot      =  randi([0 0], 1, pilotBits );

%I + j q = i=1, q=0;, all on I stream path.
pilotSymPerSec = 1-2*pilot(1:end);

osvfP = 64;

osvfD = 128;

Hp    = hadamard(osvfP);
Hd    = hadamard(osvfD);
H1    = hadamard(2);

%Channelisation code
W0   = Hp(1,:); %Pilot

% Spreading factor : 128 W1/W2
W1   = Hd(1,:); %time Multiplex traffic with pilot
W2   = Hd(2,:); %orthogonal traffic with other channel such as pilot

%2
W3   = H1(1,:);

%orthogonalData  = zeros(1, length(dschBitsSymbolMapped)*length(W1));
%orthogonalPilot = zeros(1, length(pilotSymPerSec)*length(W0));

%Generate LONG code, BASED ON GF POLY WITH??
%Each PN chip of the long code shall be generated by the modulo-2 inner product of a 42-bit
%mask and the 42-bit state vector of the sequence generator
%the mobile station shall use one of the following two long code
%masks unique to each channel: a public long code mask or a private long code mask.
%Long long mask
LongCodeComplexP =[]

%The I long code for Spreading Rate 1 shall be the long code sequence specified in 2.1.3.1.16.

%The I long code for Spreading Rate 1 shall have a chip rate of 1.2288 MHz. The Q long code
%for Spreading Rate 1 shall be the I long code delayed by one chip

%The mask used for generating the I long code for Spreading Rate 1 (or equivalently, the first
%component sequence of the I long code for Spreading Rate 3) varies depending on the
%channel type on which the mobile station is transmitting. See Figure 2.1.3.1.16
UserLongMask = 0x8181;


%Generate 1.2288e6 chips per sedecond, The Wlash spreadng also reulst 1.22288e6 chips
% per second.
ILongCodeUser = [randi([0 1], 1,chipRate)];%   longCodeGnerator(UserLongMask);
% One chip delay, I chip delay mean 0QPSK? MODULATION at final stage.
%The final OQPSK signal(S_{OQPSK}(t)) is the sum of these two components: 
%\(S_{OQPSK}(t)=I(t)cos (omega _{c}t)+ Q (t-T_{b})sin (wc*t)\)(Note:
% \(T_{b}\) is the bit period, which equals \(T_{s}/2\), the symbol half-period delay.) 
QLongCodeUser = [zeros(1,1) ILongCodeUser(1:end-1)];


%ILongCodeUser = 1 -2*ILongCodeUser;
%QLongCodeUser = 1- 2*QLongCodeUser;

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

%chipRate=1.22886e6
IchannelShortCode = randi([0 1], 1,chipRate);  %IchannelShortCodeGenerator();
QchannelShortCode = randi([0 1], 1,chipRate); %QchannelShortCodeGenerator();

%IchannelShortCode = 1-2*IchannelShortCode;
%QchannelShortCode = 1-2*QchannelShortCode;
%Scrambling with LONG CODE, Direct sequence spreading using the long code
%The in-phase spreading sequence shall be formed by a modulo-2 addition of the I-channel
%PN sequence and the I long code sequence.
%The quadrature-phase spreading sequence shall
%be formed by the modulo-2 addition of the following three terms: the
%W1 Walsh function,
%the modulo-2 addition of the I-channel PN sequence and the I long code sequence, and the
%decimated-by-2 output of the modulo-2 addition of the Q-channel PN sequence and the Q
%long code sequence. The decimator shall provide an output that is constant for the two
%chips corresponding to the two symbols of the 2 W1 Walsh function, and the value of the
%decimator output for the 2 W1 Walsh function period shall be equal to the first of the two
%symbols into the decimator in that period. The 2 W1 Walsh function time alignment shall be
%such that the first Walsh chip begins at the first chip of a frame.
%This code has info , multiply two code= modulo 2 addition/xor/spreading as
%ILongCodeComplexD = [mod((IchannelShortCode + ILongCodeUser),2)] ;
%QLongCodeComplexD = [mod((QchannelShortCode + QLongCodeUser),2)] ;

%Scrambling with LONG CODE, Direct sequence spreading using the long code
%The in-phase spreading sequence shall be formed by a modulo-2 addition of the I-channel
%PN sequence and the I long code sequence.

ILongCodeComplexD = [xor(IchannelShortCode , ILongCodeUser)] ;
QLongCodeComplexD = [xor(QchannelShortCode , QLongCodeUser)] ;

QLongCodeComplexD_Dec2   = [QLongCodeComplexD(1:2:end)];

%Walsh spread by 2;
QLongCodeComplexD_Spread = kron(QLongCodeComplexD_Dec2,W3);



ILongCodeComplexD = 1 -2*ILongCodeComplexD;
QLongCodeComplexD = 1 -2*QLongCodeComplexD_Spread;

%This whole will work per frame of encoded audio with 20ms? Call for every frame
% Every frame has 16slot, each slot 1.25 ms,

PilotSymbolPerSlot =18; %bpsk modulated
DataSymbolPerSlot  =6; % 6 bits per slot, along with pilot bits           [pilot data]
DschSymbolPerSlot  =12;  % banwidth channel 12*128 = 1536 walsh spread    [ data     ], Orthigonal

%Increment frame count
frameCount = frameCount + 1;

LongCodeLen = 24576; %1536 per slot, 16,

%18 bits of data bit time multiplex with pilot bits in slot / per slot and 16 slot per frame of 20ms, 50 frames in 1 sec
% every frame data is read and processed.
DataPerFrame20P            =[DataBitsSymbolMappedP((frameCount -1)*DataSymbolPerSlot*16+1:frameCount*DataSymbolPerSlot*16)];


pilotSymPerFrame20ms      =[pilotSymPerSec((frameCount -1)*PilotSymbolPerSlot*16+1:frameCount*PilotSymbolPerSlot*16)];

DschDataPerFrame20       =[dsch2BitsSymbolMapped((frameCount -1)*DschSymbolPerSlot*16+1:frameCount*DschSymbolPerSlot*16)];


ILongCodeComplexD20       =[ILongCodeComplexD((frameCount -1)*LongCodeLen*16+1:frameCount*LongCodeLen*16)];


QLongCodeComplexD20       =[QLongCodeComplexD((frameCount -1)*LongCodeLen*16+1:frameCount*LongCodeLen*16)];



%800 bps, so 1 bit per slot is punctured for control power bit
pbp =1; % 1 symbol is punctured, for on/off of or near or far.

powerSymbol= ones(1,1,pbp);
RepeatationFactor = 64;% needed 1536 bits per slot, 768 after desp, actual is 12,


%This will process slot by slot every 20ms frame encoded audio
% Bits Mapping-->(i/q)symbols-->walsh code spread ->   -->qudrature spreading.
% PilotCanel--->AFTER WALSH CODE SPREAD-->I channel
% traficc    -->After walsh spred------->I AND Q, all all i and q and finally
% IQ, scrambling with comples DS code.


for slotIndex=1:16

    %Crc/Convolution enoding 1/2, interleaving done  is already assumed to be done

    %data bits time multiplex with pilot in reverse link from MS to BS.
    Dsch  = [DataPerFrame20P((slotIndex-1)*DataSymbolPerSlot+1:slotIndex*DataSymbolPerSlot)];
     %Walsh spread
    DschSpreadP = kron(Dsch, W2);

    pilotSymbolPerSlot = [pilotSymPerFrame20ms((slotIndex-1)*PilotSymbolPerSlot+1:slotIndex*PilotSymbolPerSlot)];


    %Osvf spread , channelization , to make 1152 chips per slot, pilotSymbolPerSlot
    %Each symbl is multiplied by Hp, Hadmard walsh code.
    %
    pilotSymbolPerSlotSpread = kron(pilotSymbolPerSlot,W0);

    %Time multiplex   1152+384 =1156 slot with pilot and data bits
    FinalSlotData = [pilotSymbolPerSlotSpread  DschSpreadP];
    powerSymbol = [1 1];

    %Puncturing the data bits for adding power control bits
    %Rc 3, 800 bps that 1 bits
    FinalSlotDataPunctured=[ FinalSlotData(1:end-pbp)  powerSymbol(1:pbp) ];

    %As per data 9.6kbps, with 9.6 symbol per ms, each slot is 1.25ms, so 12 symbols per slot, qpsk modulated
    % In currnet slot , current symbol are spread by W1(2), Orthogonal channel at same time as pilot slot.
    DschSymbol = [DschDataPerFrame20((slotIndex-1)*DschSymbolPerSlot+1:slotIndex*DschSymbolPerSlot)];

    if slotIndex == 1
       disp(DschSymbol);
       figure;
       plot(real(DschSymbol),imag(DschSymbol), 'r.');
       title('Input: 16-QAM Constellation Diagram (Normalized)')
       xlabel('In-phase (I)');
       ylabel('Quadrature (Q)');
       grid on
      axis([-4 4 -4 4]); % Adjust axis limits based on your normalization scale

    end

    % Symbol are spread by W2, factor of 2, 24 symbol after spread but 384-24=360 space is free, so repeatatation of symbol
    % by 30 time, to match the chip rate 1.2288e6 per second.
    DschSymbolPerSlotRepeat = repmat(DschSymbol, 1, RepeatationFactor);

    %Walsh spread channelization code by W1(2);
    %This Walsh function repetition factor is the number of Walsh function sequence repetitions per
    %interleaver output symbol.
    DschSymbolPerSlotSpread = kron(DschSymbolPerSlotRepeat,W3);

    %[pilot data] (+) [ data 2 ] (+) [ data 3   ], all are otrgonally addeded, I and Q stream

    %Scrambing the pilotSymbolPerSlotSprea,1152, complex multiplication, but here its I only
    %pilotSymbolPerSlotScrambled = ILongCodeComplexD(1:1152).*pilotSymbolPerSlotSpread);

    %As per data 9.6kbps, with 9.6 symbol per ms, each slot is 1.25ms, so 12 symbols per slot, qpsk modulated
    % In currnet slot , current symbol are spread by W1(2)
    %DataSymbolPerSlot = [orthogonalDataPerFrame20((slotIndex-1)*DataSymbolPerSlot+1:slotIndex*DataSymbolPerSlot)];

    %Orthogonal addition of I and Q stream

    %Apply gain to each stream of each channel, for now its unity gain
    %IdataGain = Idata *gain; %1.2288e6 = fs, and lpf, has t0 around < 614 khz
    %QdataGain = Qdata *gain; % 16bit I , 16 bit Q, DAC


    IStream = [FinalSlotDataPunctured(1:1152) real(FinalSlotDataPunctured(1:384))]

    IOrthogonalStream = (IStream + real(DschSymbolPerSlotSpread));

    QOrthogonalStreamP = [zeros(1,1152) imag(FinalSlotDataPunctured(1:384))];

    QOrthogonalStream = (imag(DschSymbolPerSlotSpread) + QOrthogonalStreamP);

    ILongCodeComplexDTemp  = [ILongCodeComplexD20((slotIndex-1)*1536+1:slotIndex*1536)];
    QLongCodeComplexDTemp  = [QLongCodeComplexD20((slotIndex-1)*1536+1:slotIndex*1536)];
    %Complex multiplication as quadrature spreading
    IDataSymbolPerSlotScrambled = IOrthogonalStream .*ILongCodeComplexDTemp;
    QDataSymbolPerSlotScrambled = QOrthogonalStream .*QLongCodeComplexDTemp;

     %If there are mutplie SCH OR MORE TRAFFIC channel then , but now its one
    % DataSymbolPerSlotScrambled = DataSymbolPerSlotSpreadChanel2*LongCodeComplexD;

    % Fina data DataSymbolPerSlotScrambled =  DataSymbolPerSlotSpreadChanel1+DataSymbolPerSlotSpreadChanel2
    % so on, all data channel with pilot and power control bits are scambled with same scarmbing code.

    %Fow now lpf/baseband filtering is not impmeneted into this simulation.
    %I and Q are used to modulate the carrier ,OPSK, is done at during qudrature spreading.

    test=[];

    Tx(slotIndex,:) = IDataSymbolPerSlotScrambled +j*QDataSymbolPerSlotScrambled;

end


%In real implemenatation , once we slot indication,
% for every slot indication, we buffer 20ms ping buffer
% while this ping buffer is done, we process this ping buffer in
% some core  and parrelly we can buffer pong..,
% And it depends on ur processing pipline latency to process one frame
% or 1 slot.. ideally we shuld increase the cores for parrell , so that a
% all slots get far befre 20ms .
% each core given , same decsrambled data, based channnel code to be
% assigned to core t process, we can think of such design.

disp("slot index after tx");
disp( Tx(1,:));

ChipRate = 1.2288e6;
rx1 = repmat(Tx,1,50);

rx = reshape(rx1, 1, ChipRate);

%L=128, sample delay
%osf = Ts/Tc;

% sampling rate of symbol data
% Ts symbol, delay has

% be  L < (Ts/2)
%L = 128

%IRxData  = [BasebandIQ(1:1536*2)];

%To avoid intersymbol delay,
Ts = 1/dataRate;
Tc = 1/chipRate;
osf = Ts/Tc;
%1 symbol of data, and 1/2 is half symbol
MaxPathDelay = (1/2)*osf;



BasebandIQ =fadingModel(rx, MaxPathDelay,1);


disp("after fading");
disp( BasebandIQ(1:1536));



%Qrx1 =fadingModel(QDataSymbolPerSlotScrambled,1);




% When we hav channel h(k), , we can have
% 1. we can break the channel svd decompistion, to have procoding
% and channel matrix =  H = UZV, Z, we can have rank of this matrix
% based on rank, we can most L path to be used for rake reciver
% or also causing predistortion at sent so that such channel
% is made flat for all w /Ts,

% This 16 slot 20ms frame is channel

%Since the current slot cntain orthogonal data from all user
%           1152   384 chips          1556
%   slot = [pilot dataUser1] (+) [user2 data ]  (+)
% so same slot contain all the user but orthogonally , This means better orthgonallygnallity
% or long orthogonal code give more possibliy of more users!

% Generate I AND q short code PN seq PN seq of length 15, bases on currnet basestation
%

%bufferIQ = [IQ20msBuffer1  IQ20msBuffer2];
%This is cmplex baseband after , downconversion and LPF
%at 1.2288e6 chip per second.

ChipRate = 1.2288e6;

%Later this has to same as tx size with fading adding.

%BasebandIQ = randi([-15 16], 1, ChipRate) + j*randi([-15 16], 1, ChipRate); %[HardwareDownconversion(bufferIQ);

%BasebandIQ = Irx1 + j*Qrx1; %[HardwareDownconversion(bufferIQ);

% Add channel artificat/fading/awgn and fading
% Rake Receiver at base station

%long code for user leave mask? mobile term id?
%Generate same long code and i and q Pn seq, and finally long code for scrambling the
%

% Later I and Q are scrambled  with short m-sequence, the short code
% diffrenetiate teh signal emiited by diffrenet base station,
% All base station uses same code but with differnet timing offset and
% recever attemts to syncronize to this short code  strongnest corrleation
% idenetifies the nearest base stataion


% Gpss receiver seraches for presence of gpss satellite signal
% by correlting the receved signal with diffrent gold sequence
% since period, is 20ms, same signal seq is repeated in period
% of 1 ms.

% The received input is downconverted to based band using inpahse and quadrature oscillator

% In phase and  Quadrature
% In advance ,
% generate 1second worth of samples
% This generate long code 42 bits and 15 bits , quadrature scrambling code.
% Quadrature
frameCount = 1;

ILongCodeComplexD20       =[ILongCodeComplexD((frameCount -1)*LongCodeLen*16+1:frameCount*LongCodeLen*16)];
QLongCodeComplexD20       =[QLongCodeComplexD((frameCount -1)*LongCodeLen*16+1:frameCount*LongCodeLen*16)];

%pnQuadLongSeqQ = [];%generateQudraturePN(stationId);

% Note sure but need to see, pilot bits are wlash spread and scrambles
% at recver, these final symbols now matched filter
% correleated  for all possible lag at chip unit label
% or just generate decsramled sequence and match filter with input
%
% Pilot bits are known at receiver.
PilotBit  = randi([0 0], 1, 18);

PilotBitS = 1-2*PilotBit;

%scramble the pilot bits
PilotWalshSpread = kron(PilotBitS, W0);

%pnQuadLongSeqQ, is complex qurature descrambler, to
PilotScrambledRef = PilotWalshSpread .*ILongCodeComplexD20(1:1152);



% Read 2 slot worth of symbol, for code syncronization
% Descrambled =
% match filtering for frame sync
% FOR ALL slot in 20ms

%IRxData = [BasebandIQ[1:(2*Ts/Tc))]
IRxData  = [BasebandIQ(1:1536*2)];

Max=0;
offset = 0;
Engery = 0;

slotIndex=1
%for
  IRxDataS1 =  [IRxData];%IRxData((slotIndex-1)*1536+1:slotIndex*1536);
  for chipIndex=1:1151

         %Match filtering corrleation at delay chip index
         Eng = sum(IRxDataS1(chipIndex:(1151+chipIndex) ).*PilotScrambledRef(1:1152));

         Engery = (real(Eng)*real(Eng) + imag(Eng)*imag(Eng));

     if Max < Engery
        Max = Engery;
        offset = chipIndex;
     end
      IRxDataS1 = [zeros(1,chipIndex) IRxData(chipIndex:(1152+chipIndex))];

  end

%end

%Make sure Rxdata buffer is atleast holdng 40ms data, and read every
% 20ms from the buffer at offset ci, and from appliaction, this 40ms buffer is always popluated
% at every 20ms. , read offseted 20ms frame data from  cirecular buffer ? each elemnet
% is buffer of 1.25 ms worth of sample at 1.2288e6 rate.
RxCorrectedIQ = [BasebandIQ(offset:end)];



disp("after correction ");
disp( RxCorrectedIQ(1:1536));

%Need to anayze more on this channel estmation based Pilot

%  hbar = [hbar0 ...];
%  l=0 [ a0  a1  a(L-1     ]  [h0] = [hbar0]
%  l=1                        [h1] =

% estimated pilot = hbar . *yRx ; after sommtehung of h
%
% Here y = [], as recved descrambled  demultiplex pilot, from  l=0 path
% and same form other path l=2


% Fir LPF for smotthening the decsison of channel coeffent of lth path
% let say fir coeff , Need to check ,
% Adapative fir filter, the no of coeffient are based on max delay supported
% a= R
% h = SmoothingChannel(h(l,:);

% LMS UPDATE OF CHANNEL, in Rake reciver . FIR as L tapped delay , time varinat
% channel and its performance.
% It somethinging like passing or recveing the symbols with persion Ts, all consectuive symbols
% passed from such tapped delay lne channel.. as

% symbol rate or Ts, is limited by time variant channel ,
% For all pilot symbols in slot.
% this channel estimation.. in this loop based on rake rever apparoch
%
%18 pilot symbols
M=18;
L=16;

descsrambledI = zeros(L, 1536);
descsrambledQ = zeros(L, 1536);
prvSymbolI    = zeros(L,1);
prvSymbolQ    = zeros(L,1);

hsmothenI     =zeros(L,1);
hsmothenQ     =zeros(L,1);
hi            =zeros(L,1);
hq            =zeros(L,1);
estimatedPilotSymbol = zeros(L,1);
estimatedPilotSymboQ = zeros(L,1);
TempErrorQ = zeros(L,1152);
TempErrorI = zeros(L,1152);


temph = zeros(1,L);
Rhh = randi([1],L,L);
rhx = zeros(1,length(W0));
finalH = zeros(1,M);
rateFactorI =1;
errorI=0;
rateFactorQ=1;
errorQ=0;

%For every pilot symbol, channel estimat is updated. No of Pilot symbols
% are known.
for symbolIndex=1:M % M now can be seen as no of symbols in slot?

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
   %Get curret symbol spread and scrambled
  RxCorrectedI = [real(RxCorrectedIQ(((symbolIndex-1)*1152 +1):symbolIndex*1152))];
  RxCorrectedQ = [imag(RxCorrectedIQ(((symbolIndex-1)*1152 +1):symbolIndex*1152))];

  %
  for l=1:L-1

      % Descrambling wrt station code,
      % hk = conjugate(bk) . *ypk  =  *bk .*(ak*bk + wk) = *bk.*bk*ak + *bk.*wk where |bk*bk| = 1
      %  = ak + noisek, ak is channel coeffient?? which is cmplex??
      % Genrate pilot symbol, walsh sperad and then scmable with short code
      % and then same
      % first cut hk, is found to be for any path.., generally pilot
      % signal are speard with only short code .

      start = l;
      end1   = 1152 + start;
      % all delay in colmwise , l col by adding zero, as shift
      RxShiftedI = [real(RxCorrectedIQ(start:(end1-1)))];
      RxShiftedQ = [imag(RxCorrectedIQ(start:(end1-1)))];

      %descsrambledI = RxShiftedI.*real(pnQuadLongSeqQ); %Tc/Ts=N
      %descsrambledQ = RxShiftedI.*imag(pnQuadLongSeqQ);
      tempHI(l) = real(PilotScrambledRef)*RxShiftedI(:);
      tempHQ(l) = imag(PilotScrambledRef)*RxShiftedQ(:);

      temph(l) =  tempHI(l) +j*tempHQ(l);

    end

    % [1]           [1][2][3] =
    % [2]
    % [3]

    rhx = temph.*[pilotSymPerSec(1:L)];
    Rhh =  temph(:)*temph;
    Rhh = Rhh + randi([1 1],L,L);

    %LXL : L
    %Compute final channel estimate for this current symbol
    weight  = (1./Rhh)*rhx(:);
    finalH(symbolIndex) = (1/L)*temph*weight

end

  %interpolate finalH

  interplatedH = [ finalH  repmat(finalH(end),1,24-L)];



  %

  %At given instant of slot, 1152 channel estimate , we will interpolate

  %to to have 1536?

% For every user channel estimate is called as above but here once
% channel estimate for particular user is done then For all channel
equalizedSymbol=zeros(1,24);
lengthWalsh =length(W2);
SymbolPerSlot=1536/lengthWalsh;
equalizedSymbol=zeros(1,SymbolPerSlot);
%In slot 1536/W2=12 SYMBL
offset = 1536-lengthWalsh;

for symbolIndex=1:SymbolPerSlot

  %Get curret symbol spread and scrambled
  RxCorrectedI = [real(RxCorrectedIQ(((symbolIndex-1)*lengthWalsh +1):symbolIndex*lengthWalsh))];
  RxCorrectedQ = [imag(RxCorrectedIQ(((symbolIndex-1)*lengthWalsh +1):symbolIndex*lengthWalsh))];

  pnQuadLongSeqI = [ILongCodeComplexD20(((symbolIndex-1)*1536 +1):symbolIndex*1536)];
  pnQuadLongSeqQ = [QLongCodeComplexD20(((symbolIndex-1)*1536 +1):symbolIndex*1536)];

  %descsrambledI = RxShiftedI.*pnQuadLongSeqI; %Tc/Ts=N
  %descsrambledQ = RxShiftedQ.*pnQuadLongSeqQ;

  equalizedSym=zeros(1,L);
  for l=1:L-1

      %start = l;
      %end1   = lengthWalsh + start;
      % all delay in colmwise , l col by adding zero, as shift
      %RxShiftedI = [zeros(1,l)  real(RxCorrectedI(start:(end1-1)))];
      %RxShiftedQ = [imag(RxCorrectedQ(start:(end1-1)))];

      RxShiftedI = [zeros(1,l)  real(RxCorrectedI) zeros(1,(offset-l))];
      RxShiftedQ = [zeros(1,l)  real(RxCorrectedQ) zeros(1,(offset-l))];

      descsrambledI = RxShiftedI.*pnQuadLongSeqI; %Tc/Ts=N
      descsrambledQ = RxShiftedQ.*pnQuadLongSeqQ;

      %W1 is channelization code, so every sybol has orthognal symbols addeded
      % with diffrent spread orthogonal code.
      rxI = reshape(descsrambledI,length(W3),[]);
      despreadI = (W3(:)' *rxI)/length(W3);


      %equalizedI = finalH(symbolIndex)*despreadI
      rxQ = reshape(descsrambledQ,length(W3),[]);
      despreadQ = (W3(:)' *rxQ)/length(W3);

      %Rest are punctured or inored. as IT HAS BEEN REPEATED TO RATE MACTHED IN ENCODER
      rx = despreadI(1:18) +j*despreadQ(1:18);
      equalizedSym(l) = sum([finalH(1:length(rx))].*rx);
   end

   equalizedSymbol(symbolIndex) = sum(equalizedSym);

 end


   disp("slot 1 after decoding ");
   disp( equalizedSymbol);

   figure;
   plot(abs(equalizedSymbol), 'b');
   hold on;
   plot(abs(DschSymbol),'g');
   title('Output:16-QAM Constellation Diagram (Normalized)');
   hold off;

%for channelIndex=1:MaxChanell

    %Dsch
%    for





    %Summ all path, assuming channel estimate is flat for 20ms, time invraint
%//    equalizedSymb(slotIndex) = sum(equalizedSymbol(l);

    %dechannelization
    %dechanalize for ather
 %//   equalizedSymbPich = sum(equalizedSymb(slotIndex).*W0);
 %//   equalizedSymbCh1 = sum(equalizedSymb(slotIndex).*W1);
 %//   equalizedSymbCh2 = sum(equalizedSymb(slotIndex).*W2);

%//end





