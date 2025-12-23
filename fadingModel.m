
function [yCorrP] = fadingModel(data,L,type)

    %pdp, power delay profile, help to understand how is mean delay spread is
    % analyized using channel impulse response, and how once u understand
    % you pdf, it help t understand time disprevise fading or freq disvrive
    % fading and what has symbol perod wrt time dispersive figure
    % signal bw , based on freq disverpisive..
    %td = % as max delay
    %Bc = 1/tc;

    % so symbol period has to be greater td, delay spread ( mean delay, which is
    % sutied as rndom experiment?
    % all the delayed compmenet signal , will be accumulated before next symbol,
    % causing no interface...
    %ts  > td

    %Bs, where each signal tne, is flat faded, same gain, unity gain
    %Bs < Bc;

    % y(t) = h(t)*x(t),
    % fc = c/y
    % fv = (v/lamda)*cos(theta), theta angle of arrival

    % 0 <theta < pi;, so when whe theta cos theta = 1, , -1
    % max dopller shift is
    % fm  =  max dppler - min dopler
    %

    % Religh distribution, is h(t) as complex  with varying phase and ampliture
    % awgn here h(t) as awgn model, flat, contsat ampliture or fix phase whch
    % has to be corrceted at synrimzation time, and noise level is lower
    % religh as chnnel model , here constant phase tracjing is needed, good
    % channel estimation and equlization is needed.
    % L=128, sample delay
    % osf = Ts/Tc;
    % sampling rate of symbol data
    % Ts symbol, delay has to be  L < (Ts/2)
    %L = 128
    %c = 1; %? % amplitude of LOS, c + W1 + jW2
    %Rician Model:
    %sigma = .9; %? % variance of scattered componenet, if its is random
              % then what is varaince? or sigma sqr.
    %x for every L path, there is one complex random variable
    %x = randi([0 1], 1, L) + j*randi([0 1], 1, L);

    %m=1:10;
    %figure;
    %plot(m, x, '0');
    %grid on;

    %figure;


    %plot(x);

    % Asume   0 < x < 1, so x is very small
    %I0_x = 1 + (x*x/4); % higher valye of x*x*x*,, can be ignored.
    %argxc = (x*c)/sigma^2;

    %I0_x = 1 + (argxc.^2)/4; % higher valye of x*x*x*,, can be ignored.
    %argxc , aume to be < 1;, c is unity, 0.9*.9/(.18) < .1

    %argsExp = -1*(x.^2 +c^2)/(2*(sigma^2));

    %expA = exp(argsExp);


    %h = (x/(sigma^2)).*(expA.*I0_x);
    % complete avove alogorth is siuated as randn(); which is done in reliegh file.
    % compute complex h
    level = ceil(1 + log2(L))
    %level=30;
    h = Rayleigh(L);

    k_db = [-40 -100];

    for k=1:length(k_db)
       hr = Rician(k_db(k),h);
       %mag = hist(abs(hr(1,:)), level);
       %figure;
       %plot(abs(hr(1,:)))
       %plot(mag);
    end


    %y = conv(hr,data);
    yConv = [zeros(1, length(data))];

    ProcessLen = length(data) - L ;
    %Convolution.
    for m=1:ProcessLen-1

        window=[data(m:(m+L-1))] ;
        %inner product
        sum = (window *hr(:))/length(hr);
        yConv(m) =sum;
    end

    yCorrP = yConv;

    %histogram shw hw many times the same value is there , value h and its freq
    % level is bins, seems histogram , now compare this model or histogram
    % with actual reigh distritbution..
    %mag = hist(abs(hr(1,:)), level);

    figure;

    plot(abs(yCorrP(1,:)))


    xlabel('x as random value of (iid) as possible CIR impulse response');
    ylabel('h as h from L PATH');
    title('ricain model');


end

















