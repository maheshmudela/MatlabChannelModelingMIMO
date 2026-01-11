
%Fir filter using app, Optimal appraoch, for designing filter
% needed for various use in communication to have impact least on phase, and
% mantin phase linearity and contiuity.
% Instead of window appraoch , park mclellan optimal appraoch given same
% specifucation with less filter computataion, make filter operation
% fast.
% This appraoch will give best and lleast coeffient.
% Instead langrange interpolation, need to check other
% Objective is to for given specification for any filter design , how many
% Coeffient are needed to make accurate sharp filter
% For 100Mhz band , ?? Per CC componenet carrier which has max 100mhz
% if want 400mhz, 4 cc will be used to support.
%
%Ideal freq response, log mag of H(w)
% H(w0) = 20log10(1+delp)
% H(w1) = 20log10(1)=0;  1st mimima
% H(w2) = 0
% H(W3) = 0
% H(wp) = 20log(1-delp);
% H(ws) = 20log(dels);
% H(w4) = 0;
% H(w5) = 0;
% H(w6) = 0;
% H(w=pi) = 20log(-dels);
%
%     [ L  wp  ws  delp  dels  fiterCoeff[1:L]]
%M=
%

close all;
clear all;

% z=w1+jw2, complex freq eigen function,
% Creating filter:
%del1=     20log10(0.01)= -2*20=-40db
%del2 =    20log10(0.001)= -60db
samplingRate=20e6; % 2pi, fs/2 <=>pi,
%wp =   % fs/2<=>pi ==   w=2pi*f  , fs , w= 2*pi*fs/2=pi*fs , say f= fs/3
%SIGH  is unit pi/T,
 %2*pi*f
N=13;
wp=0.4*pi;  % f= 2000, fs , wp = 2*pi*(f/fs), normalizd.  fs=f===>2pi.

%  |H(jw)|  < (1+-.01), IN Linear domain..  w<wp.
ws=0.5*pi;

%----Wp
%--------Wo(cut off)
%----------Ws;
N =8
%Filter specification
w=[];
Ww=[];
Dw=[];
w0= (wp + ws)/2;

if (mod(N,2)==0)
  L=N/2 + 2;
else
L = ((N-1)/2) +2; % 2 , 6 roots , where polynomila show sign changes
end
% 0----------------------- wp----ws-------------------pi
%-z------------- -z----------
%-w0--z---w1--z---w2--z----wp----point of mxima and minim, sign altenation at wk
%---------z----------------z-
%                             z

%                                 z        z
%---------------------------------ws--w5---w6----pi--------------------------------
%                                     z          z

%Consider these extrmal frew, with equriripple say
% L roots, resolution for ripple width= pi/6
%
delp=0.01
dels=0.0001;
K=delp/dels;
Wp=1/K;
Ws=1;

%L=  32;
w = [1:L+2]'*(pi/(L+2));

W = Wp*(w<=wp) + Ws*(w>=ws); %Weight
D = 1*(w<=w0);
shift=0;
k=1:L+2
resPi = pi/(L+2);
nextFreqSet = resPi/L;
wk = [resPi .*k'];
xk = [cos(wk)];
xi = [cos(wk)];
Max=1
b = randi([0 1], 1, L+2);
% K=delp/dels
% W = 1/K for  0<w < wp
% W = 1  for   w>ws to pi
%Iterative approach to compute freq response based on optimal and lagrange
%inuterpolation as used by Mcclellan INSTEAD of chebyshev approximation.s
% We can try?
for iteration=0:Max

    %using trignomtericl polynomail interpolation instead of solvin linear eqtion
    % to comput   P(w) = a(k)*cos(kw) k=0 to L+1  k=0 where w=0 and k=L+!, w=pi
    % E(w) = W(w)(H(w)-P(w)), max error at any point in freq resonse
    % Max E(w) = delp  , dels in stop band in linear range delp=0.01, dels=0.001
    % at any point E(w)
    %     L+2
    %bk = PROD( 1/(xk-xi)
    %     i
    % x = cos(wk)
    % compue P(w)  ??
    %Compute Mdel

    % E(W) = W(w)(H(w) - P(w), if for these set of all freq
    % E(w) <= Mdel, then optimal solution is found..
    for i=1:L+2
      p=1;
      for j=1:L+2
         if i ~=j
           %xk, are cosne basis vector.
           b1 = 1/(xk(i)-xi(j))
           %figure;
           %plot(wk, b);

           p = p*b1;
         end
       end
       b(i)=p;
    end
    % K=delp/dels
    % W = 1/K for  0<w < wp
    % W = 1  for   w>ws to pi

    %Desired, Maxima points where sign vary
    %  H(w=0) = 20*log10(1+delp)
    %  H(w=1) = 20*log10(1-delp)
    %  H(W=2) = 20*log10(1+delp)
    MaxDelp = 20*log10(1+delp);
    MinDelp = 20*log10(1-delp);

    MaxDels = 20*log10(dels);
    MinDels = 20*log10(-dels);



    %desired
    %repmat(A, n, m): Replicates the entire array \(A\) to create an \(n\times m\) grid of that array.
    Hp= repmat([MaxDelp MinDelp],1, 2);
    Hs= repmat([MaxDels MinDels],1, 2);
    H=[Hp Hs];
   % plot(abs(power(10,H*.05)))

    % plot([1:length(H)],H)
    bk=[b];
    %Matrix Inner Product	num = bk' * H;
    num    = sum(H .*bk) %1 to L+2; maxima +
    weight = b./ W';
    %sign = [-1 1 -1 1]
    sign = repmat([1 -1],1, 4);
    denm = weight .*sign;
    denom1 = sum(weight .*sign)

    Mdel = num / denom1;
    %Logic: It scales each weight bk by the distance of its corresponding extremal
    % frequency from the final extremal point xk(end).
    %This is a standard normalization step in the barycentric interpolation
    %formula to ensure numerical stability.
    %Fix: Ensure both are the same
    %orientation: dk = bk(:) .* (xk(:)' - xk(end)); (this forces both to be compatible).
    dk= bk .* [xk' - xk(end)];

    ck = H - denm*Mdel;

    %finally estimating t=he p(w)...Langrange interpolation at any x.
    % w =  as any point in  0<w <pi
    %Use the linspace function to generate exactly 8 points, starting at 0 and ending at \(\pi \)
    %w = rand([0 pi], 1, L+2);
    % 0 : pi/7 : pi
    w = linspace(0, pi, 8);

    %dk1 = dk./w-xk
    %num = ck.*dk1;
    %estimated
    %pw = sum(num)/sum(dk1);

    %these lines perform the barycentric evaluation of
    %the frequency response across a dense grid to find new extremal frequencies
    %dk1 = dk./w-xk
    %Ensure w is a column and xk is a row to create a matrix of distances
    dist = w - xk';
    dk1  = dk(:)' ./ dist;

    %num = ck.*dk1;
    %dk1 = dk(:)' ./ dist;
    %pw = (dk1 * ck(:)) ./ sum(dk1, 2);
    pw = (dk1 .*ck) ./ sum(dk1, 2);

    figure;
    plot(abs(power(10,H*.05)))
    figure;
    plot(abs(power(10,pw*.05)))

    %estimated
    %pw = sum(num)/sum(dk1);
    %if status < L+2, it means not staified
     %for all w P(w) < Mdel
    %Status = sum(sign(pw(:)-Mdel) )

    %Move to next set of freq where maxima could be found.
    %if (Status)
    %  wnk = nexFreqSet + wk';
    %  wk = [wnk]
    %end

    Desired=[H];
    Weighting=[weight];

    % Calculate the error function across the dense grid
    Error = Weighting .* (Desired - pw);

    % Find indices where the error function has local maxima/minima
    % This replaces the 'Status' logic to ensure we find L+2 points
    [~, locs] = findpeaks(abs(Error));

    % Check if we have found enough extremal points
    Status = length(locs);

    if Status >= (L + 2)
        % Update the extremal frequencies with the new peak locations
        % Sort them to maintain the required ascending order
        wk = sort(w_grid(locs(1:L+2)));
    else
        % If not enough peaks found, adjust the grid or step size
        wnk = nextFreqSet + wk';
        wk = sort(wnk);
    end

     xk = [cos(wk)];
     xi = [cos(wk)];

end

%Check the resonse for given L, other wise increase L ,
% Make sure what is best L

