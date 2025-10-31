
%How , power is spread across, M path , it guven list of
% of offset, phi, with equal power
%Rms agular spread,
function [phi] = PasSubArry(sigma)

  %sigma = 2;
  % AS is 2 deg, differnt offset from table, but also you can
  % calcluate based on formaule ,
  if sigma==2
    % M= -10  to 10, will have eqalu power
    % 0 to M-1/2;
    phi = [0.894 .2826 .4984 .7431 1.0257 1.3594  1.17688 2.2961 3.0389  4.3101 ];
    disp(phi);

  elseif sigma==5
   phi = [0.894 .2826 .4984 .7431 1.0257 1.3594  1.17688 2.2961 3.0389  4.3101 ];


  elseif sigma==35
   phi = [0.894 .2826 .4984 .7431 1.0257 1.3594  1.17688 2.2961 3.0389  4.3101 ];

  else
   disp('errr');

  end


end







