
function [H]= Rayleigh(L)
  %rand normal distribution, H is magnitute response which is compplex.
  mean = 0;
  variance = 10; % Avg Noise power

  %sigma sqare is noise power

  %Random probablity, what has to be complex ..value.
  pdf1 =randi([1 10],1,L)/10;


  temp1 = sqrt(2*3.14*variance);
  %dispersion = sqr((x-variance))/(2*variance);
  %pdf = temp1* e(-dispersion)
  %log(pdf) = log(temp1) -dispersion;

  invPdf = 1./pdf1;

  x1 = sqrt(log(temp1*invPdf)*2*variance) + variance;

  pdf2 =randi([1 10],1,L)/10;

  invPdf2 = 1./pdf2;
  x2 = sqrt(log(temp1*invPdf2)*2*variance) + variance;

  H = x1+j*x2;
end



