# FCPCA
Combined analysis of Amplitude and Phase variations in Functional Data

Matlab software for 
  “Combined Analysis of Amplitude and Phase Variations in Functional Data”
  by Sungwon Lee and Sungkyu Jung

Preparation: 
  Requires the installation of Matlab package “SRVF_FDA” by 
  Statistical Shape Analysis & Modeling Group at Florida State University. 
  For this, download .zip file at http://ssamg.stat.fsu.edu/downloads/SRVF_FDA.zip,
  unzip and add paths to Matlab. 

Brief descriptions of main routines:
  FCPCA.m: Performs the Functional Combined Principal Component Analysis
  FCCCA.m: Performs the Functional Combined Canonical Correlation Analysis
  FCPCAvis.m: Visualize the major modes of variations estimated by FCPCA

Example codes:
%%%%%%%%%%%%%%
clear;
load berkeley_male_velocity; % discretized function = f, time points = t 
pc = FCPCA(f,t);

for jj = 1:4;
figure(100+jj);
FCPCAvis(pc,jj)
print('-dpng',['berkeley_FCPCA' num2str(jj)])
print('-dpsc',['berkeley_FCPCA' num2str(jj)]) 
end

cc = FCCCA(pc,15,10);  
%%%%%%%%%%%%%%%
