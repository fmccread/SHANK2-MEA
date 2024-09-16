%% ENTROPY ALGORITHM 
% This algorithm calculates the correlated spectral entropy (CorSE) based functional connectivity between
% two time series. For more information, <http://journal.frontiersin.org/article/10.3389/fncom.2016.00112/full Kapucu et al.(2016)>. Please cite the paper if you have used or inspired by the method.

%% Introduction
% For this version, algorithm first calculates time variant spectral entropies (SE)
% from all the channels. Multichannel data, e.g., MEA data, has formatted
% as DataCell.

%% Parameters
fs = 12500; %sampling freq.
WS = fs*1; %window size to calculate Spectral Entropy in samples fs for 1s
overlap=0.5;%other than 0.5 makes it more complex and code needs to be changed
TotalDataSample = length(DataCell{1,7});%length data in samples
SElength = round(2*(TotalDataSample / WS) - 1);
chNO = height(DataCell);%number of channels to process
    
% 60 Hz notch filter to remove power line noise
%d = designfilt('bandstopiir','FilterOrder',2, ...
%                   'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
%                   'DesignMethod','butter','SampleRate',fs);

[b,a]=rico(0.01,6,60,1/fs); 
               
A=zeros(chNO, SElength-1);


    parfor c = 1:chNO
        
           x = DataCell{c,7};
    
           s1 = CorSE_body(x,fs,WS,overlap,SElength,b,a);
    
           A(c,:)=s1;
    
    end
    

    