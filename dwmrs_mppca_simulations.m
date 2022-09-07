%% Simulation of diffusion-weighted MRS data for the evaluation of MP-PCA 
% Author: Jessie Mosso 
% Last update: 08/2022
% Reference: https://doi.org/10.48550/arXiv.2204.04266
% Copyright © 2022 jessie-mosso

%%
clear; clc; close all;

%% Simulation parameters to choose: 
decaytype="stickmodel"; %diffusion model: 
%"stickmodel": S/S0=sqrt(pi/(4bD))erf(sqrt(bD))(Ref: Callaghan et al. 
%Biophys J 28, 1979)
%"monoexp": S/S0=exp(-bD)
phasedist=false; %true/false: adds a random phase to each spectrum from a 
%[0,30]deg uniform distribution
b0drift=false; %true/false: adds a random frequency drift to each spectrum 
%from a [-15,-15]Hz uniform distribution
addwater=false; %true/false: adds a residual water signal at 4.7ppm, 16Hz 
%Lor broadening,random phase in a Gaussian [-pi,pi] distribution, 
%mono exp. diffusion decay with D=0.2um^2/ms
singleavSNR=13; %time-domain SNR on a single shot
MM=true; %true/false: adds the macromolecules to the basis-set
na=100; %number of shots per b-value
np=2048; %number of complex points in the FID (Note: the basis-set has 
%4096 points, cut here to 'np')
nMC=1; %number of Monte-Carlo generations of the diffusion matrix
pathfidraw=[pwd '\basis-set\']; %basis-set folder location 
dt=200e-6; %dwell time of the basis-set, sec
sw=5000; %spectral width of the basis-set, Hz
B0field=400; %B0 field of the basis-set, MHz
doplot=true; 
doMPPCA=true;
pathdncodes=[pwd '\mp-pca_VeraartNeuroimage2016\'];%MP code folder location 


%% Load metabolites FID (.RAW)
metnames=["ala","asc","asp","cr","pcr","gaba","gln","glu","gsh","ins",...
    "lac","naa","scyllo","tau","glc","naag","pe","gpc","pcho"]; 
if MM
    metnames=[metnames,"MM_full_fit_new"]; 
end
nbptsrawfiles=4096; %number of points in the basis-set
metfids=zeros(nbptsrawfiles*2,length(metnames));

for nmet=1:length(metnames)
    fid = fopen([pathfidraw + metnames(nmet) + ".RAW"]);
    C = textscan(fid, '%f', 'headerlines', 2);
    metfids (:,nmet)= C{1,1};
    fclose(fid);
end

%% Simulation parameters
% metnames=["ala","asc","asp","cr","pcr","gaba","gln","glu","gsh","ins",...
%     "lac","naa","scyllo","tau","glc","naag","pe","gpc","pcho"]; 
conc=[0.8,1.5,2,4,4.5,1.6,3.0,10,1.5,6.5,0.8,9,0.1,4.5,1.7,0.3,0.5,0.8,0.2];
scalefactor=0.07; 
if MM
    conc=[conc,1.3];
end
conc=conc.*scalefactor;

%b-values
bval=[0.4,1.5,3.4,6, 7.6,13.4, 15.7, 20.8, 25.2,33.3];
%diffusion coefficients
Ddiff=5*[0.0539,0.0623,0.134,0.1,0.1,0.0756,0.0768,0.1,0.0531,0.09,0.13,...
    0.08,0.0761,0.11,0.114,0.08,0.0636,0.09,0.09,0.001];
%Ddiff tCr: 0.1
%Ddiff tCho: 0.09
%Ddiff set to 5 times the values from Ligneul et al., Neuroimage 2018 
%- when reported, otherwise choosen randomly between 0.05 and 0.15 um^2/ms

%% Construction of the DWMRS matrix
Mnoiseless_t=zeros(length(bval),np);
for MCit=1:nMC % for Monte-Carlo studies: generates the matrix nMC times
    disp(["MCit=" + num2str(MCit)])
    %Diffusion matrix: will contain complex FIDs
    Mdiff_t=zeros(na*length(bval),np);  
           
    %% Random fluctuations definition
    randomphases=rand(1,length(bval)*na).*deg2rad(30); 
    %random phase from a [0,30]deg uniform distribution
    randomshifts=-0.02+rand(1,length(bval)*na).*0.04; 
    %random frequency drift from a [-15,-15]Hz uniform distribution

    for b=1:length(bval)
        Mav=zeros(na,np);  %submatrix for 1 bval, na shots
        
        %% Construct the 1D spectrum for each b-value
        temp0=0;
        for c=1:size(metfids,2)
            
            metbs=[metfids(1,c)/2; metfids(2,c)/2; metfids(3:np*2,c)]; 
            %divide the first real and imag point by 2 to remove 
            %the frequency offset caused by fft in Matlab
            %cut to half 2*2048 (tot: 2*4096)
            %NB: in the .RAW files, real and imag are interleaved
            %special processing for the MM
            if c==20 %flip the MM spectrum
                metbs_r=metbs(1:2:end);
                metbs_i=metbs(2:2:end);
                flipmetbs_r=metbs_r;
                flipmetbs_i=-metbs_i;
                metbs_new=zeros(np*2,1);
                metbs_new(1:2:end)=flipmetbs_r;
                metbs_new(2:2:end)=flipmetbs_i;
                metbs=[metbs_new(1)/2;metbs_new(2)/2;metbs_new(3:end)];
            end
            
            %decay model
            if decaytype=="stickmodel"
                temp0=temp0+metbs*conc(c)*sqrt(pi./(4.*bval(b).*Ddiff(c))).*erf(sqrt(bval(b).*Ddiff(c))); 
            elseif decaytype=="monoexp"
                temp0=temp0+metbs*conc(c)*exp(-bval(b)*Ddiff(c)); 
            end
        end

        temp=temp0';
        rr=1:2:(length(temp)-1); %real
        ii=2:2:length(temp); %imag
        
        spectrum_c_t=temp(rr)+1i*temp(ii); %complex FID at 1 bval
        
        %% add water residual
        tt=dt:dt:np*dt;
        ph_water=-pi+randn(1)*2*pi;
        fid_water=exp(-tt.*16*pi).*exp(i*ph_water).*exp(-bval(b)*0.2);
        waterscaling=11.4;
        if addwater
            spectrum_c_t=spectrum_c_t+waterscaling.*fid_water;
        end
        
        %% additional 5Hz Lorentzian LB applied to all spectra
        spectrum_c_t=spectrum_c_t.*exp(-tt*5*3.14);
        

        %% Define noise level from the first shell, first MC iteration 
        if b==1 && MCit==1 
            sigmanoise=abs(spectrum_c_t(1))/singleavSNR; 
        end
        
        %% construct the na shots per b-value with different noise 
        %% generations and distortions
        for n=1:na
            spectrum_c_tnoisy=spectrum_c_t+sigmanoise*(randn(1,np)+1i*randn(1,np)); 

            %% apply B0 drift distortions
            if b0drift
                shift_spec=randomshifts((b-1)*na+n);
                spectrum_c_tnoisy=spectrum_c_tnoisy.*exp(i.*shift_spec.*[0:1:np-1]);
            else
                spectrum_c_tnoisy=spectrum_c_tnoisy;
            end
                
           %% add phase0 distortions
            if phasedist
                ph0=randomphases((b-1)*na+n);
                spectrum_c_tnoisy=spectrum_c_tnoisy.*exp(1i*ph0);
            else
                spectrum_c_tnoisy=spectrum_c_tnoisy;
            end
            Mav(n,:)=spectrum_c_tnoisy;
        end
        
        %% Store noiseless matrix
        if MCit==1 
            Mnoiseless_t(b,:)=spectrum_c_t;
        end

        Mdiff_t((b-1)*na+1:b*na,:)=Mav; %store the shell in the diffusion matrix
    
    end
    Mdiff_t_mc(:,:,MCit)=Mdiff_t; %For Monte-Carlo studies: store the 
    %matrix at each MCit iteration
end 

if doplot
    Mplot=Mdiff_t_mc(:,:,1);
    specmat=fftshift(fft(Mplot,[],2),2);
    figure; 
    fmax=sw/2;
    f=[fmax:-2*fmax/(np-1):-fmax];
    ppmscale=f/B0field+4.7;

    for b=1:length(bval) %plot all b-values overlapped
        plot(ppmscale,real(sum(specmat((b-1)*na+1:b*na,:))))
        hold on 
    end 
    set(gca,'xdir','reverse')
    xlabel('ppm')
    title('Original')
end

if doMPPCA
    addpath(pathdncodes)

    for MCit=1:nMC
        disp(["MPPCA on MCit=" + num2str(MCit)])
        centering=true; %do centering before MP-PCA
        %denoising on temporal-domain matrix, concatenation of Re and Im
        Mfordn_t=[real(Mdiff_t_mc(:,:,MCit));imag(Mdiff_t_mc(:,:,MCit))];
        [dn_Mdiff_t_concat,sigma,npars]=MP(Mfordn_t,50,centering);

        %denoised temporal matrix
        for n=1:na*length(bval)
             dn_Mdiff_t(n,:)=dn_Mdiff_t_concat(n,:)+1i*dn_Mdiff_t_concat(n+na*length(bval),:);
        end
        
        dn_Mdiff_t_mc(:,:,MCit)=dn_Mdiff_t;  %For Monte-Carlo studies: 
        %store the matrix at each MCit iteration
    end 
    if doplot
        Mplot=dn_Mdiff_t(:,:,1);
        specmat=fftshift(fft(Mplot,[],2),2);
        figure; 
        fmax=sw/2;
        f=[fmax:-2*fmax/(np-1):-fmax];
        ppmscale=f/B0field+4.7;

        for b=1:length(bval) %plot all b-values overlapped
            plot(ppmscale,real(sum(specmat((b-1)*na+1:b*na,:))))
            hold on 
        end 
        set(gca,'xdir','reverse')
        xlabel('ppm')
        title('MP-PCA')
    end
end 

    


