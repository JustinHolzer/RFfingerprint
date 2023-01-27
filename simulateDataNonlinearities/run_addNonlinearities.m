clear all
close all
num_Txrs = 15;
flag_debug = 0;
samples_in_frame = 160;
num_frames_per_txr = 10; % =10 if NN is already trained; = 5000 if used to train NN
training_frac = 0.8;
text_frac = 1 - training_frac;


%X_training = zeros(samples_in_frame,num_frames_per_txr*training_frac*numTxrs)


%% input parameters for nonlinearities: phaseOffset, freqOffset, phNzFreqOff, ampGain, sampRate
% input parameter:   Range assumed reasonable (units)
% phaseOffset:       0-90; (degrees)
% freqOffset:        TBD maybe 1-100; must be>0(Hz) Note: freqOffset smears constellation around unit circle (constant modulus signals)
% phNzFreqOff(1)     TBD (maybe -35 to -60 <0 phNzFreqOff(1) [dBc/Hz] Phase noise level,  NOTE: more negative = less error introduced
% phNzFreqOff(2)     1-1000 [Hz] Frequency offset 
% Ia;                0-TBD 3; %[dB] Amplitude imbalance
% Ip;                0-TBD 3; %[degrees] Phase imbalance
% IDC;               0-TBD 0.01; %[?] DC offset applied to I
% QDC;               0-TBD -0.02; %[?] DC offset applied to IQ
% ampGain:           0-50 [dB]  0=do nothing, 50= noticable change in scatter plot


%% (digital) Transmit Signal Definition
% note: will replace this with 5G signal
M = 4;
data = randi([0 (M - 1)],samples_in_frame*num_frames_per_txr,1);
pindBm = 10; % (dBm) signal power level
pin = 10.^((pindBm-30)/10); % Convert dBm to linear Watts
x = qammod(data,M,'UnitAveragePower',true)*sqrt(pin); %modulate data and set input power level (pin)
if(flag_debug)
    scatterplot(x); title('Input Signal')
end

%% define limits to parameters 
phaseOffset_limit = [0,30]; %[degrees]
freqOffset_limit = [1,100]; %50; %[Hz]
phNzFreqOff1_limit = [-1,-50];%-50;
phNzFreqOff2_limit = [1,1e3];
Ia_limit = [1, 10];%3; %[dB] Amplitude imbalance
Ip_limit = [1, 10];%3; %[degrees] Phase imbalance
IDC_limit = [0.0001, 0.1];%0.01; %[?] DC offset applied to I
QDC_limit = [-0.0001, -0.1];%-0.02; %[?] DC offset applied to IQ
ampGain_vec = [0, 50]; %20;

phaseOffset = [];
freqOffset  = []; 
phNzFreqOff1 = []; 
phNzFreqOff2 = []; 
Ia          = [];
Ip          = [];
IDC         = [];
QDC         = [];
ampGain     = [];
X_nl        = {};
%%  randomize units within limits set above
for i = 1:num_Txrs
    %% find parameters for nonlinearities from rand number generator
    %  x_out = randRange([xmin, xmax])
    phaseOffset = [phaseOffset randRange(phaseOffset_limit)];
    freqOffset = [freqOffset randRange(freqOffset_limit)];
    phNzFreqOff1 = [phNzFreqOff1 randRange(phNzFreqOff1_limit)];
    phNzFreqOff2 = [phNzFreqOff2 randRange(phNzFreqOff2_limit)];
    Ia = [Ia randRange(Ia_limit)];
    Ip = [Ip randRange(Ip_limit)];
    IDC = [IDC randRange(IDC_limit)];
    QDC = [QDC randRange(QDC_limit)];
    ampGain = [ampGain randRange(ampGain_vec)];
    
    IQ_imbalance = [Ia(i), Ip(i), IDC(i), QDC(i)];
    phNzFreqOff_vec = [phNzFreqOff1(i) phNzFreqOff2(i)];
    
    %% add nonlinearities to signal
    flag_debug=1;
    [x_nl] =addNonlinearities(x, phaseOffset(i), freqOffset(i), phNzFreqOff_vec, IQ_imbalance, ampGain(i), flag_debug);
    
    %% add awgn ??? (TBD)
    
    %% save off each signal appropriately
    X_nl{i} = x_nl;
    x=x_nl;
    filename = ['tx',num2str(i),'.mat']
    tx_phaseOffset = phaseOffset(i);
    tx_freqOffset = freqOffset(i); 
    tx_phNzFreqOff_vec = phNzFreqOff_vec; 
    tx_IQ_imbalance = IQ_imbalance; 
    tx_ampGain = ampGain(i);
    
    save(filename,'x','tx_phaseOffset', 'tx_freqOffset','tx_phNzFreqOff_vec', 'tx_IQ_imbalance','tx_ampGain')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_out] =addNonlinearities(x, phaseOffset, freqOffset, phNzFreqOff, IQ_balance, ampGain, flag_plot)
%% inputs: x - input signal to modify e.g. QAM modulated signal
%%         phaseOffset [degrees], freqOffset[Hz], 
%%         phNzFreqOff(1) [dBc/Hz] Phase noise level,  phNzFreqOff(2) [Hz] Frequency offset 
%%         IQ_balance 
%%         ampGain [dB] amplifier gain (dB)
%%         flag_plot
%%         phNzLevel % in dBc/Hz
%%         phNzFreqOff % in Hz
%%          ampGain - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~exist('x'))
    M = 4;
    data = randi([0 (M - 1)],1000,1);
    pindBm = 10; % (dBm) signal power level
    pin = 10.^((pindBm-30)/10); % Convert dBm to linear Watts
    x = qammod(data,M,'UnitAveragePower',true)*sqrt(pin);
    if(flag_plot)
        scatterplot(x); title('Input Signal')
    end
end
if(~exist('phaseOffset'))
    phaseOffset = 30;
end
if(~exist('freqOffset'))
    freqOffset = 50;
end
if(~exist('phNzFreqOff'))
    phNzFreqOff = -50;
end
if(~exist('IQ_balance'))
    IQ_balance = [3, 3, 0.01, -0.02];
end
if(~exist('ampGain'))
    ampGain = 20;
end
%% Phase Offset  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference: https://www.mathworks.com/help/comm/ref/comm.phasefrequencyoffset-system-object.html
po = comm.PhaseFrequencyOffset('PhaseOffset',phaseOffset);
x_po = po(x); %%%%%%%%%%%%<------------------
if(flag_plot)
    scatterplot(x_po); title('Phase Offset')
end

%% Freq Offset
fo = comm.PhaseFrequencyOffset('FrequencyOffsetSource','Property','FrequencyOffset',freqOffset);
x_fo = fo(x_po); %%%%%%%%%%%%<------------------
if(flag_plot)
    scatterplot(x_fo); title('Freq Offset')
end

%% freq Offset; phase noise  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reference: https://www.mathworks.com/help/comm/ref/comm.phasenoise-system-object.html


phznoise = comm.PhaseNoise('Level',phNzFreqOff(1),'FrequencyOffset',phNzFreqOff(2));
x_po_phz = phznoise(x_fo); %%%%%%%%%%%%<------------------
if(flag_plot)
    scatterplot(x_po_phz); title('Phase Noise Freq Offset')
end


%% I/Q imbalance  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference https://www.mathworks.com/help/comm/ref/iqimbalance.html#mw_c204123f-e047-4b26-9198-707597ad0063

Ia = IQ_balance(1); %3; %[dB] Amplitude imbalance
Ip = IQ_balance(1); %3; %[degrees] Phase imbalance
IDC = IQ_balance(3); %0.01; %[?] DC offset applied to I
QDC = IQ_balance(4); %-0.02; %[?] DC offset applied to IQ
x_iq_in = x_po_phz; %%%%%%%%%%%%<------------------
tmp_real = real(x_iq_in);
tmp_imag = imag(x_iq_in);
yAI = 10^(0.5*Ia/20)*tmp_real + 1i*10^(-0.5*Ia/20)*tmp_imag; %yAmplitudeImbalance
yPI = exp(1i*(-0.5*pi*Ip/180))*real(yAI) + exp(1i*(pi/2+0.5*pi*Ip/180))*imag(yAI); %yPhaseImbalance
y_DCoffset = (real(yPI) + IDC) + 1i*(imag(yPI) + QDC); % add DC offset to I and Q individually
x_IQimbal = y_DCoffset;  %%%%%%%%%%%%<------------------
if(flag_plot)
    scatterplot(x_IQimbal); title('IQ Imbalance') %figure; plot(real(x_IQimbal),imag(x_IQimbal))
end

%% Amplifier %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% memoryless nonlinearity
% reference: https://www.mathworks.com/help/comm/ref/comm.memorylessnonlinearity-system-object.html
if(0)
    amplifier = comm.MemorylessNonlinearity("Method","Cubic polynomial", "LinearGain",ampGain,"AMPMConversion",0);%,"ReferenceImpedance",50);
elseif(0)
    %% nonlinear amplifier gain compression
    linear_gain = 1;
    amplifier = comm.MemorylessNonlinearity('IIP3',30,'LinearGain',linear_gain,'Method','Rapp model','OutputSaturationLevel',2);
elseif(0)
    IIP3_level = ampGain;
    amplifier = comm.MemorylessNonlinearity('IIP3',ampGain);
else
    linear_gain = ampGain; %[dB]
    amplifier = comm.MemorylessNonlinearity('LinearGain',linear_gain,'Method','Rapp model','OutputSaturationLevel',2);
end
    
x_amp = amplifier(x_IQimbal);
if(flag_plot)
    scatterplot(x_amp); title('Amplifier Output')
end
if(0) %debug
    snr = 25;
    noisyLinOut = awgn(x_amp(:,1),snr,"measured");
end

if(0)  %debug 
    x_amp_test = amplifier(x);
    scatterplot(x);
    scatterplot(x_amp_test)
    Fs = 100;
    plotFFT(x,Fs)
    plot(x_amp_test,Fs)
end

x_out = x_amp;
if(flag_plot)
    scatterplot(x_out)
    title('Scatter plot - x-out')
end
end



