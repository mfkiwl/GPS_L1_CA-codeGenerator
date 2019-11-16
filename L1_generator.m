%I watched a video series about C/A code generation
%Polynomials 1 & 2 Combined
%https://www.youtube.com/watch?v=EZV98DgFNz4

%%
%%STEP 1 - Generate PRN(C/A) Codes
%"Polynomials"
pol_1_array = ones(1,10);
pol_2_array = ones(1,10);
res_SV_1 = getCombined_SV1(pol_1_array,pol_2_array); %%PRN 1
pol_1_array = ones(1,10);
pol_2_array = ones(1,10);
res_SV_2 = getCombined_SV2(pol_1_array,pol_2_array); %%PRN 2

prn_1 = repmat(res_SV_1,1,40); %%PRN 1 for 40ms
prn_2 = repmat(res_SV_2,1,40); %%PRN 2 for 40ms

prn_1_kk=kron(prn_1, ones(1,5));  %To be able to analyse the signal well
prn_2_kk=kron(prn_2, ones(1,5));

%%
%%STEP 2 - Generate NAV DATA (Fake)
%50 Hz -> 50e-3 
time= linspace(0,40920/1023,40920*5);
nav_data_signal = square(pi*50e-3*time);
nav_data_signal(nav_data_signal==-1)=0; %Converts the signal into a binary signal

%%
%%STEP 3.1 - NAV + C/A (SV 1)
nav_and_ca_code = prn_1_kk + nav_data_signal;
nav_and_ca_code=mod(nav_and_ca_code,2);

figure('name','L1 Signal of SV 1 (40ms)','numbertitle','off');
subplot(3,1,1);
plot(time,nav_data_signal);
xlabel('Time'); 
ylabel('Amplitude'); 
title('NAV Data(Fake)');
h = zoom;
h.Motion = 'horizontal';
h.Enable = 'on';

subplot(3,1,2);
plot(time,prn_1_kk);
xlabel('Time'); 
ylabel('Amplitude'); 
title('PRN Data');
h = zoom;
h.Motion = 'horizontal';
h.Enable = 'on';

subplot(3,1,3);
plot(time,nav_and_ca_code);
xlabel('Time'); 
ylabel('Amplitude'); 
title('PRN + NAV Data');
h = zoom;
h.Motion = 'horizontal';
h.Enable = 'on';

%%
%%Step 4. BPSK Modulation -------------------------------------------------
%%https://www.mathworks.com/matlabcentral/fileexchange/30582-binary-phase-shift-keying
[PSK_signal,org_data,time_bpsk] = BPSK_Modulation(prn_1(1:1023));

figure('name','L1 Modulated Signal of SV 1 (1ms)','numbertitle','off');
subplot(2,1,2);
plot(time_bpsk,PSK_signal,'LineWidth',2);
xlabel('Time (bit period)');
ylabel('Amplitude');
title('PSK Signal with two Phase Shifts');
axis([0 time_bpsk(end) -1.5 1.5]);
h = zoom;
h.Motion = 'horizontal';
h.Enable = 'on';
grid  on;

% Plot the Original Digital Signal
subplot(2,1,1);
plot(time_bpsk,org_data,'r','LineWidth',2);
xlabel('Time (bit period)');
ylabel('Amplitude');
title('Original Digital Signal');
axis([0 time_bpsk(end) -0.5 1.5]);
h = zoom;
h.Motion = 'horizontal';
h.Enable = 'on';
grid on;

%%Spectrum of L1 signal that modulated using BPSK--------------------------
sampling_freq_range = linspace(-50e6,50e6,1001);
L1_spectrum = spectrum_BPSK(1.023e6,sampling_freq_range);
figure('name','PSD of Modulated L1 Signal','numbertitle','off');
plot(sampling_freq_range/1e6,10*log10(L1_spectrum));
ylabel('PSD (dBW/Hz)');
xlabel('Frequency (MHz)');
grid on;

%%
%%-----------------SV 2 Signals -------------------------------------------
%%-------------------------------------------------------------------------
%%STEP 3.2 - NAV + C/A (SV 2)
%%
nav_and_ca_code_sv2 = prn_2_kk + nav_data_signal;
nav_and_ca_code_sv2=mod(nav_and_ca_code_sv2,2);

figure('name','L1 Signal of SV 2 (40ms)','numbertitle','off');
subplot(3,1,1);
plot(time,nav_data_signal);
xlabel('Time'); 
ylabel('Amplitude'); 
title('NAV Data(Fake)');
h = zoom;
h.Motion = 'horizontal';
h.Enable = 'on';

subplot(3,1,2);
plot(time,prn_2_kk);
xlabel('Time'); 
ylabel('Amplitude'); 
title('PRN Data');
h = zoom;
h.Motion = 'horizontal';
h.Enable = 'on';

subplot(3,1,3);
plot(time,nav_and_ca_code_sv2);
xlabel('Time'); 
ylabel('Amplitude'); 
title('PRN + NAV Data');
h = zoom;
h.Motion = 'horizontal';
h.Enable = 'on';

%%
%%Step 4.2 BPSK Modulation ------------------------------------------------
%%https://www.mathworks.com/matlabcentral/fileexchange/30582-binary-phase-shift-keying
[PSK_signal_sv2,org_data_sv2,time_bpsk] = BPSK_Modulation(prn_2(1:1023));

figure('name','L1 Modulated Signal of SV 2 (1ms)','numbertitle','off');
subplot(2,1,2);
plot(time_bpsk,PSK_signal_sv2,'LineWidth',2);
xlabel('Time (bit period)');
ylabel('Amplitude');
title('PSK Signal with two Phase Shifts');
axis([0 time_bpsk(end) -1.5 1.5]);
h = zoom;
h.Motion = 'horizontal';
h.Enable = 'on';
grid  on;

% Plot the Original Digital Signal
subplot(2,1,1);
plot(time_bpsk,org_data_sv2,'r','LineWidth',2);
xlabel('Time (bit period)');
ylabel('Amplitude');
title('Original Digital Signal');
axis([0 time_bpsk(end) -0.5 1.5]);
h = zoom;
h.Motion = 'horizontal';
h.Enable = 'on';
grid on;

%%Spectrum of L1 signal that modulated using BPSK--------------------------
sampling_freq_range = linspace(-50e6,50e6,1001);
L1_spectrum = spectrum_BPSK(1.023e6,sampling_freq_range);
figure('name','PSD of Modulated L1 Signal','numbertitle','off');
plot(sampling_freq_range/1e6,10*log10(L1_spectrum));
ylabel('PSD (dBW/Hz)');
xlabel('Frequency (MHz)');
grid on;
%%

%%
%%Correlation
figure('name','PRN Codes Correlations','numbertitle','off');
prn_1_corr = xcorr(res_SV_1,'normalized');
% prn_1_corr = autocorr(res_SV_1,1);
subplot(3,1,1);
plot(prn_1_corr);
title('Auto Correlation - SV 1');
k = zoom;
k.Motion = 'horizontal';
k.Enable = 'on';

prn_2_corr = xcorr(res_SV_2,'normalized');
subplot(3,1,2);
plot(prn_2_corr);
title('Auto Correlation - SV 2');
k = zoom;
k.Motion = 'horizontal';
k.Enable = 'on';

prn_1_prn_2_corr = xcorr(res_SV_2,res_SV_1,'normalized');
subplot(3,1,3);
plot(prn_1_prn_2_corr);
ylim([1 2]);
title('Cross Correlation');
%%

%%Carrier Freq L1 = 1575.42 MHz but we use IF instead of RF
%%Rf to IF conversion make the progress faster

%%--------------------------------------%%---------------------------------
%%FUNCTIONS
function [output, pol_array] = pol_1(pol_array)
    %EQ - POL 1 = 1 + x(3) + x(10)
    %shift and get 1, 3, 10th elements of the array
    %and sum them and assign it to the first element of the array in the
    %next iteration       
    output = pol_array(end);
    first_element_calculated = mod(pol_array(3)+pol_array(10),2);
    pol_array=[first_element_calculated,pol_array(1:end-1)];
end

function [output, pol_array] = pol_2_SV_No_1(pol_array)
    %EQ - POL 2 = 1 + x(2) + x(3) + x(6) + x(8) + x(9) + x(10)
    %This function is a unique for each satellite, in this example we will
    %focus on SV 1 (Satellite identifier), C/A Code Selection for this sat 2 and 3rd bits. So we
    %shift the array and get 2nd and 3rd bits to sum and assign the summing
    %result as an output. For the first element of the array in the next
    %iteration will be sum of 2, 3, 6, 8, 9 and 10th elements.
    
    %This function can be generalized by adding one parameter showing the selection bits
    %for using the same funtion for all satellites.
    
    output = mod(pol_array(2)+pol_array(6),2);
    first_element_calculated = mod(pol_array(2)+pol_array(3)+pol_array(6)+pol_array(8)+pol_array(9)+pol_array(10),2);
    pol_array=[first_element_calculated,pol_array(1:end-1)];
end

function [output, pol_array] = pol_2_SV_No_2(pol_array)
    %EQ - POL 2 = 1 + x(2) + x(3) + x(6) + x(8) + x(9) + x(10)
    %This function is a unique for each satellite, in this example we will
    %focus on SV 1 (Satellite identifier), C/A Code Selection for this sat 3 and 7th bits. So we
    %shift the array and get 2nd and 3rd bits to sum and assign the summing
    %result as an output. For the first element of the array in the next
    %iteration will be sum of 2, 3, 6, 8, 9 and 10th elements.
    
    %This function can be generalized by adding one parameter showing the selection bits
    %for using the same funtion for all satellites.
    
    output = mod(pol_array(3)+pol_array(7),2);
    first_element_calculated = mod(pol_array(2)+pol_array(3)+pol_array(6)+pol_array(8)+pol_array(9)+pol_array(10),2);
    pol_array=[first_element_calculated,pol_array(1:end-1)];
end

function result_array = getCombined_SV1(pol_1_ar,pol_2_ar)
    result_array=zeros(1,10);
    for n=1:1023
        [output_pol_1_ar, pol_1_ar] = pol_1(pol_1_ar);
        [output_pol_2_ar, pol_2_ar] = pol_2_SV_No_1(pol_2_ar);
        result_array(n)=mod(output_pol_1_ar+output_pol_2_ar,2);
    end
end

function result_array = getCombined_SV2(pol_1_ar,pol_2_ar)
    result_array=zeros(1,10);
    for n=1:1023
        [output_pol_1_ar, pol_1_ar] = pol_1(pol_1_ar);
        [output_pol_2_ar, pol_2_ar] = pol_2_SV_No_2(pol_2_ar);
        result_array(n)=mod(output_pol_1_ar+output_pol_2_ar,2);
    end
end

function out = autocorr(indata, tn)

if nargin < 2
tn = 1;
end

ln = length(indata);
out = zeros(1,ln*tn);

for ii=0:ln*tn-1
out(ii+1) = sum(indata.*shift(indata,ii,0));
end

end

function outregi = shift(inregi,shiftr,shiftu)
[h, v] = size(inregi);
outregi = inregi;

shiftr = rem(shiftr,v);
shiftu = rem(shiftu,h);

if shiftr > 0
outregi(:,1 :shiftr) = inregi(:,v-shiftr+1:v );
outregi(:,1+shiftr:v ) = inregi(:,1 :v-shiftr);
elseif shiftr < 0
outregi(:,1 :v+shiftr) = inregi(:,1-shiftr:v );
outregi(:,v+shiftr+1:v ) = inregi(:,1 :-shiftr);
end

inregi = outregi;

if shiftu > 0
outregi(1 :h-shiftu,:) = inregi(1+shiftu:h, :);
outregi(h-shiftu+1:h, :) = inregi(1 :shiftu,:);
elseif shiftu < 0
outregi(1 :-shiftu,:) = inregi(h+shiftu+1:h, :);
outregi(1-shiftu:h, :) = inregi(1 :h+shiftu,:);
end

end

function out = crosscorr(indata1, indata2, tn)

if nargin < 3
tn = 1;
end

ln = length(indata1);
out = zeros(1,ln*tn);

for ii=0:ln*tn-1
out(ii+1) = sum(indata1.*shift(indata2,ii,0));

end

end

%%BPSK spectrum
%I have a lot of lack of knowledge about this part, after I study on that
%topic maybe I need to update this part.
function S = spectrum_BPSK(fc,F)
%==========================================================================
% BPSK of chipping rate fc.
%--------------------------------------------------------------------------
% Author: Daniel Pascual (daniel.pascual [at] protonmail.com) 
% Copyright 2017 Daniel Pascual
% License: GNU GPLv3
%==========================================================================

    S = sinc((1/fc)*F); 
    aux = sum(abs(S).^2)/length(S); 
    S = S/sqrt(aux);
    S = abs(S).^2;
end

function [mod_sig, org_sig, time_out] = BPSK_Modulation(data)
    %%For BPSK Modulation
    %%https://www.mathworks.com/matlabcentral/fileexchange/30582-binary-phase-shift-keying
    % Enter the two Phase shifts - in Radians
    % Phase for 0 bit
    P1 = 0; 
    % Phase for 1 bit
    P2 = pi;
    freq_L1=20; %Instead of using L1 real frequency, we use lower frequency
    %to make the calculation easier and make the details more appear on the
    %graph
    freq_Samp=freq_L1*1e1; %Sampling frequency have to be at least double the
    %carrier frequency. Higher sampling, more resolution.
    f=freq_L1;
    t=0:1/freq_Samp:1;
    time=[];
    PSK_signal = [];
    Digital_signal = [];

    for ii = 1: 1: length(data)

        % The FSK Signal
        PSK_signal = [PSK_signal (data(ii)==0)*cos(2*pi*f*t + P1)+...
            (data(ii)==1)*cos(2*pi*f*t + P2)];

        % The Original Digital Signal
        Digital_signal = [Digital_signal (data(ii)==0)*...
            zeros(1,length(t)) + (data(ii)==1)*ones(1,length(t))];

        time = [time t];
        t =  t + 1;        

    end
    
    mod_sig = PSK_signal;
    org_sig = Digital_signal;
    time_out = time;
end