clear all
close all
clc
% % % % %
% Wireless Receivers: algorithms and architectures
% Audio Transmission Framework 
%
%
%   3 operating modes:
%   - 'matlab' : generic MATLAB audio routines (unreliable under Linux)
%   - 'native' : OS native audio system
%       - ALSA audio tools, most Linux distrubtions
%       - builtin WAV tools on Windows 
%   - 'bypass' : no audio transmission, takes txsignal as received signal

%fsym = [100,400,800,1200,1600,2000,]%100:600:2000
%for jj = 1:length(fsym);

    % Configuration Values
    conf.audiosystem = 'matlab'; % Values: 'matlab','native','bypass'

    conf.do_phase_estim = 1;
    conf.do_phase_track = 1;
    
    conf.n_carriers  = 256; % Number of OFDM carriers, multiple of 2 for efficient FFT implementation
    conf.f_s     = 48000;   % sampling rate, f_sampling = N/T with T the OFDM symbol length without prefix? 
    conf.f_spacing = 5; % In hz, f_spacing = 1/T, with T defined as above
    conf.os_factor  = ceil(conf.f_s/(conf.f_spacing*conf.n_carriers)); %os factor for OSIFFT & IFFT also used for preamble
    conf.f_spacing = conf.f_s/( conf.os_factor*conf.n_carriers)
   % conf.f_sym   = conf.n_carriers * conf.f_spacing;     % symbol rate, f_spacing = 1/T 
    conf.f_bw = ceil(( conf.n_carriers +1 )/2 )*conf.f_spacing;
    conf.nframes = 1;       % number of frames to transmit - dit was gegeven en is verwarrend
    conf.modulation_order = 2; % BPSK:1, QPSK:2
    conf.f_c     = 8000;
    conf.offset = 0;
    conf.nbits = 256*200;
    
    conf.data_length = conf.nbits/conf.modulation_order;
    conf.n_data_frames = ceil(conf.data_length/conf.n_carriers);
    
    % f_spacing must be divisor of f_sampling
    
    conf.npreamble  = 100;
    conf.bitsps     = 16;   % bits per audio sample
    conf.cpref_length = 0.5; % cyclic prefix length is half of the OFDM symbol length
    
    

    % Training symbols
    conf.training_symbols = (1 - 2 * lfsr_framesync(conf.n_carriers));
    conf.training_interval = 200;
    
    % Init Section
    % all calculations that you only have to do once
    if mod(conf.os_factor,1) ~= 0
       disp('WARNING: Sampling rate must be a multiple of the symbol rate'); 
    end
    conf.nsyms      = ceil(conf.nbits/conf.modulation_order);

    % Initialize result structure with zero
    res.biterrors   = zeros(conf.nframes,1);
    res.rxnbits     = zeros(conf.nframes,1);

    % To speed up your simulation pregenerate data you can reuse
    % beforehand.
    % Create RRC pulse
    tx_filterlen = 10;
    rolloff_factor = 0.22;
    conf.h = rrc(conf.os_factor, rolloff_factor, conf.os_factor*tx_filterlen);


    % Results


    for k=1:conf.nframes

        % Generate random data
        txbits = randi([0 1],conf.nbits,1);
        %txbits = ones(conf.nbits, 1);
        % TODO: Implement tx() Transmit Function
        [txsignal, conf] = tx(txbits,conf,k);
        
        % % % % % % % % % % % %
        % Begin
        % Audio Transmission
        %

        % normalize values
        peakvalue       = max(abs(txsignal));
        normtxsignal    = txsignal / (peakvalue + 0.3);

        
        % create vector for transmission
        rawtxsignal = [ zeros(conf.f_s,1) ; normtxsignal ;  zeros(conf.f_s,1) ]; % add padding before and after the signal
        rawtxsignal = [  rawtxsignal  zeros(size(rawtxsignal)) ]; % add second channel: no signal
        txdur       = length(rawtxsignal)/conf.f_s; % calculate length of transmitted signal

        % wavwrite(rawtxsignal,conf.f_s,16,'out.wav')   
        audiowrite('out.wav',rawtxsignal,conf.f_s)  

        % Platform native audio mode 
        if strcmp(conf.audiosystem,'native')

            % Windows WAV mode 
            if ispc()
                disp('Windows WAV');
                wavplay(rawtxsignal,conf.f_s,'async');
                disp('Recording in Progress');
                rawrxsignal = wavrecord((txdur+1)*conf.f_s,conf.f_s);
                disp('Recording complete')
                rxsignal = rawrxsignal(1:end,1);

            % ALSA WAV mode 
            elseif isunix()
                disp('Linux ALSA');
                cmd = sprintf('arecord -c 2 -r %d -f s16_le  -d %d in.wav &',conf.f_s,ceil(txdur)+1);
                system(cmd); 
                disp('Recording in Progress');
                system('aplay  out.wav')
                pause(2);
                disp('Recording complete')
                rawrxsignal = wavread('in.wav');
                rxsignal    = rawrxsignal(1:end,1);
            end

        % MATLAB audio mode
        elseif strcmp(conf.audiosystem,'matlab')
            disp('MATLAB generic');
            playobj = audioplayer(rawtxsignal,conf.f_s,conf.bitsps);
            recobj  = audiorecorder(conf.f_s,conf.bitsps,1);
            record(recobj);
            disp('Recording in Progress');
            playblocking(playobj)
            pause(0.5);
            stop(recobj);
            disp('Recording complete')
            rawrxsignal  = getaudiodata(recobj,'int16');
            rxsignal     = double(rawrxsignal(1:end))/double(intmax('int16')) ;

        elseif strcmp(conf.audiosystem,'bypass')
            rawrxsignal = rawtxsignal(:,1);
            rxsignal    = rawrxsignal;
        end

        % Plot received signal for debgging
        figure;
        plot(rxsignal);
        title('Received Signal')

        %
        % End
        % Audio Transmission   
        % % % % % % % % % % % %

        [rxbits, conf]       = rx(rxsignal,conf);

       % plot([txbits rxbits])
        
        res.rxnbits(k)      = length(rxbits);
        res.biterrors(k)    = sum(rxbits ~= txbits);

    end

   per = sum(res.biterrors > 0)/conf.nframes
   ber = sum(res.biterrors)/sum(res.rxnbits)

%end
 
%  figure('Name', 'BER');
%  semilogy(fsym, ber);
%  xlabel('Symbol rate');
%  ylabel('BER [%]');
%  title('BER with time offset at receiver equal to 0.5Hz');
