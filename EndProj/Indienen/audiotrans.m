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

    % Configuration Values
    conf.audiosystem = 'bypass'; % Values: 'matlab','native','bypass'

    conf.do_phase_estim = 1; % Do phase estimation?
    conf.do_phase_track = 1; % Do phase tracking?
    
    conf.n_carriers  = 256; % Number of OFDM carriers, multiple of 2 for efficient FFT implementation
    conf.f_s     = 48000;   % sampling rate, f_sampling = Nos/T
    conf.f_spacing = 5; % will be modified below so it becomes a value close to this estimation but one that makes OS factor integer
    conf.os_factor  = ceil(conf.f_s/(conf.f_spacing*conf.n_carriers)); %normally, formula for the os factor is without the ceil, but must be integer
    conf.f_spacing = conf.f_s/( conf.os_factor*conf.n_carriers) % update fspacing so oversampling factor is the integer defined abvoe
    conf.f_bw = ceil(( conf.n_carriers +1 )/2 )*conf.f_spacing; % baseband bandwidth of OFDM spectrum
    conf.corner_f = conf.f_bw*1.2; % corner frequency of lowpass
    
    conf.modulation_order = 2; % BPSK:1, QPSK:2
    conf.f_c     = 8000; % carrier frequency
    conf.offset = 0; % simulate oscillator frequency offset in receiver
    
    conf.nframes = 1;    % number of frames to transmit
    conf.nbits = 256*2*20; % number of bits to transmit
    
    conf.npreamble  = 100; % length of preamble
    conf.bitsps     = 16;   % bits per audio sample
    conf.cpref_length = 128/256; % cyclic prefix length is half of the OFDM symbol length
    
    % Training symbols
    conf.training_symbols = (1 - 2 * lfsr_framesync(conf.n_carriers));
    conf.training_interval = 10;
    
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
   
    conf.data_length = conf.nbits/conf.modulation_order;
    conf.n_data_symbols = ceil(conf.data_length/conf.n_carriers);
    
    tx_filterlen = 10; % Create RRC pulse 
    rolloff_factor = 0.22;
    conf.h = rrc(conf.os_factor, rolloff_factor, conf.os_factor*tx_filterlen);


    % Results
    for k=1:conf.nframes

        % Generate random data
        txbits = randi([0 1],conf.nbits,1);
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
        
        %figure, plot(rxbits ~= txbits), title('rxbits ~= txbits');
        
        res.rxnbits(k)      = length(rxbits);
        res.biterrors(k)    = sum(rxbits ~= txbits);

    end

   per = sum(res.biterrors > 0)/conf.nframes
   ber = sum(res.biterrors)/sum(res.rxnbits)