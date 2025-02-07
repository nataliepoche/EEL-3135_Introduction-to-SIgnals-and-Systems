%% ACKNOWLEDGEMENTS / REFERENCES: 
% This code uses functions written by Ken Schutte in 2019, which is used to
% read and decode midi files. The code is under a GNU General
% Public License, enabling us to run, study, share, and modify the
% software. 
%
% More info can be found at: http://www.kenschutte.com/midi


%% INITIAL SETUP
clear
close all
clc

%% DEFINE MUSIC

% INTITIAL VARIABLES
Fs = 44100;             % ==> This is the cyclic frequency sample that is the audio signal <==

%  ===> 
% readmidi loads the gym.mid file and midiInfo processes those results, 
% 0 and two specify how it should be decoded. The results come out with Notes 
% being the parsed data, and endtime holds the ened time of the data 
% <===
[Notes, endtime] = midiInfo(readmidi('gym.mid'), 0, 2);
L = size(Notes,1);      % ==> 
                        % Calculates the number of rows in the Notes array 
                        % that is tken from the MIDI file, which means that it's 
                        % the total number of notes in the MIDI file 
                        % <==

%  ===> 
% Calls the build_song function where
%   ones(L,1) is the amplitude determined by the number of notes in the MIDI File 

%   Notes (:,) is the key numberfor the notes in the MIDI file to determine 
%   the corresponding key note

%   Notes(:,6)-Notes(;,5) is how long each note is by calculating the difference 
%   between end time in column 6 and start time in column 5 of the notes array

%   Fs is the cyclic sampling frequency that defines the rate the audio signal is sampled 
% <===
x = build_song(ones(L,1), Notes(:,3), Notes(:,6)-Notes(:,5), Fs);

%  ===> 
% Calculates the total samples needed to represent the entire audio signal
% The total duration of the song is obtained by adding all the durations of 
% the notes taken by subtracting end time with start time, and multiplying 
% that sum by the Fs or samples per second (cyclic sample frequency)
% <===
tot_samples = ceil(sum(Notes(:,6)-Notes(:,5))*Fs);

%  ===> 
% This line creates the time vector from 0 to the total duration of the song 
% based on total number of smaples (tot_samples) and the cyclic sampling frequency (Fs)
% <===
t = 0:1/Fs:(tot_samples-1)/Fs;  % ==> Step size is 1/Fs so that the time vector 
                                % corresponds with each audio signal sample. The 
                                % resulting t will be used to plot audio waveform 
                                % <==

figure(1);                      % Creates a ficure to put the image in
subplot(211)                    % Creates a subplot in the figure(1)
plot(t, x);                     % Plots audio signal x against time vector t the subplot(211)
xlabel('Time [s]')              % Labels x-axis of subplot(211)
ylabel('Amplitude')             % Labels y-axis of subplot(211)
subplot(212)                    % Creaes a second subplot under the first in the figure(1)
plot(t, x);                     % Plots audio signal x against time vector t the subplot(212)
xlabel('Time [s]')              % Labels x-axis of subplot(212)
ylabel('Amplitude')             % Labels y-axis of subplot(212)
axis([0 0.1 -1 1])              % ==> Sets parameters of the x-axis to 0.1 
                                % seconds and the y-axis to the range [-1,1], 
                                % it basically zooms in on the waveform for 
                                % better visualization of the start of the song 
                                % <==

%  ===> 
% The line pauses to wait for the user to press a button before playing 
% a sound inorder to examine the waveform before hearing the audio 
% <===
input('Click any button to play sound')
soundsc(x, Fs);                 % Plays the audio signal



% ========
% YOU DO *NOT* NEED TO DESCRIBE THESE LINES (your free to figure it out though)
W = 0.1;    % Window size
tic;
for mm = 1:ceil(tot_samples/Fs/W)
    % PAUSE UNTIL NEXT FRAME
    xlim([(mm-1)*W+[0 W]]); % Set limits of plot
    tm = toc;                        % Check current time
    if mm*W < tm, disp(['Warning: Visualization is ' num2str(mm*W-tm)  's behind']); end
    drawnow; pause(mm*0.1-tm);       % Synchronize with clock
end
% =======



%%
% =========================================
% SUPPORTING FUNCTIONS FOUND BELOW
% Add comments appropriately below
% =========================================


function x = key_to_note(A, key, dur, fs)
% key_to_note: ========> Creates the sinusoid waveform of a single note <=========
%
% Input Args:
%     A: complex amplitude
%   key: number of the note on piano keyboard
%   dur: duration of each note (in seconds)
%    fs: A scalar sampling rate value
%
% Output:
%     x: sinusoidal waveform of the note
    
    %  ===> Takes apart the components of the note based on the MIDI key number <===
    N    = floortol(dur*fs);            % Calculates the number of samples needed for each note using floortol to prevent floating points
    t    = (0:(N-1)).'/fs;              % Creates a time vector that spans the duration of the note, t corresponds to each sample point in 
                                        % note duration
    freq = (440/32)*2^((key-9)/12);     % Calulates frequency of note based on MIDI key number, 440 Hz corresponds to note A4 or key number 49, 2^(key - 9) / 12) 
                                        % shifts frequency based on MIDI key number, using equal-tempered scale where each key is a half-step apart
    
    %  ===> Generates sinusoidal waveform for the note <===
    x    = real(A*exp(1j*2*pi*freq*t));   


end


function x = build_song(As, keys, durs, fs)
% build_song:  ========> 
% Creates the full audio signal by placing each individual note int he correct 
% position based on the start time and duration after calling key_to_note to generate 
% waveform for each note 
% <=========
%
% Input Args:
%	  As: A length-N array of complex amplitudes for building notes
%	keys: A length-N array of key numbers (which key on a keyboard) for building notes
%   durs: A length-N array of durations (in seconds) for building notes
%     fs: A scalar sampling rate value
%
% Output Args: 
%      x: A length-(N*fs) length raw audio signal
%
    %  ===> 
    % Initializes the audio signal with zeros the same size as all the 
    % song notes by summing all note durations times the cyclic sampling rate (Fs) 
    % <===
    x = zeros(ceil(sum(durs)*fs), 1);      
    for k = 1:length(keys) 
        
        %  ===> 
        % note generates the waveform of the single note, start_time calculates 
        % the total start time for the current note by calculating all the durations 
        % of previous notes (basically calculates when the note should begin in the final audio 
        % <===
        note       = key_to_note(As(k), keys(k), durs(k), fs);  
        start_time = sum(durs(1:k-1));
        
        %  ===> 
        % This line calulates sample indices of the current note, n1 is 
        % start and n2 is end, flortol rounds to avoid float point errors 
        % <===
        n1         = floortol(start_time*fs) + 1;               
        n2         = floortol(start_time*fs) + floortol(durs(k)*fs);
        x(n1:n2)   = x(n1:n2) + note;                                       % This line places the generated waveform note into the correct position 
                                                                            % in the audio signal x. The note is added to x between the indices n1 and 
                                                                            % 2 so each note is placed at the correct time in the final audio signal              
        
    end

end

function x = floortol(x)
%FLOORTOL Apply floor operation after adding 0.5 to ensure no
%   floating-point rounding errors that unintendedly decrease the 
%   value
%

    x = floor(x+0.5);

end