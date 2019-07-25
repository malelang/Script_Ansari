function r = wabp(Araw, Offset,Scale, Fs)
% r = wabp(Araw,Offset,Scale, Fs);
% Input: Araw (125Hz sampled) waveform in wfdb-MIT format, 
%        Offset, Scale
% Output: The onset times of the input ABP waveform
% Defaults are Offset = 1600; Scale=20; Fs=125; 
% If you pass 0 as the scale then above defaults are invoked
% Default Fs=125Hz unless you pass a fourth argument.
%
% Gnu Public License Applies
% 
% James Sun Feb 09 2005 with some changes from Gari Clifford
% based upon wabp.c by Wei Zong (www.physionet.org)

% if the signal is not at 125Hz, then resample it.
if nargin<4
Fs=125;
end

if nargin < 3
Scale = 20;
end

if nargin<2
Offset = 1600;
end

if Scale==0
Scale = 20;
Offset = 1600;
end

% if the sample frequency is not 125, resample to 125
if Fs~=125
Q=round(Fs);
P=round(125);
Araw = resample(Araw, P, Q);
end



%%%%%%%%%%%%%%%%FILTRO ORIGINAL QUE VIENE POR DEFECTO
% %%LPF 
% A = filter([1 0 0 0 0 -2 0 0 0 0 1],[1 -2 1],Araw)/24+30;
% A = (A+Offset)/Scale;
% A = A(4:end);  % Takes care of 4 sample group delay
% %Slope-sum function ... not used?
% x = zeros(size(A));
%%%%%%%%%%%%%%  SIN NIGUN TIPO DE FILTRO
 % A = (Araw+Offset)/Scale;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db8 level 3%%%%%%%%%%%%%
% [C,L] = wavedec(Araw,3,'db8'); 
% A3 = wrcoef('a',C,L,'db8',3); % mejor linea base
% cleanedSignal = detrend(A3);
% A = (cleanedSignal+Offset)/Scale;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db6 level 4%%%%%%%%%%%%%
% [C,L] = wavedec(Araw,4,'db6'); 
% A3 = wrcoef('a',C,L,'db6',4); % mejor linea base
% cleanedSignal = detrend(A3);
% A = (cleanedSignal+Offset)/Scale;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db10 level 4%%%%%%%%%%%%%
[C,L] = wavedec(Araw,4,'db10'); 
A3 = wrcoef('a',C,L,'db10',4); % mejor linea base
cleanedSignal = detrend(A3);
A = (cleanedSignal+Offset)/Scale;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets db6 level 5%%%%%%%%%%%%%
% [C,L] = wavedec(Araw,5,'db6'); 
% A3 = wrcoef('a',C,L,'db6',5); % mejor linea base
% cleanedSignal = detrend(A3);
% A = (cleanedSignal+Offset)/Scale;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets sym4 level 3%%%%%%%%%%%%%
% [C,L] = wavedec(Araw,3,'sym4'); 
% A3 = wrcoef('a',C,L,'sym4',3); % mejor linea base
% cleanedSignal = detrend(A3);
% A = (cleanedSignal+Offset)/Scale;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelets sym6 level 3%%%%%%%%%%%%%
% [C,L] = wavedec(Araw,3,'sym6'); 
% A3 = wrcoef('a',C,L,'sym6',3); % mejor linea base
% cleanedSignal = detrend(A3);
% A = (cleanedSignal+Offset)/Scale;
%%%%%%%%%%%%%%%%%%%%%%%%%%Savitzky-Golay Smoothing Filter%%%%%
%   A3=sgolayfilt(Araw,3,41);
%   cleanedSignal = detrend(A3);
%   A = (cleanedSignal+Offset)/Scale;
 %%%%%%%%%%%%%%%%%%%%%%%%%%Mov median Filter%%%%%        
%     A3  = movmedian(Araw,11);
%     cleanedSignal = detrend(A3);
%     A = (cleanedSignal+Offset)/Scale;
%%%%%%%%%%%%%%%%%%%%%%%%%%% wavelets with thresholding %%%%%%%%%
%     [C,L] = wavedec(Araw,5,'db6'); 
%    [thr,sorh,keepapp]=ddencmp('den','wv',Araw);
%    A3=wdencmp('gbl',C,L,'db6',5,thr,sorh,keepapp);
%    cleanedSignal = detrend(A3);
%    A = (cleanedSignal+Offset)/Scale;
%%%%%%%%%%%%%%%%%%%%%%%%%%% EMD  %%%%%%%%%%%%%%%%%%%%%
%   cleanedSignal = emd_dfadenoising(Araw);
%   cleanedSignal = cleanedSignal';
%   A3=detrend(cleanedSignal); 
%   A = (A3+Offset)/Scale;

dyneg = [A' 0] - [0 A'];
dyneg(find(dyneg>0)) = 0;

dypos = [A' 0] - [0 A'];
dypos(find(dypos<0)) = 0;
h = ones([16 1]);
ssf = conv(h,dypos);
ssf = [0 0 ssf]'; %'
%plot(ssf);


% Decision rule
avg0 = sum(ssf(1:1000))/1000;   % average of 1st 8 seconds (1000 samples) of SSF 
Threshold0 = 3*avg0;            % initial decision threshold 

% ignoring "learning period" for now
lockout = 0;    % lockout >0 means we are in a refractory period
timer = 0;
BeatTime = 0;
z=zeros(1,100000);
z1=z;
counter=0;

for t= 50:length(ssf)-17
    lockout = lockout -1;
    timer = timer + 1;      % Timer used for counting time after previous ABP pulse
    
    if (lockout < 1) & (ssf(t) > avg0+5)       % Not in refractory and SSF has exceeded threshold here
        timer = 0;
        maxSSF = max(ssf(t:t+16));  % Find local max of SSF
        minSSF = min(ssf(t-16:t));  % Find local min of SSF
        if maxSSF > (minSSF + 10)
            onset = 0.01*maxSSF ;  % Onset is at the time in which local SSF just exceeds 0.01*maxSSF
            
            tt = t-16:t;
            dssf = ssf(tt) - ssf(tt-1);
            BeatTime = max(find(dssf < onset))+t-17;
            counter = counter+1;

            if isempty(BeatTime)
                counter = counter-1;
            else
            z(counter) = BeatTime;
        end

            Threshold0 = Threshold0 + 0.1*(maxSSF - Threshold0);  % adjust threshold
            avg0 = Threshold0 / 3;        % adjust avg
            
            lockout = 32;   % lock so prevent sensing right after detection (refractory period)
        end
    end
    
    if timer > 312  % Lower threshold if no pulse detection for a while
        Threshold0 = Threshold0 -1;
        avg0 = Threshold0 / 3;
    end
end

r = z(find(z))'; %'
r = r-2;         % adjust for filter delay?