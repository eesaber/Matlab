score =[1,1,5,5,6,6,5,4,4,3,3,2,2,1];
beat = [1,1,1,1,1,1,2,1,1,1,1,1,1,2];
name = 'twinkle';

%Baseline of the requirement
%getmusic(score, beat, name) 

%Adding feature
%getmusic(score, beat, name, 'Balance', [100,0],'Tempo', 240, 'Chord',4) 


function getmusic(score, beat, name, varargin)
	% Usage: 
	% 1. getmusic(score, beat, name) : Generate the base requirement of the  homwwork.
	% 2. getmusic(score, beat, name, 'Tempo', X, 'Balance', [L,R]): Generate
	% the music at tempo X. X should be a non-negative number and default X=60. 
	% Also you can tune the  ratio of the volumme of left and right audio
	% channel. For example [100, 0], right channel will mute.
	% 3. getmusic(score, beat, name, 'Chord', X) : 0 <= X <= 4, 
	%	X=0: no chord	X=1: Piano	X=2: flute	X=3: Violin	
	%	X=4: random spectrum. 

	%% Parse the argument 
	parse_ = inputParser;
	validationFcn_1_ = @(x) validateattributes(x,{'double'},{'nonnegative'}); 
	validationFcn_2_ = @(x) validateattributes(x,{'double'},{'nonnegative'},2); 
	validationFcn_3_ = @(x) validateattributes(x,{'numeric'},{'scalar'});
	addParameter(parse_,'Tempo',60,validationFcn_1_);
	addParameter(parse_,'Balance',[1,1],validationFcn_2_);
	addParameter(parse_,'Chord',0,validationFcn_3_);
	parse(parse_,varargin{:})
	%% Generate part
	f_s = 20000;
	tempo = parse_.Results.Tempo;
	qq = parse_.Results.Balance;
	n = 2*2*pi*[0,262,294,330,349,392,440,494];
	n_sharp = 2*2*pi*[0,277,311,0,370,415,466];
	pu = n(int32(score.* ((mod(score,1) == 0)) + 1)) + n_sharp(int32(floor(score .* (mod(score,1)>0)) + 1));
	instru = [poisspdf(3:11,5); 0 10 10 1 2 1.8 0.1 0.1 0.1;...
		0 10 6 6 7 4 2 4 1; rand(1,9)];
	m_len = sum(int32(beat * 60 / tempo * f_s)) ;
	y = zeros(1,m_len);
	pt = 1;
	f_ch = [0.5, 1:8].';
	if parse_.Results.Chord 
		for i = 1 : length(score)
			t = linspace(0,60/tempo * beat(i),int32(beat(i) * 60 / tempo * f_s)); 
			y(pt:pt+length(t)-1) = instru(parse_.Results.Chord,:) * sin( pu(i) * f_ch * t) ;
			pt = pt + length(t);
		end
	else 
		for i = 1 : length(score)
			t = linspace(0,60/tempo * beat(i),int32(beat(i) * 60 / tempo * f_s)); 
			y(pt:pt+length(t)-1) = sin(pu(i) * t);
			pt = pt + length(t);
		end
	end
	
	y = [qq(1)*y; qq(2)*y ];
	y = y / max(max(y));
	
	audiowrite([name '.wav'],y.',f_s);

end