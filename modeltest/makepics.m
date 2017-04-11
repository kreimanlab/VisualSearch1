NBOBJS = 40;
MAKENATIMAGES=0;
NBCLUT=0;
NBUNPRED =40;
NBDOUB=0;
NBQUAD=0;
NBQUADSIZE=0;
NBCONT=0;
NBARRAY=40;
STIMSIZE = 124;  
MARGIN = 60; %40;
CLUTSTIMSIZEBIG = 220; % Stimuli are larger for cluttered background images
CSIZE = 256;  
ctr = CSIZE/2;

odir = './images/';
natdir = './naturalimages/';
imagedatadir = './imagedata/';

EXTENSION = 'ot';
%fid = fopen('./listOTHER.txt');
%list = textscan(fid, '%s');
%fclose(fid);

zz = dir([imagedatadir '*.mat']);
list = {zz.name};
zz = dir([natdir '*.pgm']);
listnat = {zz.name};

%list(1:10) = [];

%fid = fopen([natdir 'list.txt']);
%listnat = textscan(fid, '%s');
%fclose(fid);

if NBUNPRED > 0
	funpredid= fopen([odir 'unpredpos.txt'], 'w');
	fclose(funpredid);
end;
if NBQUAD > 0
	fquadid= fopen([odir 'quad.txt'], 'w');
	fclose(fquadid);
end;
if NBQUADSIZE > 0
	fsizeid= fopen([odir 'quadsize.txt'], 'w');
	fclose(fsizeid);
end;
if NBARRAY > 0
	farrayid = fopen([odir 'array.txt'], 'w');
	fclose(farrayid);
end;



fullpic=[];
row=[];
nbims = 0; totalims=0;
pics={}; pic=[];

%if MAKENATIMAGES
%	for i=1:numel(listnat)
%
%		disp(['Natural image ' num2str(i)]);
%		fname=listnat{i};
%		background=imread(fname);
%		if ((size(background, 1) ~= CSIZE) || (size(background, 1) ~= CSIZE)) error('Wrong image size!'); end; %Checks if 'natural images' used for background have the right size
%		
%
%		% Contrast maximization
%		z = double(background);
%		z = z-min(z(:)); z = z ./ max(z(:));
%		background=z;
%
%		
%		fname=[odir 'nat' num2str(i) '.png'];	
%		canvas = background;
%		imwrite(canvas, fname); 
%		%imwrite(medfilt2(canvas), fname);     % MEDIAN FILTERING !!!!!
%
%		t=imread(fname);
%		fid = fopen([fname '.raw'], 'w');
%		fwrite(fid, t');
%		fclose(fid);
%
%	end;
%end;


for i=1:NBOBJS
%for i=1:numel(list)
    fname = [imagedatadir list{i}];
	disp(fname);
	z=load(fname); 
	z = z.image;
	maxz = max(abs(z(:)));
	z = z ./ (1e-10+maxz);  % Maximizing contrast within the [-1,1] range. Note that the background, which is at 0, is not affected.
    pic = (round(pic*100)) / 100;
	z = .5 + .5 * z;
	z(1:MARGIN,:)=[]; z(end-MARGIN+1:end,:)=[];
	z(:,1:MARGIN)=[]; z(:,end-MARGIN+1:end)=[];
	origpic = z;
	stim = imresize(z, [STIMSIZE STIMSIZE], 'bicubic');

	fname=[odir num2str(i) '.normal.' EXTENSION '.png'];	
	canvas = .5*ones(CSIZE); 
	canvas(ctr-STIMSIZE/2:ctr+STIMSIZE/2-1, ctr-STIMSIZE/2:ctr+STIMSIZE/2-1) = stim;
	imwrite(canvas, fname); 
	%imwrite(medfilt2(canvas), fname);     % MEDIAN FILTERING !!!!!
	t=imread(fname);
    fid = fopen([fname '.raw'], 'w');
    fwrite(fid, t');
    fclose(fid);

	sstim = imresize(z, [STIMSIZE/2 STIMSIZE/2]);
	fname=[odir num2str(i) '.small.' EXTENSION '.png'];	
	canvas = .5*ones(CSIZE); 
	canvas(ctr-STIMSIZE/4:ctr+STIMSIZE/4-1, ctr-STIMSIZE/4:ctr+STIMSIZE/4-1) = sstim;
	imwrite(canvas, fname); 
	%imwrite(medfilt2(canvas), fname);     % MEDIAN FILTERING !!!!!
	t=imread(fname);
    fid = fopen([fname '.raw'], 'w');
    fwrite(fid, t');
    fclose(fid);

	for j=1:NBCLUT

		%Checks if 'natural images' used for background have the right size
		fname = [natdir listnat{ceil(rand() * numel(listnat))}];
		background=imread(fname);
		if ((size(background, 1) ~= CSIZE) || (size(background, 1) ~= CSIZE)) error('Wrong image size!'); end;

		% Contrast maximization
		z = double(background);
		z = z-min(z(:)); z = z ./ max(z(:));
		background=z;

		if (mod(j, 4) == 1) xp = 2; yp=2; end;
		if (mod(j, 4) == 2) xp = 130; yp=2; end;
		if (mod(j, 4) == 3) xp = 2; yp=130; end;
		if (mod(j, 4) == 0) xp = 130; yp=130; end;
		
		fname=[odir num2str(i) '.clut' num2str(j) '.' EXTENSION '.png'];	
		canvas = .5*ones(CSIZE); 
		canvas(yp:yp + STIMSIZE-1, xp:xp + STIMSIZE-1) = stim;
		alphapix=find(abs(canvas - .5)<=.01);
		canvas(alphapix) = background(alphapix);
		imwrite(canvas, fname); 

		t=imread(fname);
		fid = fopen([fname '.raw'], 'w');
		fwrite(fid, t');
		fclose(fid);
	end;


	fposid = fopen([odir 'unpredpos.txt'], 'a');
	for j=1:NBUNPRED

		%fname=listnat{ceil(rand() * numel(listnat))};
		fname=[natdir listnat{ceil(rand() * 250)}]; % Using only half of all natural images, so other half can be used for training on backgrounds alone
		background=imread(fname);
		%Checks if 'natural images' used for background have the right size
		if ((size(background, 1) ~= CSIZE) || (size(background, 1) ~= CSIZE)) error('Wrong image size!'); end;

			%NOTE : THIS DOESN'T DO WHAT I THINK IT DOES - or something goes
			%weird. What actually happens if an image is too large is that only
			%the top-left sectino of the image is used. Which I guess is OK, as
			%long as we're aware of it.


		% Contrast maximization
		z = double(background);
		z = z-min(z(:)); z = z ./ max(z(:));
		background=z;

		%USIZE = round(STIMSIZE * .75 * (1 + rand()));
		%USIZE = round(STIMSIZE * .6666);
		USIZE = STIMSIZE / 2;
		ustim = imresize(stim, [USIZE USIZE]);
		%ustim = imresize(.5 + .5 * (stim - .5), [USIZE USIZE]);

		xp = 2 + floor(rand() * (CSIZE - USIZE - 4));
		yp = 2 + floor(rand() * (CSIZE - USIZE - 4));
		fprintf(fposid, '%d %d %d %d %d\n', i, j, xp, yp, USIZE); 
		
		fname=[odir num2str(i) '.unpred' num2str(j) '.' EXTENSION '.png'];	
		canvas = .5*ones(CSIZE); 
		canvas(yp:yp + USIZE-1, xp:xp + USIZE-1) = ustim;
		alphapix=find(abs(canvas - .5)<=.005);
		canvas(alphapix) = background(alphapix);
		imwrite(canvas, fname); 
		%imwrite(medfilt2(canvas), fname);     % MEDIAN FILTERING !!!!!

		t=imread(fname);
		fid = fopen([fname '.raw'], 'w');
		fwrite(fid, t');
		fclose(fid);

	end;
	fclose(fposid);

	for j=1:NBDOUB

		% Pick a random other object
		%fname = list{ceil(rand() * NBOBJS)};
		fname = [imagedatadir list{j}];
		z=load(fname); 
		z = z.image;
		maxz = max(abs(z(:)));
		z = z ./ (1e-10+maxz);  % Maximizing contrast within the [-1,1] range. Note that the background, which is at 0, is not affected.
		z = .5 + .5 * z;
		z(1:MARGIN,:)=[]; z(end-MARGIN+1:end,:)=[];
		z(:,1:MARGIN)=[]; z(:,end-MARGIN+1:end)=[];
		otherstim = imresize(z, [STIMSIZE STIMSIZE], 'bicubic');

		if (mod(j, 2) == 1) xp1 = 2; yp1 = 2; xp2 = 130; yp2 = 130; end;
		if (mod(j, 2) == 0) xp2 = 2; yp2 = 2; xp1 = 130; yp1 = 130; end;
		
		fname=[odir num2str(i) '.doublet' num2str(j) '.' EXTENSION '.png'];	
		canvas = .5*ones(CSIZE); 
		canvas(yp1:yp1 + STIMSIZE-1, xp1:xp1 + STIMSIZE-1) = stim;
		canvas(yp2:yp2 + STIMSIZE-1, xp2:xp2 + STIMSIZE-1) = otherstim;
		imwrite(canvas, fname); 

		t=imread(fname);
		fid = fopen([fname '.raw'], 'w');
		fwrite(fid, t');
		fclose(fid);
	end;

	farrayid = fopen([odir 'array.txt'], 'a');
	for j=1:NBARRAY
		INSMARG=20; %Insertion margin
		otherstim={};
		RSIZE = floor(CSIZE / 3 - 2) - 2*INSMARG;
		% Pick 9 random other objects - one of which will be obliterated by the current object
		for numo=1:9
			othero=i;
			while(othero == i)
				othero = ceil(rand() * NBOBJS);
				%othero = ceil(NBOBJS+rand() * NBOBJS);
			end;
			fname = [imagedatadir list{othero}];  
			z=load(fname); 
			z = z.image;
			maxz = max(abs(z(:)));
			z = z ./ (1e-10+maxz);  % Maximizing contrast within the [-1,1] range. Note that the background, which is at 0, is not affected.
			z = .5 + .5 * z;
			z(1:MARGIN,:)=[]; z(end-MARGIN+1:end,:)=[];
			z(:,1:MARGIN)=[]; z(:,end-MARGIN+1:end)=[];

			z = imresize(z, [RSIZE RSIZE]);
			otherstim{numo} = z;
		end;

		rstim = imresize(stim, [RSIZE RSIZE]);

		canvas = .5*ones(CSIZE); 
		for pos=1:9
			posc = mod(pos-1, 3)+1;
			posl = floor((pos -1)/ 3)+1;
			%canvas(1+ (posl-1)*RSIZE:posl*RSIZE, 1+ (posc-1)*RSIZE:posc*RSIZE) = otherstim{pos};
%			canvas((posl-1)*(RSIZE+2*INSMARG)+RSIZE:posl*(RSIZE+2*INSMARG)-RSIZE+1,(posc-1)*(RSIZE+2*INSMARG)+RSIZE:posc*(RSIZE+2*INSMARG)-RSIZE+1)  = otherstim{pos};
canvas(1+(posl-1)*(RSIZE+2*INSMARG)+INSMARG:posl*(RSIZE+2*INSMARG)-INSMARG,1+(posc-1)*(RSIZE+2*INSMARG)+INSMARG:posc*(RSIZE+2*INSMARG)-INSMARG) = otherstim{pos};

		end;
		posstim = ceil(rand()* 9);
		posc = mod(posstim-1, 3)+1;
		posl = floor((posstim -1)/ 3)+1;
		%canvas(1+ (posl-1)*RSIZE:posl*RSIZE, 1+ (posc-1)*RSIZE:posc*RSIZE) = rstim;
			canvas(1+(posl-1)*(RSIZE+2*INSMARG)+INSMARG:posl*(RSIZE+2*INSMARG)-INSMARG,1+(posc-1)*(RSIZE+2*INSMARG)+INSMARG:posc*(RSIZE+2*INSMARG)-INSMARG)  = rstim;
		
		fprintf(farrayid, '%d %d %d %d %d\n', i, j, 1+(posc-1)*(RSIZE+2*INSMARG)+INSMARG, 1+(posl-1)*(RSIZE+2*INSMARG)+INSMARG, RSIZE); 
		
		
		fname=[odir num2str(i) '.array' num2str(j) '.' EXTENSION '.png'];	
		imwrite(canvas, fname); 

		t=imread(fname);
		fid = fopen([fname '.raw'], 'w');
		fwrite(fid, t');
		fclose(fid);
	end;
	fclose(farrayid);

	fsizeid = fopen([odir 'quadsize.txt'], 'a');
	for j=1:NBQUADSIZE
		otherstim={};
		% Pick 4 random other objects - one of which will be obluiterated by the current object
		for numo=1:4
			othero=i;
			while(othero == i)
				othero = ceil(rand() * NBOBJS);
				%othero = ceil(NBOBJS+rand() * NBOBJS);
			end;
			fname = [imagedatadir list{othero}];
			z=load(fname); 
			z = z.image;
			maxz = max(abs(z(:)));
			z = z ./ (1e-10+maxz);  % Maximizing contrast within the [-1,1] range. Note that the background, which is at 0, is not affected.
			z = .5 + .5 * z;
			z(1:MARGIN,:)=[]; z(end-MARGIN+1:end,:)=[];
			z(:,1:MARGIN)=[]; z(:,end-MARGIN+1:end)=[];

			%RSIZE = ceil((.5 + .5 * rand()) * STIMSIZE) - 4;
			RSIZE = ceil((.3333 + .3333 * rand()) * STIMSIZE) ;
			t = imresize(z, [RSIZE RSIZE]);
			sp = ceil((STIMSIZE - RSIZE) / 2);
			z = .5*ones(STIMSIZE);
			z(sp:sp+RSIZE-1, sp:sp+RSIZE-1) = t;

			otherstim{numo} = z;
		end;

		RSIZE = ceil((.3333 + .3333 * rand()) * STIMSIZE) ;
		t = imresize(stim, [RSIZE RSIZE]);
		sp = ceil((STIMSIZE - RSIZE) / 2);
		rstim  = .5*ones(STIMSIZE);
		rstim(sp:sp+RSIZE-1, sp:sp+RSIZE-1) = t;

		canvas = .5*ones(CSIZE); 
		canvas(2:2 + STIMSIZE-1, 2:2 + STIMSIZE-1) = otherstim{1};
		canvas(130:130 + STIMSIZE-1, 2:2 + STIMSIZE-1) = otherstim{2};
		canvas(2:2 + STIMSIZE-1, 130:130 + STIMSIZE-1) = otherstim{3};
		canvas(130:130 + STIMSIZE-1, 130:130 + STIMSIZE-1) = otherstim{4};
		
		if (mod(j, 4) == 1) xp1 = 2; yp1 = 2;  end;
		if (mod(j, 4) == 2) xp1 = 130; yp1 = 2;  end;
		if (mod(j, 4) == 3) xp1 = 2; yp1 = 130;  end;
		if (mod(j, 4) == 0) xp1 = 130; yp1 = 130;  end;
		
		fprintf(fsizeid, '%d %d %d %d %d\n', i, j, xp1+sp, yp1+sp, RSIZE); 
		
		canvas(yp1:yp1 + STIMSIZE-1, xp1:xp1 + STIMSIZE-1) = rstim;
		
		fname=[odir num2str(i) '.quadsize' num2str(j) '.' EXTENSION '.png'];	
		imwrite(canvas, fname); 

		t=imread(fname);
		fid = fopen([fname '.raw'], 'w');
		fwrite(fid, t');
		fclose(fid);
	end;
	fclose(fsizeid);
	
	fquadid = fopen([odir 'quad.txt'], 'a');
	for j=1:NBQUAD
		otherstim={};
		% Pick 4 random other objects - one of which will be obluiterated by the current object
		for numo=1:4
			othero=i;
			while(othero == i)
				othero = ceil(rand() * NBOBJS);
				%othero = ceil(NBOBJS+rand() * NBOBJS);
			end;
			fname = [imagedatadir list{othero}];
			z=load(fname); 
			z = z.image;
			maxz = max(abs(z(:)));
			z = z ./ (1e-10+maxz);  % Maximizing contrast within the [-1,1] range. Note that the background, which is at 0, is not affected.
			z = .5 + .5 * z;
			z(1:MARGIN,:)=[]; z(end-MARGIN+1:end,:)=[];
			z(:,1:MARGIN)=[]; z(:,end-MARGIN+1:end)=[];

			RSIZE = ceil(.6667 * STIMSIZE);
			t = imresize(z, [RSIZE RSIZE]);
			sp = ceil((STIMSIZE - RSIZE) / 2);
			z = .5*ones(STIMSIZE);
			z(sp:sp+RSIZE-1, sp:sp+RSIZE-1) = t;

			otherstim{numo} = z;
		end;

		RSIZE = ceil(.6667 * STIMSIZE);
		t = imresize(stim, [RSIZE RSIZE]);
		sp = ceil((STIMSIZE - RSIZE) / 2);
		rstim  = .5*ones(STIMSIZE);
		rstim(sp:sp+RSIZE-1, sp:sp+RSIZE-1) = t;

		canvas = .5*ones(CSIZE); 
		canvas(2:2 + STIMSIZE-1, 2:2 + STIMSIZE-1) = otherstim{1};
		canvas(130:130 + STIMSIZE-1, 2:2 + STIMSIZE-1) = otherstim{2};
		canvas(2:2 + STIMSIZE-1, 130:130 + STIMSIZE-1) = otherstim{3};
		canvas(130:130 + STIMSIZE-1, 130:130 + STIMSIZE-1) = otherstim{4};
		
		if (mod(j, 4) == 1) xp1 = 2; yp1 = 2;  end;
		if (mod(j, 4) == 2) xp1 = 130; yp1 = 2;  end;
		if (mod(j, 4) == 3) xp1 = 2; yp1 = 130;  end;
		if (mod(j, 4) == 0) xp1 = 130; yp1 = 130;  end;
		
		fprintf(fquadid, '%d %d %d %d %d\n', i, j, xp1+sp, yp1+sp, RSIZE); 
		
		canvas(yp1:yp1 + STIMSIZE-1, xp1:xp1 + STIMSIZE-1) = rstim;
		
		fname=[odir num2str(i) '.quad' num2str(j) '.' EXTENSION '.png'];	
		imwrite(canvas, fname); 

		t=imread(fname);
		fid = fopen([fname '.raw'], 'w');
		fwrite(fid, t');
		fclose(fid);
	end;
	fclose(fquadid);
	
	for j=1:NBCONT
		otherstimlc={};
		% Pick 4 random other objects - one of which will be obluiterated by the current object
		for numo=1:4
			fname = [imagedatadir list{ceil(rand() * NBOBJS)}];
			z=load(fname); 
			z = z.image;
			maxz = max(abs(z(:)));
			z = z ./ (1e-10+(maxz*2.0));  % Maximizing contrast within the [-.5,.5] range - i.e. low contrast
			z = .5 + .5 * z;
			z(1:MARGIN,:)=[]; z(end-MARGIN+1:end,:)=[];
			z(:,1:MARGIN)=[]; z(:,end-MARGIN+1:end)=[];
			otherstimlc{numo} = imresize(z, [STIMSIZE STIMSIZE], 'bicubic');
		end;

		canvas = .5*ones(CSIZE); 
		canvas(2:2 + STIMSIZE-1, 2:2 + STIMSIZE-1) = otherstimlc{1};
		canvas(130:130 + STIMSIZE-1, 2:2 + STIMSIZE-1) = otherstimlc{2};
		canvas(2:2 + STIMSIZE-1, 130:130 + STIMSIZE-1) = otherstimlc{3};
		canvas(130:130 + STIMSIZE-1, 130:130 + STIMSIZE-1) = otherstimlc{4};
		
		if (mod(j, 4) == 1) xp1 = 2; yp1 = 2;  end;
		if (mod(j, 4) == 2) xp1 = 130; yp1 = 2;  end;
		if (mod(j, 4) == 3) xp1 = 2; yp1 = 130;  end;
		if (mod(j, 4) == 0) xp1 = 130; yp1 = 130;  end;
		
		canvas(yp1:yp1 + STIMSIZE-1, xp1:xp1 + STIMSIZE-1) = stim;
		
		fname=[odir num2str(i) '.cont' num2str(j) '.' EXTENSION '.png'];	
		imwrite(canvas, fname); 

		t=imread(fname);
		fid = fopen([fname '.raw'], 'w');
		fwrite(fid, t');
		fclose(fid);
	end;
end;


