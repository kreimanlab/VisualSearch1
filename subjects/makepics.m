NBOBJS=40;
PSIZE= 320; % 500; %320;  %430;
STIMSIZE= 44; %48 with diam 11; % 56;  % 88 120; % 44; %88 %120; smallest is best!
ctr=PSIZE/2;
ARRAYSIZE=6;
pic=128*ones(PSIZE);

%ECC = params.ECCENTRICITY;
ECC=  110;  % 180; %110;  %150

datadir = './subjectsdata/';

%Note : These files need to be the same as in attentionMapC
z=load([datadir 'ResultsSaccades_AndrewT_blk1.mat']);
z2=load([datadir 'ResultsSaccades_AlissaT_blk2.mat']);
z3=load([datadir 'ResultsSaccades_Fetemeh_blk3.mat']);
z4=load([datadir 'ResultsSaccades_Mario_blk4.mat']);
z5=load([datadir 'ResultsSaccades_Vivian_blk5.mat']);
z6=load([datadir 'ResultsSaccades_Sam_blk6.mat']);
z7=load([datadir 'ResultsSaccades_Claudio_blk7.mat']);
z8=load([datadir 'ResultsSaccades_Alexandru_blk8.mat']);
z9=load([datadir 'ResultsSaccades_Mona_blk9.mat']);
z10=load([datadir 'ResultsSaccades_HC_blk10.mat']);
results=[z.results z2.results z3.results z4.results z5.results z6.results z7.results z8.results z9.results z10.results]; params = z.params; 
sequence=[z.sequence z2.sequence z3.sequence z4.sequence z5.sequence z6.sequence z7.sequence z8.sequence z9.sequence z10.sequence];
BLOCKSIZE = numel(z.sequence);
oshown=[];

for numpic=1:ARRAYSIZE
		 angle = pi/2 + (.5 + numpic-1) * 2*pi / ARRAYSIZE;
		 myx = ctr + round(cos(angle) * ECC);
		 myy = ctr + round(sin(angle) * ECC);
		 mrect=[myy-STIMSIZE/2, myy+STIMSIZE/2-1, myx-STIMSIZE/2, myx+STIMSIZE/2-1];
		 destrect{numpic} = mrect;
end;

save('destrect.mat', 'destrect', 'PSIZE');

img={}; textures={}; imgbig={};
for i=1:NBOBJS
	if params.SHINED == 1
		img{i} = imresize(imread(['./objectimages/SHINEd_' num2str(i) '.normal.png']), [STIMSIZE, STIMSIZE], 'nearest') ;
		imgbig{i} = imresize(imread(['./objectimages/SHINEd_' num2str(i) '.normal.png']), [2*STIMSIZE, 2*STIMSIZE], 'nearest') ;
	else
		img{i} = imresize(imread(['./objectimages/' num2str(i) '.normal.png']), [STIMSIZE, STIMSIZE], 'nearest');
		imgbig{i} = imresize(imread(['./objectimages/' num2str(i) '.normal.png']), [2*STIMSIZE, 2*STIMSIZE], 'nearest');
	end;
	isolatedimg=128*ones(PSIZE); ctri=size(isolatedimg, 1)/2; isolatedimg(ctri-STIMSIZE/2:ctri+STIMSIZE/2-1, ctri-STIMSIZE/2:ctri+STIMSIZE/2-1) = img{i};
	isolatedimgBig=128*ones(PSIZE); ctri=size(isolatedimgBig, 1)/2; isolatedimgBig(ctri-STIMSIZE:ctri+STIMSIZE-1, ctri-STIMSIZE:ctri+STIMSIZE-1) = imgbig{i};
	fname = ['./images/obj' num2str(i) '.png'];
	imwrite(isolatedimg./256, fname);
	t=imread(fname);
    fid = fopen([fname '.raw'], 'w');
    fwrite(fid, t');
    fclose(fid);
	%fname = ['./images/obj' num2str(i) '.big.png'];
	%imwrite(isolatedimgBig./256, fname);
	%t=imread(fname);
    %fid = fopen([fname '.raw'], 'w');
    %fwrite(fid, t');
    %fclose(fid);
end;

for numtrial= 1:length(results) 

	disp(numtrial);
	for numpic=1:ARRAYSIZE
			 pic(destrect{numpic}(1):destrect{numpic}(2), destrect{numpic}(3):destrect{numpic}(4)) = img{results{numtrial}.probes(numpic)};
	end;
	fname =  ['./images/arraySub_blk' num2str(1+floor((numtrial-1) / BLOCKSIZE)) '_num' num2str(sequence(numtrial)) '.png'];
	imwrite(pic./256, fname);
	t=imread(fname);
    fid = fopen([fname '.raw'], 'w');
    fwrite(fid, t');
    fclose(fid);
end;








