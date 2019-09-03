parpool
%% load data set
dir = 'L:\users\Aaron\151124_BMWR30\DenseRun1binaryImages';
conf = loadjson([ dir,  '\conf.json']);
%%
files = rdir([dir '\*.bin']);
%%
x = conf.dims(1);
y = conf.dims(2);
z = conf.dims(3);
n = length(files);
%%
imgs = zeros([x, y, z, n], 'uint16');
%%
tic
parfor i = 1:length(files)
    disp(i)
	fid = fopen(files(i).name,'r');
	imgs(:,:,:,i) = reshape(fread(fid, 'uint16'), x, y, z);
	fclose(fid);
end
toc
