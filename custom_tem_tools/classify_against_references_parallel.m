%% description
%% startup
close all, clear all, clc
dalpha = 10;

cd('/Users/jonasfunke/MATLAB/MY_MATLAB_TOOLS/custom_tem_tools/private')
%% load stack of particles
%[fname, pname] = uigetfile('*.img', 'Select the stack of images to be classified');
%img = ReadImagic([pname, fname(1:end-4)]); % read imagic file
floc = '/Users/jonasfunke/Dropbox/Nucleosome_Dropbox/Figure4/TEM/Single-particle angle measurement/No-int/Fob25/150531_PK_GR_13__B2_FRETst_noint_stack.img';
img = ReadImagic(floc(1:end-4)); % read imagic file


%% load reference classes
%[fname_ref, pname_ref] = uigetfile('*.img', 'Select the stack of reference classes.');
%ref = ReadImagic([pname_ref, fname_ref(1:end-4)]); % read imagic file

f_sapcer = '/Users/jonasfunke/Dropbox/Nucleosome_Dropbox/Figure4/TEM/Angle_references/Class_averages_selection/FS_CA_200x200.img';
f_no = '/Users/jonasfunke/Dropbox/Nucleosome_Dropbox/Figure4/TEM/Angle_references/no_interaction/iter940/FS_no-interaction_200x200.img';
tmp1 = ReadImagic(f_sapcer(1:end-4)); % read imagic file
tmp2 = ReadImagic(f_no(1:end-4)); % read imagic file

ref = cat(3, tmp1, tmp2);

N_img = size(img,3);
N_ref = size(ref,3);
N = size(img,1);
%% filter references

r_filter = 20;
f_filter = fspecial('gaussian', r_filter , r_filter); % gaussian filter
ref_filtered = zeros(size(ref,1),size(ref,2),N_ref,'single');
subt  = zeros(size(ref,1),size(ref,2),1,N_ref); %, 'single');
for r=1:N_ref
    ref_filtered(:,:,r) = double(ref(:,:,r))-double(imfilter(ref(:,:,r), f_filter, 'same'));
    subt(:,:,1,r) = double(imfilter(ref(:,:,r), f_filter, 'same'));
end

% plot references
bla = zeros(size(ref,1),size(ref,2),1,N_ref, 'single');
for r=1:N_ref
    bla(:,:,1,r) = ref(:,:,r);
end
%figure(1)
%montage(bla, 'DisplayRange', [min(bla(:)) max(bla(:))])

figure(1)
%montage(ref_filtered, 'DisplayRange', [min(ref_filtered(:)) max(ref_filtered(:))])
ref_filt_plot = reshape(ref_filtered, 200,200,1,N_ref);
montage(ref_filt_plot, 'DisplayRange', [min(ref_filt_plot(:)) max(ref_filt_plot(:))])

%% Classify

% start cluster/parallel pool
my_pool = parpool(2);


mycluster = parcluster('SharedCluster');
n_worker=63;
j = batch(mycluster, @classify_images, {img, ref_filtered, dalpha }, 1,'CaptureDiary',true, 'CurrentDirectory', '.','AdditionalPaths', {'/nfs/matlabuser/jonasfunke/MATLAB/parallel_test'}, 'Pool',  n_worker(n));
wait(j)
result = j.fetchOutputs{1};% Get results into a cell array

% close parallel pool
%% write output

