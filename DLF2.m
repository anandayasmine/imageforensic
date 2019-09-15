function [ mask, outData ] = DLF2( mpf_y,mpf_x,param )
% dlf_fig = figure();
% regularize offsets field by median filtering
[DD_med, NN_med] = genDisk(param.rd_median); 
NN_med = (NN_med+1)/2;
[ mpf_y,  mpf_x] = MPFregularize(mpf_y,mpf_x,DD_med,NN_med);
% subplot(2,4,1); imshow(uint8(mpf_y));

% Compute the squared error of dense linear fitting
DLFerr  =  DLFerror(mpf_y,mpf_x,param.rd_dlf); % ini yg penting
mask    = (DLFerr<=param.th2_dlf);
outData.maskDLF =  mask;
% subplot(2,4,2); imshow((mask));

% removal of close couples
% dist2 = MPFspacedist2(mpf_y,mpf_x);
% mask  = mask & (dist2>=param.th2_dist2);
% subplot(2,4,3); imshow((mask));

% morphological operations
% mask  = bwareaopen(mask,param.th_sizeA,1); 
se1 = strel('disk',4,4);
mask = imclose(mask,se1);
mask = imerode(mask,se1);


mask  = imdilate(mask,se1);
% subplot(2,4,4); imshow((mask));
outData.maskMPF = mask;
mask  = bwareaopen(mask,param.th_sizeA); 

% mask  = MPFdual(mpf_y,mpf_x,mask); % mirroring of detected regions
mask  = bwareaopen(mask,param.th_sizeB,8);
% subplot(2,4,5); imshow((mask));

mask  = imdilate(mask,strel('disk',param.rd_dil));
% subplot(2,4,6); imshow((mask));

% % final threshold
% mask  = bwareaopen(mask,param.th_sizeC,8);



% put the borders
raggioU = ceil((param.diameter-1)/2);
raggioL = floor((param.diameter-1)/2);
mask = padarray_both(mask,[raggioU,raggioU,raggioL,raggioL],false()); 
DLFerr = padarray_both(DLFerr,[raggioU,raggioU,raggioL,raggioL],0); 
% subplot(2,4,7); imshow((mask));

% miroring
% figure, imshow(mask);
% for p = 1:size(mpf_y,1)
%     for q = 1:size(mpf_y,2)
%         if(mask(p,q) == 1)
%             coup_y = mpf_y(p,q);
%             coup_x = mpf_x(p,q);
%             mask(coup_y,coup_x)=1;
%         end
%     end
% end

outData.cnn_end = cat(3,mpf_x,mpf_y);
outData.DLFerr  = DLFerr;
% outData.dlf_fig = dlf_fig; 

%% Remake post processing
% mask_raw = double(repmat(outData.maskDLF,[1,1,1]));
% mask_tidy = bwareaopen(mask_raw,200);
% se1 = strel('line',20,20);
% mask_tidy = imclose(mask_tidy,se1);
% se2 = strel('line',10,10);
% mask_tidy = imerode(mask_tidy,se2);
% se3 = strel('line',10,10);
% mask_tidy = imdilate(mask_tidy,se3);
% mask_tidy = bwareaopen(mask_tidy,200);
% mask = mask_tidy;

end

