function [ debug, mask_pp] = PM_DLF_GUI_OTO(feat, param, handles,I)
%% Inisialisasi Awal
tic();

% Pencatatan waktu
fprintf('Initalizing... ');
tic();
fprintf(' \n >')

% Ukuran Img sebelum di padding
pjg = size(feat(:,:,1),1);
lbr = size(feat(:,:,1),2);
NNF=[];
% Padding Img
targetImg = feat(:,:,:);
targetImg_pad = padarray(targetImg,[param.M.midPatch-1 param.M.midPatch-1],0);%padding img
% pjg_pad = size(targetImg_pad,1);
% lbr_pad = size(targetImg_pad,2);
offsets = inf(pjg,lbr);

% Pembuatan NNF secara random
% if (lbr == 320)
%     load NNF_I;
% elseif (lbr == 640)
%     load NNF_II;
% elseif (lbr == 128)
%     load NNF_128;
% end

% Inisialisasi offsets
offsets_raw=[];
param.M.iter_ini=5;
for i = 1:pjg %row
    for j = 1:lbr %col
        for k=1:param.M.iter_ini
            NNF(i,j,1) = randi([1 ,pjg],[1,1]); 
            NNF(i,j,2) = randi([1 ,lbr],[1,1]);
            temp = UpdateOffsets4_ini( i, j,param.M.minNNF, param.M.midPatch, NNF, offsets, targetImg_pad, param.M.sizePatch );
            offsets_raw(k,1:3) = [temp NNF(i,j,1) NNF(i,j,2)];
        end
        [~,row_min] = min(offsets_raw(:,1));
        offsets(i,j) = offsets_raw(row_min,1);
        NNF(i,j,1) = offsets_raw(row_min,2);
        NNF(i,j,2) = offsets_raw(row_min,3);
    end
end

debug=struct;
debug.offsets_ini=offsets;
debug.NNF_ini=NNF;
debug.NNF=[];

mask_pp=struct;
mask_pp.pp_dlf=struct;
mask_pp.mask=struct;
debug.NNF_prp_rs.nnf = struct;
debug.offsets_prp_rs.ofs = struct;
diff_NNF=[];

toc();

%% Main Iteration

for iterasi = 1:param.M.main_iteration    
ganjil = mod(iterasi,2)==1;

if ganjil
    disp(['Iterasi ke - ',num2str(iterasi),' (raster scan order)']);
    ii = 1:pjg;
    jj = 1:lbr;
else
    disp(['Iterasi ke - ',num2str(iterasi),' (reverse raster scan order)']);
    ii = pjg:-1:1;
    jj = lbr:-1:1;
end

% Pencatatan waktu
tic();  
% Sub Iteration
for i = ii
    for j = jj
        %% Propagation
        
        if ganjil % propagate dari kiri atas
            
            % tengah - baris atas - kolom kiri - diagonal kiri atas - antidiagonal kanan atas
            
            % Zero Order
            ofs_prp(1) = offsets(i,j); % tengah
            ofs_prp(2) = offsets(max(1,i-1),j); % 1 baris atas
            ofs_prp(3) = offsets(i,max(1,j-1)); % 1 kolom kiri
            ofs_prp(4) = offsets(max(1,i-1),max(1,j-1)); % 1 diagonal
            ofs_prp(5) = offsets(max(1,i-1),min(lbr,j+1)); % 1 antidiagonal
            
            % First Order
            ofs_prp(6) = offsets(max(1,i-2),j); % 2 baris atas
            ofs_prp(7) = offsets(i,max(1,j-2)); % 2 kolom kiri
            ofs_prp(8) = offsets(max(1,i-2),max(1,j-2)); % 2 diagonal
            ofs_prp(9) = offsets(max(1,i-2),min(lbr,j+2)); % 2 antidiagonal
            
            % Cari offset dengan nilai terkecil dari ke 9 kandidat ofs tsb
            [~,idx] = min(ofs_prp);
            if(ofs_prp(idx)<param.M.minOFF)
            
            
            switch idx
                case 2 % propagate dari baris atas
                    NNF(i,j,1) = max(NNF(max(i-1,1),j,1)-1,1);
                    NNF(i,j,2) = NNF(max(i-1,1),j,2);
                    offsets = UpdateOffsets4( i, j, param.M.minNNF, param.M.midPatch, NNF, offsets, targetImg_pad, param.M.sizePatch );
                    
                case 3 % propagate dari kolom kiri
                    NNF(i,j,1) = NNF(i,max(j-1,1),1);
                    NNF(i,j,2) = max(NNF(i,max(j-1,1),2)-1,1);
                    offsets = UpdateOffsets4( i, j, param.M.minNNF, param.M.midPatch, NNF, offsets, targetImg_pad, param.M.sizePatch );
                    
                case 4 % propagate dari diagonal kiri atas
                    NNF(i,j,1) = max(NNF(max(i-1,1),max(j-1,1),1)-1,1);
                    NNF(i,j,2) = max(NNF(max(i-1,1),max(j-1,1),2)-1,1);
                    offsets = UpdateOffsets4( i, j, param.M.minNNF, param.M.midPatch, NNF, offsets, targetImg_pad, param.M.sizePatch );
                     
                case 5 % propagate dari antidiagonal kanan atas
                    NNF(i,j,1) = max(NNF(max(i-1,1),min(j+1,lbr),1)-1,1);
                    NNF(i,j,2) = min(NNF(max(i-1,1),min(j+1,lbr),2)+1,pjg);
                    offsets = UpdateOffsets4( i, j, param.M.minNNF, param.M.midPatch, NNF, offsets, targetImg_pad, param.M.sizePatch );                  
                    
                case 6 % propagate dari baris atas
                    NNF(i,j,1) = max(NNF(max(i-2,1),j,1)-2,1);
                    NNF(i,j,2) = NNF(max(i-2,1),j,2);
                    offsets = UpdateOffsets4( i, j, param.M.minNNF, param.M.midPatch, NNF, offsets, targetImg_pad, param.M.sizePatch );
                    
                case 7 % propagate dari kolom kiri
                    NNF(i,j,1) = NNF(i,max(j-2,1),1);
                    NNF(i,j,2) = max(NNF(i,max(j-2,1),2)-2,1);
                    offsets = UpdateOffsets4( i, j, param.M.minNNF, param.M.midPatch, NNF, offsets, targetImg_pad, param.M.sizePatch );
                    
                case 8 % propagate dari diagonal kiri atas
                    NNF(i,j,1) = max(NNF(max(i-2,1),max(j-2,1),1)-2,1);
                    NNF(i,j,2) = max(NNF(max(i-2,1),max(j-2,1),2)-2,1);
                    offsets = UpdateOffsets4( i, j, param.M.minNNF, param.M.midPatch, NNF, offsets, targetImg_pad, param.M.sizePatch );
                     
                case 9 % propagate dari antidiagonal kanan atas
                    NNF(i,j,1) = max(NNF(max(i-2,1),min(j+1,lbr),1)-2,1);
                    NNF(i,j,2) = min(NNF(max(i-2,1),min(j+1,lbr),2)+2,pjg);
                    offsets = UpdateOffsets4( i, j, param.M.minNNF, param.M.midPatch, NNF, offsets, targetImg_pad, param.M.sizePatch );                  
            end    
            end
        else % propagate dari kanan bawah
            
            % tengah - baris bawah - kolom kanan - diagonal kanan bawah - antidiagonal kiri bawah
            
            % Zero Order
            ofs_prp(1) = offsets(i,j); % tengah
            ofs_prp(2) = offsets(min(i+1,pjg),j); % 1 baris bawah
            ofs_prp(3) = offsets(i,min(j+1,lbr)); % 1 kolom kanan
            ofs_prp(4) = offsets(min(i+1,pjg),min(j+1,lbr)); % 1 diagonal
            ofs_prp(5) = offsets(max(i-1,1),min(j+1,lbr)); % 1 antidiagonal
            
            % First Order
            ofs_prp(6) = offsets(min(i+2,pjg),j); % 2 baris bawah
            ofs_prp(7) = offsets(i,min(j+2,lbr)); % 2 kolom kanan
            ofs_prp(8) = offsets(min(i+2,pjg),min(j+2,lbr)); % 2 diagonal
            ofs_prp(9) = offsets(max(i-1,1),min(j+2,lbr)); % 2 antidiagonal
            
            % Cari offset dengan nilai terkecil dari ke 9 kandidat ofs tsb
            [~,idx] = min(ofs_prp);
            if(ofs_prp(idx)<param.M.minOFF)
                
            switch idx
                case 2 % propagate dari baris bawah
                    NNF(i,j,1) = min(NNF(min(i+1,pjg),j,1)+1,pjg);
                    NNF(i,j,2) = NNF(min(i+1,pjg),j,2);
                    offsets = UpdateOffsets4( i, j, param.M.minNNF, param.M.midPatch, NNF, offsets, targetImg_pad, param.M.sizePatch );
                    
                case 3 % propagate dari kolom kanan
                    NNF(i,j,1) = NNF(i,min(j+1,lbr),1);
                    NNF(i,j,2) = min(NNF(i,min(j+1,lbr),2)+1,lbr);
                    offsets = UpdateOffsets4( i, j, param.M.minNNF, param.M.midPatch, NNF, offsets, targetImg_pad, param.M.sizePatch );
                    
                case 4 % propagate dari diagonal kanan bawah
                    NNF(i,j,1) = min(NNF(min(i+1,pjg),min(j+1,lbr),1)+1,pjg);
                    NNF(i,j,2) = min(NNF(min(i+1,pjg),min(j+1,lbr),2)+1,lbr);
                    offsets = UpdateOffsets4( i, j, param.M.minNNF, param.M.midPatch, NNF, offsets, targetImg_pad, param.M.sizePatch );
                     
                case 5 % propagate dari antidiagonal kiri bawah
                    NNF(i,j,1) = min(NNF(min(i+1,pjg),max(j-1,1),1)+1,pjg);
                    NNF(i,j,2) = max(NNF(min(i+1,pjg),max(j-1,1),2)-1,1);
                    offsets = UpdateOffsets4( i, j, param.M.minNNF, param.M.midPatch, NNF, offsets, targetImg_pad, param.M.sizePatch );                  
                    
                case 6 % propagate dari baris bawah
                    NNF(i,j,1) = min(NNF(min(i+2,pjg),j,1)+2,pjg);
                    NNF(i,j,2) = NNF(min(i+2,pjg),j,2);
                    offsets = UpdateOffsets4( i, j, param.M.minNNF, param.M.midPatch, NNF, offsets, targetImg_pad, param.M.sizePatch );
                    
                case 7 % propagate dari kolom kanan
                    NNF(i,j,1) = NNF(i,min(j+2,lbr),1);
                    NNF(i,j,2) = min(NNF(i,min(j+2,lbr),2)+2,lbr);
                    offsets = UpdateOffsets4( i, j, param.M.minNNF, param.M.midPatch, NNF, offsets, targetImg_pad, param.M.sizePatch );
                    
                case 8 % propagate dari diagonal kanan bawah
                    NNF(i,j,1) = min(NNF(min(i+2,pjg),min(j+2,lbr),1)+2,pjg);
                    NNF(i,j,2) = min(NNF(min(i+2,pjg),min(j+2,lbr),2)+2,lbr);
                    offsets = UpdateOffsets4( i, j, param.M.minNNF, param.M.midPatch, NNF, offsets, targetImg_pad, param.M.sizePatch );
                     
                case 9 % propagate dari antidiagonal kiri bawah
                    NNF(i,j,1) = min(NNF(min(i+2,pjg),max(j-2,1),1)+2,pjg);
                    NNF(i,j,2) = max(NNF(min(i+2,pjg),max(j-2,1),2)-2,1);
                    offsets = UpdateOffsets4( i, j, param.M.minNNF, param.M.midPatch, NNF, offsets, targetImg_pad, param.M.sizePatch );
            end 
            end
        end
        
        
        %% Random Search
        
        % Param Random Search
        radius = size(feat,1)/((1/4)*(1/iterasi)); 
        Radius = round(radius*0.5.^(0:(-floor(log(radius)/log(0.5)))));        
        
        % Membuat radius 
        i_min = max(1,NNF(i,j,1)-Radius(:));
        i_max = min(NNF(i,j,1)+Radius(:),pjg);
        j_min = max(1,NNF(i,j,2)-Radius(:));
        j_max = min(NNF(i,j,2)+Radius(:),lbr);

        % Cari kordinat acak di dalam radius
        iNNF = floor(rand(length(Radius),1).*(i_max(:)-i_min(:)+1)) + i_min(:);
        jNNF = floor(rand(length(Radius),1).*(j_max(:)-j_min(:)+1)) + j_min(:);

        nns(:,1) = NNF(i,j,:);
        nns(:,2:length(Radius)+1) = [iNNF';jNNF'];
        
        ofs_rs(1) = offsets(i,j);
        
        for k = 1:length(Radius)
            % jarak patch ini harus > TD1
            if( norm([iNNF(k),jNNF(k)]-[i,j]) > param.M.minNNF ) 
               a= targetImg_pad(...
                    i : i+param.M.sizePatch-1,...
                    j : j+param.M.sizePatch-1,:);
               b= targetImg_pad(...
                    iNNF(k) : iNNF(k)+param.M.sizePatch-1,...
                    jNNF(k) : jNNF(k)+param.M.sizePatch-1,:);
                
               % Mean Square Error
               temp1 = a-b;
               temp = temp1(~isnan(temp1(:)));
               temp3 = sum(temp.^2)/length(temp);
               ofs_rs(k+1)=temp3;
               
               % Normalisasi Euclidean Distance
%                 diff_R=norm(a(:,:,1)-b(:,:,1));
%                 diff_G=norm(a(:,:,2)-b(:,:,2));
%                 diff_B=norm(a(:,:,2)-b(:,:,2));
%                 offsets(k+1)=((diff_R+diff_G+diff_B)/3);

             else
                temp3 = 10;
                ofs_rs(k+1)=temp3;
            end
        end
        
        % offsets diganti jadi yang paling minimum dari kandidat random search 
        [~,idx] = min(ofs_rs);
        offsets(i,j) = ofs_rs(idx); 
        NNF(i,j,:) = nns(:,idx);
    end
end

axes(handles.axes2);imshow(uint8(NNF(:,:,1)));
axes(handles.axes4);imshow(uint8(NNF(:,:,2)));drawnow 

%% Post Processing: DLF
mpf_y = double((NNF(:,:,2,1)));
mpf_x = double((NNF(:,:,1,1)));
[mask, pp_DLF] = DLF2(mpf_y,mpf_x,param.PP);

mask_pp.mpfy(iterasi).row= mpf_y ;
mask_pp.mpfx(iterasi).col= mpf_x ;
mask_pp.mask(iterasi).pp_dlf = mask;
mask_pp.pp_dlf(iterasi).pp_dlf = pp_DLF;

%% Display offsets
mask2=imresize(mask,[size(feat,1) size(feat,2)]);
mat_off=[];r=1;

for p = 1:size(feat,1)
    for q = 1:size(feat,2)
        if(mask2(p,q) == 1)
            mask_off(p,q) = offsets(p,q);
            mat_off(r,1)  = p;
            mat_off(r,2)  = q;
            mat_off(r,3)  = NNF(p,q,1);
            mat_off(r,4)  = NNF(p,q,2);
            mat_off(r,5)  = offsets(p,q);
            r=r+1;
        end
    end
end
param_off=[];
if(max(max(mask2))~=0)
    mat_off=sortrows(mat_off,5);
    param_off = randi([1,round(size(mat_off,1)*0.1)], [1,20]);
    
    for l=1:1000
        mask(mat_off(l,3),mat_off(l,4))=1;
    end
    mask  = bwareaopen(mask,param.PP.th_sizeA); 
    se1 = strel('disk',4,4);
    mask = imclose(mask,se1);
    mask = imerode(mask,se1);
end   
mat_off_all=struct;
debug.mat_off_all(iterasi).disp_ofs = mat_off;
debug.mask(iterasi).mask=mask;
debug.mask2(iterasi).mask2=mask2;
debug.param_off(iterasi).param_off=param_off;
%% Display Process
mpf_y_pre = double(NNF(:,:,2));
mpf_x_pre = double(NNF(:,:,1));
mpf_y     = double(pp_DLF.cnn_end(:,:,2));
mpf_x     = double(pp_DLF.cnn_end(:,:,1));

dist2 = MPFspacedist2(mpf_y_pre,mpf_x_pre);
max_dist = max(sqrt(dist2(:)));
axes(handles.axes5);imshow(sqrt(dist2),[0, max_dist]); colormap(jet());drawnow 

dist2 = MPFspacedist2(mpf_y,mpf_x);
axes(handles.axes6);imshow(sqrt(dist2),[0, max_dist]); colormap(jet());drawnow  %colorbar();

DLF_db = 10*log10(pp_DLF.DLFerr); DLF_db(DLF_db<-50) = -50;
axes(handles.axes7);imshow(DLF_db,[]); colormap(jet());drawnow 

% === offset ===
axes(handles.axes1);imshow(uint8(I),[]);
if(max(max(mask2))~=0)
    hold on
    for i =param_off
       plot(mat_off(i,2), mat_off(i,1), 'ro');
       plot(mat_off(i,4), mat_off(i,3), 'ro');

       Y=[mat_off(i,1) mat_off(i,3)];
       X=[mat_off(i,2) mat_off(i,4)];
       line(X,Y,'Color','y');
    end
end

axes(handles.axes9);imshow(double(repmat(pp_DLF.maskDLF,[1,1,3])));drawnow 
axes(handles.axes8);imshow(double(repmat(pp_DLF.maskMPF,[1,1,3])));drawnow 
axes(handles.axes10);imshow(double(repmat(mask,[1,1,3])));drawnow 
set(handles.slider1, 'Value', iterasi);drawnow 
set(handles.text19, 'Value', iterasi);drawnow
set(handles.text19, 'String', iterasi);drawnow 

debug.NNF_prp_rs(iterasi).nnf = NNF;
debug.offsets_prp_rs(iterasi).ofs = offsets;
toc();
fprintf('\n >')

if (iterasi>param.M.iter_auto_min )
   diff_mask = norm(mask_pp.mask(iterasi).pp_dlf - mask_pp.mask(iterasi-1).pp_dlf)
   if (diff_mask < param.M.difer_mask)
      break; 
   end
end

end % end of iterasi
set(handles.text16,'String',iterasi);
set(handles.slider1, 'Max', iterasi);drawnow 

debug.NNF=NNF;
debug.NNF_row = NNF(:,:,1);
debug.NNF_col = NNF(:,:,2);
debug.offsets=offsets;

set(handles.text18,'Enable','on');
if(sum(sum(mask)) > 0)
    set(handles.text18,'String','Terdapat Area Duplikasi');
else
    set(handles.text18,'String','Tidak Terdapat Area Duplikasi');
end

end

