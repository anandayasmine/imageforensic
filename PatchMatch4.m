function [ debug ] = PatchMatch4(feat, param)
%% Inisialisasi Awal
tic();


% Pencatatan waktu
fprintf('Initalizing... ');
tic();
fprintf(' \n >')

% Ukuran Img sebelum di padding
pjg = size(feat(:,:,1),1);
lbr = size(feat(:,:,1),2);

% Padding Img
targetImg = feat(:,:,:);
targetImg_pad = padarray(targetImg,[param.midPatch-1 param.midPatch-1],0);%padding img
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
param.iter_ini=5;
for i = 1:pjg %row
    for j = 1:lbr %col
        for k=1:param.iter_ini
            NNF(i,j,1) = randi([1 ,pjg],[1,1]); 
            NNF(i,j,2) = randi([1 ,lbr],[1,1]);
            temp = UpdateOffsets4_ini( i, j,param.minNNF, param.midPatch, NNF, offsets, targetImg_pad, param.sizePatch );
            offsets_raw(k,1:3) = [temp NNF(i,j,1) NNF(i,j,2)];
        end
        [~,row_min] = min(offsets_raw(:,1));
        offsets(i,j) = offsets_raw(row_min,1);
        NNF(i,j,1) = offsets_raw(row_min,2);
        NNF(i,j,2) = offsets_raw(row_min,3);
        
%         offsets = UpdateOffsets4( i, j,param.minNNF, param.midPatch, NNF, offsets, targetImg_pad, param.sizePatch );
           
    end
end

debug.offsets_ini=offsets;
debug.NNF_ini=NNF;
toc();

%% Main Iteration
debug.NNF_prp_rs.nnf = struct;
debug.offsets_prp_rs.ofs = struct;
diff_NNF=[];

for iterasi = 1:param.main_iteration    
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
            if(ofs_prp(idx)<param.minOFF)
            
            
            switch idx
                case 2 % propagate dari baris atas
                    NNF(i,j,1) = max(NNF(max(i-1,1),j,1)-1,1);
                    NNF(i,j,2) = NNF(max(i-1,1),j,2);
                    offsets = UpdateOffsets4( i, j, param.minNNF, param.midPatch, NNF, offsets, targetImg_pad, param.sizePatch );
                    
                case 3 % propagate dari kolom kiri
                    NNF(i,j,1) = NNF(i,max(j-1,1),1);
                    NNF(i,j,2) = max(NNF(i,max(j-1,1),2)-1,1);
                    offsets = UpdateOffsets4( i, j, param.minNNF, param.midPatch, NNF, offsets, targetImg_pad, param.sizePatch );
                    
                case 4 % propagate dari diagonal kiri atas
                    NNF(i,j,1) = max(NNF(max(i-1,1),max(j-1,1),1)-1,1);
                    NNF(i,j,2) = max(NNF(max(i-1,1),max(j-1,1),2)-1,1);
                    offsets = UpdateOffsets4( i, j, param.minNNF, param.midPatch, NNF, offsets, targetImg_pad, param.sizePatch );
                     
                case 5 % propagate dari antidiagonal kanan atas
                    NNF(i,j,1) = max(NNF(max(i-1,1),min(j+1,lbr),1)-1,1);
                    NNF(i,j,2) = min(NNF(max(i-1,1),min(j+1,lbr),2)+1,pjg);
                    offsets = UpdateOffsets4( i, j, param.minNNF, param.midPatch, NNF, offsets, targetImg_pad, param.sizePatch );                  
                    
                case 6 % propagate dari baris atas
                    NNF(i,j,1) = max(NNF(max(i-2,1),j,1)-2,1);
                    NNF(i,j,2) = NNF(max(i-2,1),j,2);
                    offsets = UpdateOffsets4( i, j, param.minNNF, param.midPatch, NNF, offsets, targetImg_pad, param.sizePatch );
                    
                case 7 % propagate dari kolom kiri
                    NNF(i,j,1) = NNF(i,max(j-2,1),1);
                    NNF(i,j,2) = max(NNF(i,max(j-2,1),2)-2,1);
                    offsets = UpdateOffsets4( i, j, param.minNNF, param.midPatch, NNF, offsets, targetImg_pad, param.sizePatch );
                    
                case 8 % propagate dari diagonal kiri atas
                    NNF(i,j,1) = max(NNF(max(i-2,1),max(j-2,1),1)-2,1);
                    NNF(i,j,2) = max(NNF(max(i-2,1),max(j-2,1),2)-2,1);
                    offsets = UpdateOffsets4( i, j, param.minNNF, param.midPatch, NNF, offsets, targetImg_pad, param.sizePatch );
                     
                case 9 % propagate dari antidiagonal kanan atas
                    NNF(i,j,1) = max(NNF(max(i-2,1),min(j+1,lbr),1)-2,1);
                    NNF(i,j,2) = min(NNF(max(i-2,1),min(j+1,lbr),2)+2,pjg);
                    offsets = UpdateOffsets4( i, j, param.minNNF, param.midPatch, NNF, offsets, targetImg_pad, param.sizePatch );                  
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
            if(ofs_prp(idx)<param.minOFF)
                
            switch idx
                case 2 % propagate dari baris bawah
                    NNF(i,j,1) = min(NNF(min(i+1,pjg),j,1)+1,pjg);
                    NNF(i,j,2) = NNF(min(i+1,pjg),j,2);
                    offsets = UpdateOffsets4( i, j, param.minNNF, param.midPatch, NNF, offsets, targetImg_pad, param.sizePatch );
                    
                case 3 % propagate dari kolom kanan
                    NNF(i,j,1) = NNF(i,min(j+1,lbr),1);
                    NNF(i,j,2) = min(NNF(i,min(j+1,lbr),2)+1,lbr);
                    offsets = UpdateOffsets4( i, j, param.minNNF, param.midPatch, NNF, offsets, targetImg_pad, param.sizePatch );
                    
                case 4 % propagate dari diagonal kanan bawah
                    NNF(i,j,1) = min(NNF(min(i+1,pjg),min(j+1,lbr),1)+1,pjg);
                    NNF(i,j,2) = min(NNF(min(i+1,pjg),min(j+1,lbr),2)+1,lbr);
                    offsets = UpdateOffsets4( i, j, param.minNNF, param.midPatch, NNF, offsets, targetImg_pad, param.sizePatch );
                     
                case 5 % propagate dari antidiagonal kiri bawah
                    NNF(i,j,1) = min(NNF(min(i+1,pjg),max(j-1,1),1)+1,pjg);
                    NNF(i,j,2) = max(NNF(min(i+1,pjg),max(j-1,1),2)-1,1);
                    offsets = UpdateOffsets4( i, j, param.minNNF, param.midPatch, NNF, offsets, targetImg_pad, param.sizePatch );                  
                    
                case 6 % propagate dari baris bawah
                    NNF(i,j,1) = min(NNF(min(i+2,pjg),j,1)+2,pjg);
                    NNF(i,j,2) = NNF(min(i+2,pjg),j,2);
                    offsets = UpdateOffsets4( i, j, param.minNNF, param.midPatch, NNF, offsets, targetImg_pad, param.sizePatch );
                    
                case 7 % propagate dari kolom kanan
                    NNF(i,j,1) = NNF(i,min(j+2,lbr),1);
                    NNF(i,j,2) = min(NNF(i,min(j+2,lbr),2)+2,lbr);
                    offsets = UpdateOffsets4( i, j, param.minNNF, param.midPatch, NNF, offsets, targetImg_pad, param.sizePatch );
                    
                case 8 % propagate dari diagonal kanan bawah
                    NNF(i,j,1) = min(NNF(min(i+2,pjg),min(j+2,lbr),1)+2,pjg);
                    NNF(i,j,2) = min(NNF(min(i+2,pjg),min(j+2,lbr),2)+2,lbr);
                    offsets = UpdateOffsets4( i, j, param.minNNF, param.midPatch, NNF, offsets, targetImg_pad, param.sizePatch );
                     
                case 9 % propagate dari antidiagonal kiri bawah
                    NNF(i,j,1) = min(NNF(min(i+2,pjg),max(j-2,1),1)+2,pjg);
                    NNF(i,j,2) = max(NNF(min(i+2,pjg),max(j-2,1),2)-2,1);
                    offsets = UpdateOffsets4( i, j, param.minNNF, param.midPatch, NNF, offsets, targetImg_pad, param.sizePatch );
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
            if( norm([iNNF(k),jNNF(k)]-[i,j]) > param.minNNF ) 
               a= targetImg_pad(...
                    i : i+param.sizePatch-1,...
                    j : j+param.sizePatch-1,:);
               b= targetImg_pad(...
                    iNNF(k) : iNNF(k)+param.sizePatch-1,...
                    jNNF(k) : jNNF(k)+param.sizePatch-1,:);
                
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

debug.NNF_prp_rs(iterasi).nnf = NNF;
debug.offsets_prp_rs(iterasi).ofs = offsets;
toc();
fprintf('\n >')
% debug.iterasi = iterasi;

% Kondisi Iterasi
% if (iterasi > 5)
%     diff_NNF = (norm(debug.offsets_prp_rs(iterasi-1).ofs(:,:)-debug.offsets_prp_rs(iterasi).ofs(:,:)))
%     if(  diff_NNF < param.minItrNNF) 
%         debug.diff_NNF = [debug.diff_NNF; diff_NNF];
%         debug.iterasi = iterasi;
%         break;
%     end
% end

debug.NNF=NNF;
debug.NNF_row = NNF(:,:,1);
debug.NNF_col = NNF(:,:,2);
debug.offsets=offsets;

end

