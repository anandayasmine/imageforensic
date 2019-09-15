function [ output ] = UpdateOffsets4_ini( i, j, minNNF, midPatch, NNF, offsets, targetImg_pad, sizePatch )
% offset yang NNF nya telah di update Akan di match dengan kordinat NNF baru
    iNNF = (NNF(i,j,1));
    jNNF = (NNF(i,j,2));

% jarak patch ini harus > TD1
if( norm([iNNF,jNNF]-[i,j]) > minNNF )   
    a=  targetImg_pad(...
        i : i+sizePatch-1,...
        j : j+sizePatch-1,:);
    b= targetImg_pad(...
        iNNF : iNNF+sizePatch-1,...
        jNNF : jNNF+sizePatch-1,:);
    
    % Mean Square Error
    temp1 = a-b;
    temp2  = temp1(~isnan(temp1(:)));    
    temp3 = sum(temp2.^2)/length(temp2);
%     offsets(i,j)=temp3;
    
    % Normalisasi Euclidean Distance
%     diff_R=norm(a(:,:,1)-b(:,:,1));
%     diff_G=norm(a(:,:,2)-b(:,:,2));
%     diff_B=norm(a(:,:,2)-b(:,:,2));
%     offsets(i,j)=((diff_R+diff_G+diff_B)/3);

else
    temp3 = 1000;
%     offsets(i,j)=temp3;
end

% output = offsets;  
output = temp3;  
end

