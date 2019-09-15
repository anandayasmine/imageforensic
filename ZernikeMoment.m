function [ feat ] = ZernikeMoment( I,param )

if isscalar(param.ORDER),
    NM = zeros(0,2); 
    for n = 0:param.ORDER, 
        for m = 0:n, 
            if mod(n-abs(m),2)==0,
                NM = cat(1, NM , [n m]);
            end 
        end
    end
else
    NM  = param.ORDER;
end;
len   = size(NM,1);
param.ORDER = max(NM(:,1));
F = factorial(0:param.ORDER);

coef = zeros(len,param.ORDER);
WF   = zeros(len,1);
for index = 1:len
    n = NM(index,1);
    m = NM(index,2);
    tu = floor((n+abs(m))/2);
    td = floor((n-abs(m))/2);
    for s = 0:td, 
       coef(index,s+1) = ( ((-1).^s) .* F(1+n-s) )./ (F(1+s).*F(1+tu-s).*F(1+td-s)); 
    end; 
    WF(index) = sqrt((n+1));
end;

gf = @(x,y) max(0,1-abs(x)).*max(0,1-abs(y)); % bi-linear

[X,Y] = meshgrid(1:param.SZ,1:param.SZ);
Y = (2.*Y-param.SZ-1)./2;
X = (2.*X-param.SZ-1)./2;

radiusMax = (param.SZ-1)/2;
radiusMin = 0;
radiusPoints = linspace(radiusMin,radiusMax,param.radiusNum);

BF = zeros(param.SZ,param.SZ,len);

for indexA = 1:param.anglesNum;
    A  = (indexA-1)*360./param.anglesNum;

    for indexR = 1:param.radiusNum;
        R = radiusPoints(indexR);
        pX = X - R*cosd(A);
        pY = Y - R*sind(A);
        J = gf(pX,pY);
        R = 2*R/param.SZ;
        for index = 1:len
            n = NM(index,1);
            m = NM(index,2);
            Rad = 0;
            td = floor((n-abs(m))/2);
            for s = 0:td, 
                Rad = Rad + coef(index,s+1).*(R.^(n-2*s));
            end;
            Rad = Rad.*(cosd(m*A) + 1i*sind(m*A));
            BF(:,:,index) = BF(:,:,index)+(J*Rad)*R;
        end;
    end;
end;

WF = WF./sqrt(sum(sum(abs(BF(:,:,1)).^2)));

bfdata                = struct;
bfdata.number         = size(NM,1);
bfdata.param.ORDERs   = NM;
bfdata.bf             = BF;
bfdata.factor         = WF;

% Filtering
len = bfdata.number;
y   = zeros(cat(2,size(I),len));
for index = 1:len
    ac=imfilter(I, conj(bfdata.bf(:,:,index)));
    y(:,:,:,index) = (bfdata.factor(index)) .* ac;
end
feat = abs(y);

end

