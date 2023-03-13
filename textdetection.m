
colorImage = imread('113.jpg');
I = im2gray(colorImage);

% MSER bölgelerini tespit etmek
[mserRegions, mserConnComp] = detectMSERFeatures(I, ... 
    'RegionAreaRange',[200 8000],'ThresholdDelta',4);

figure
imshow(I)
hold on
plot(mserRegions, 'showPixelList', true,'showEllipses',false)
title('MSER regions')
hold off
% MSER ÖZELLİKLERİNİ Ölçmek için regionprops kullanıyoruz
mserStats = regionprops(mserConnComp, 'BoundingBox', 'Eccentricity', ...
    'Solidity', 'Extent', 'Euler', 'Image');

% sınırlayıcı box verilerini kullanarak en ve boy oranını hesaplıyoruz
bbox = vertcat(mserStats.BoundingBox);
w = bbox(:,3);
h = bbox(:,4);
aspectRatio = w./h;

% hangi bölgelerinin kaldırılacapını belirlemek için verileri eşikliyoruz
% diğer görünütler için yuzdelik ayarı yapıyoruz
filterIdx = aspectRatio' > 3; 
filterIdx = filterIdx | [mserStats.Eccentricity] > .995 ;
filterIdx = filterIdx | [mserStats.Solidity] < .3;
filterIdx = filterIdx | [mserStats.Extent] < 0.2 | [mserStats.Extent] > 0.9;
filterIdx = filterIdx | [mserStats.EulerNumber] < -4;

% bölgeleri kaldırıyoruz
mserStats(filterIdx) = [];
mserRegions(filterIdx) = [];

%geomerik ozekliklere dayalı metin olmayan bölgeleri kaldırdıktan sonra
figure
imshow(I)
hold on
plot(mserRegions, 'showPixelList', true,'showEllipses',false)
title('After Removing Non-Text Regions Based On Geometric Properties')
hold off
% bir bölgenin ikili goruntusunu alıyoruz ve sınır etkilerinden kaçınmak
% için onu dolduruyoruz
% kontur genişliği hesaplanması sırasında.
regionImage = mserStats(6).Image;
regionImage = padarray(regionImage, [1 1]);

% kontur genişliği görüntüsünü hesaplayın.
distanceImage = bwdist(~regionImage); 
skeletonImage = bwmorph(regionImage, 'thin', inf);

strokeWidthImage = distanceImage;
strokeWidthImage(~skeletonImage) = 0;

% Bölge resmini kontur genişliği resminin yanında göster. 
figure
subplot(1,2,1)
imagesc(regionImage)
title('Region Image')

subplot(1,2,2)
imagesc(strokeWidthImage)
title('Stroke Width Image')
% kontur genişliği değişimini hesaplıyoruz 
strokeWidthValues = distanceImage(skeletonImage);   
strokeWidthMetric = std(strokeWidthValues)/mean(strokeWidthValues);
% eşik belirliyoruz kontur genişliği için
strokeWidthThreshold = 0.4;
strokeWidthFilterIdx = strokeWidthMetric > strokeWidthThreshold;
% kalan bölgeleri işliyoruz
for j = 1:numel(mserStats)
    
    regionImage = mserStats(j).Image;
    regionImage = padarray(regionImage, [1 1], 0);
    
    distanceImage = bwdist(~regionImage);
    skeletonImage = bwmorph(regionImage, 'thin', inf);
    
    strokeWidthValues = distanceImage(skeletonImage);
    
    strokeWidthMetric = std(strokeWidthValues)/mean(strokeWidthValues);
    
    strokeWidthFilterIdx(j) = strokeWidthMetric > strokeWidthThreshold;
    
end

% kontur genişliğine göre bölgelri kaldırıyoruz
mserRegions(strokeWidthFilterIdx) = [];
mserStats(strokeWidthFilterIdx) = [];

% kalan bölgeleri gösteriyoruz
figure
imshow(I)
hold on
plot(mserRegions, 'showPixelList', true,'showEllipses',false)
title('After Removing Non-Text Regions Based On Stroke Width Variation')
hold off
% tüm bölgelr için sınırlayıcı boxlar alıyoruz
bboxes = vertcat(mserStats.BoundingBox);

% [x y width height] e [xmin ymin xmax ymax]  dönüşüm.
xmin = bboxes(:,1);
ymin = bboxes(:,2);
xmax = xmin + bboxes(:,3) - 1;
ymax = ymin + bboxes(:,4) - 1;

% sınırlayıcı kutuları az iktarda genişletiyoruz.
expansionAmount = 0.02;
xmin = (1-expansionAmount) * xmin;
ymin = (1-expansionAmount) * ymin;
xmax = (1+expansionAmount) * xmax;
ymax = (1+expansionAmount) * ymax;

% sınırlayıcı kutuları kırpıyoruz
xmin = max(xmin, 1);
ymin = max(ymin, 1);
xmax = min(xmax, size(I,2));
ymax = min(ymax, size(I,1));

% genişletilmiş sınırlayıcı kutularını gosteriyoruz
expandedBBoxes = [xmin ymin xmax-xmin+1 ymax-ymin+1];
IExpandedBBoxes = insertShape(colorImage,'rectangle',expandedBBoxes,'LineWidth',3);

figure
imshow(IExpandedBBoxes)
title('Expanded Bounding Boxes Text')
% örtüşme oranını hesaplıyoruz
overlapRatio = bboxOverlapRatio(expandedBBoxes, expandedBBoxes);

%Grafik gösterimini basitleştirmek için sınırlayıcı kutu 
% ile kendisi arasındaki örtüşme oranını sıfıra ayarlayın.
n = size(overlapRatio,1); 
overlapRatio(1:n+1:n^2) = 0;

% grafik oluşturuyoruz
g = graph(overlapRatio);

% grafik içindeki bağlantılı metini bul
componentIndices = conncomp(g);
% Kutuları minimum ve maksimimum boyutlara göre birleştiriyoruz.
xmin = accumarray(componentIndices', xmin, [], @min);
ymin = accumarray(componentIndices', ymin, [], @min);
xmax = accumarray(componentIndices', xmax, [], @max);
ymax = accumarray(componentIndices', ymax, [], @max);

%  [x y width height] biçimini kullanarak birleştirilmiş sınırlayıcı
%  kutuları oluşturuyoruz
textBBoxes = [xmin ymin xmax-xmin+1 ymax-ymin+1];
% yalnızca bir metin bölgesi içeren sınıralyıcı kutuları kaldırıyoruz
numRegionsInGroup = histcounts(componentIndices);
textBBoxes(numRegionsInGroup == 1, :) = [];

% show.
ITextRegion = insertShape(colorImage, 'rectangle', textBBoxes,'LineWidth',3);

figure
imshow(ITextRegion)
title('Detected Text')
ocrtxt = ocr(I, textBBoxes);
[ocrtxt.Text]