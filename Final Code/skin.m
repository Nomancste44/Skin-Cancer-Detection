
%% Image Processing For Skin Cancer Features Extraction.
function varargout = skin(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @skin_OpeningFcn, ...
                   'gui_OutputFcn',  @skin_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before skin is made visible.
function skin_OpeningFcn(hObject, eventdata, handles, varargin)
%
handles.output = hObject;
clc

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes skin wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = skin_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%Read the Image
% --- Executes on button press in Loadfile.
function Loadfile_Callback(hObject, eventdata, handles)

[FileName,PathName] = uigetfile({'*.*'},'Load Image File');

if (FileName==0) % cancel pressed
    return;
end


handles.fullPath = [PathName FileName];
[a, b, Ext] = fileparts(FileName);
availableExt = {'.bmp','.jpg','.jpeg','.tiff','.png','.gif','.pmp'};
for (i=1:length(availableExt))
    if (strcmpi(Ext, availableExt{i}))
        break;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Implementation of Otsu?s method

%Read the Image
RGB = imread(handles.fullPath);
axes(handles.axes1); cla; imshow(RGB),title('Orginal Image');
% Convert the RGB image into graytlevel intensity image
Igray = rgb2gray(RGB);
axes(handles.axes2); cla; imshow(Igray),title('Gray level Image');



Ibw = im2bw(RGB,graythresh(RGB));
axes(handles.axes3); cla; imshow(Ibw),title('Otsu Segmentation method');

%dilation and erosion morphological operations
se=strel('square',5);
f0=imopen(Ibw,se);
foc=imclose(f0,se);
%Fractal Dimension calculation of the segmented images.
[ D ] = hausDim( foc );
set(handles.text16,'string',[ D ]);

figure
subplot(2,2,[1 3]),imshow(Igray),xlabel('(a)');
%subplot(2,2,2),imshow(Igray),xlabel('(b)');
%subplot(2,2,3),imshow(Ibw),xlabel('(c)');
subplot(2,2,[2 4]),imshow(~foc),xlabel('b)');
set(handles.text21,'string','Skin Lesion(Melanoma) Asymmetry Index Range is From 10% to 25%');
set(handles.text22,'string','Skin Lesion(Melanoma) Compactness Index Range is From 1.50 to 2.00');
set(handles.text23,'string','Skin Lesion(Melanoma) Fractal Dimension Range is From 0.80 to 0.91');
set(handles.text24,'string','Skin Lesion(Melanoma) Edge Abruptness Range is From 0.020 to 0.08');
set(handles.text25,'string','Skin Lesion(Melanoma) Color Variance Range is From 0.100 to 0.155');
set(handles.text28,'string','Skin Lesion(Melanoma) Diameter Range is From 90 to 130 pixels');



% Boundary tracing operation for center finding
dim = size(Ibw);
IN=ones(dim(1),dim(2));
BW=xor(Ibw,IN);  %inverting
%Finding of initial point
row = round(dim(1)/2);
col = min(find(BW(row,:)));
%Tracing
boundary = bwtraceboundary(BW,[row, col],'w');
%subplot(2,2,4),imshow(I), title(?Traced?);
%hold on;
%Display traced boundary
figure
plot(boundary(:,2),boundary(:,1),'g','LineWidth',2);
hold off
% figure
% plot(boundary(:,2),boundary(:,1),?black?,'LineWidth?,2);

nn=size(boundary);
KM=zeros(dim(1),dim(2));
 ii=0;
 %Create new matrix with boundary points. there fore we can get rid off
 %other distortions outside boundaries
 while ii<nn(1)
     ii=ii+1;
    KM(boundary(ii,1),boundary(ii,2))=1;
end
 figure
 subplot(2,2,1),plot(boundary(:,2),boundary(:,1),'black','LineWidth',2);
 subplot(2,2,2),imshow(KM)
%Fill inner boundaries where lesion is located
KM2 = imfill(KM,'holes');
subplot(2,2,3),imshow(KM2)
KM1=xor(KM2,IN);
subplot(2,2,4),imshow(KM1)
%Geometrical center
IVx=[1:dim(2)];
IVy=[1:dim(1)];
IMx=ones(dim(1),1)*IVx;
IMy=ones(dim(2),1)*IVy;
IMy = imrotate(IMy,-90);
Koordx=IMx.*KM2;
Koordy=IMy.*KM2;
xmean=mean(Koordx,2);
yc=round(sum(xmean.*IMy(:,1))/sum(xmean));
ymean=mean(Koordy);
xc=round(sum(ymean.*IVx)/sum(ymean));
figure
plot(boundary(:,2),boundary(:,1),'green','LineWidth',2);
hold on
plot(xc,1:dim(1),'red','LineWidth',2);
plot(1:dim(2),yc,'red','LineWidth',2);
hold off


ll=length(boundary);
d=0;
for ii=1:ll
    d =d+sqrt((xc-boundary(ii,1))^2+(yc-boundary(ii,2))^2);
end
 
d_mean= d/ll;
d_sq=0;
for ii=1:ll
    d_sq=d_sq+((sqrt((xc-boundary(ii,1))^2+(yc-boundary(ii,2))^2))-d_mean)^2;
end

axes(handles.axes4); cla; imshow(~foc);

hold on

BW=bwlabel(foc);
% Area, Diamter and Perimeter Calculation
s  = regionprops(BW, 'Area','EquivDiameter', 'Perimeter');
BW1 = bwlabel(Igray);
s2  = regionprops(BW1, 'Area');
dia = s.EquivDiameter;
set(handles.text27,'string',dia);

A_tot = s2.Area;
A_ob = s.Area;
PL=s.Perimeter;

%Asymmetryy Index Calculation
Asymetry_index=(0.5*(A_tot-A_ob)/A_tot)*100;
set(handles.text1,'string',Asymetry_index);
%Compactness Index Calculation

Compact_index=(PL)^2/(4*pi*A_ob);
set(handles.text4,'string',Compact_index);


% Barrier Irregularity Calculation
B_irr= d_sq/(PL*(d_mean^2));
set(handles.text19,'string',B_irr);


%Color Varience Calculation
im = rgb2hsv(RGB);
V=var(im(:));
set(handles.text13,'string',V);

%Diameter Calculation
dia = s.EquivDiameter;
set(handles.text27,'string',dia);
%TDS
if (Asymetry_index>25)
    Asymetry_index=25;
end
    A1 = (3.8/25)*(Asymetry_index)*(1.3);

if (B_irr>0.08)
    B_irr=0.08;
end
if (Compact_index>2.00)
    Compact_index=2.00;
end
if (D>0.91)
    D=0.91;
end
B12=(B_irr+Compact_index+D)/3;
B1=(6.2/0.97266)*(B12*0.1);


if (V>0.155)
    V=0.155;
end
C1=((3/0.155)*(V)*(0.5));
if (dia>130)
    dia=130;
end
D1=((2.5/130)*(dia)*(0.5));
TDS=A1+B1+C1+D1;

if (TDS>=5.45)
    fprintf('Asymmetry\tBoarder_Irregularities\tColor_Variation\t\tDiameter\t\tTDS\t\tResult');
    fprintf('\n%.2f\t\t\t%.2f\t\t\t\t\t%.2f\t\t\t%.2f\t\t\t%.2f',A1,B1,C1,D1,TDS);
    fprintf('\tMELANOMA\n\n\n\n\n\n');
  figure,imshow(RGB),text(2,5,'This is Melanoma','Color',[1 0 0],'FontName','TimesNewRoman','FontSize',15);
else
   fprintf('Asymmetry\tBoarder_Irregularities\tColor_Variation\t\tDiameter\t\tTDS\t\tResult');
    fprintf('\n%.2f\t\t\t%.2f\t\t\t\t\t%.2f\t\t\t%.2f\t\t\t%.2f',A1,B1,C1,D1,TDS);
    fprintf('\tBENIGN\n\n\n\n\n\n');
    figure,imshow(RGB),text(2,5,'This is Benign','Color',[0 1 0],'FontName','TimesNewRoman','FontSize',15);
end  
%% Implementation of Colored-based Segmentation using k-means Clustering.
% RGB color space to L*a*b* color space
C = makecform('srgb2lab');
lab = applycform(RGB,C);
ab  =double(lab(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab  = reshape(ab,nrows*ncols,2);
ncolors =3;

%Classify the colors in ?a*b*? space using k-mean clustering condition
[cluster_idx cluster_center]= kmeans(ab,ncolors,'distance','sqEuclidean','Replicates',3);

%Label every pixel in the image using the results from k-means

pixel_labels = reshape(cluster_idx,nrows,ncols);
axes(handles.axes7); cla; imshow(pixel_labels,[]), title('image labeled by cluster index');

%Create images that segments the image by color
segmented_images = cell(1,3);
rgb_label = repmat(pixel_labels,[1 1 3]);

for k = 1 : ncolors
    color = RGB;
    color(rgb_label ~= k) = 0;
    segmented_images{k} = color;
end;
axes(handles.axes8); cla; imshow(segmented_images{1}), title('objects in cluster 1'); 
axes(handles.axes9); cla; imshow(segmented_images{2}), title('objects in cluster 2'); 
axes(handles.axes10); cla; imshow(segmented_images{3}), title('objects in cluster 3'); 

%Segment the nuclei into separate image.
mean_cluster_value = mean(cluster_center,2);
[tmp, idx] = sort(mean_cluster_value);
blue_cluster_num = idx(1);

L = lab(:,:,1);
blue_idx = find(pixel_labels == blue_cluster_num);
L_blue = L(blue_idx);
is_light_blue = im2bw(L_blue,graythresh(L_blue));
nuclei_labels = repmat(uint8(0),[nrows ncols]);
nuclei_labels(blue_idx(is_light_blue==false)) = 1;
nuclei_labels = repmat(nuclei_labels,[1 1 3]);
blue_nuclei = RGB;
blue_nuclei(nuclei_labels ~= 1) = 0;
axes(handles.axes11); cla; imshow(blue_nuclei),title('blue nuclei'); 


se=strel('square',10);
f03=imdilate(blue_nuclei,se);
foc3=imerode(f03,se);
axes(handles.axes12); cla; imshow(foc3),title('blue nuclei after Morphological operation');
figure
subplot(2,2,1),imshow(RGB),xlabel('(a)');
subplot(2,2,2),imshow(segmented_images{1}),xlabel('(b)');
subplot(2,2,3),imshow(segmented_images{2}),xlabel('(c)');
subplot(2,2,4),imshow(segmented_images{3}),xlabel('(d)');




%% Implementation of GVF Method
[u,v] = GVF(Igray, 0.25, 100);
axes(handles.axes5); cla; imshow(u), title('Gradient Vector Flow(GVF)Segmentation method');


se=strel('square', 7);
f01=imdilate(u,se);
foc1=imerode(f01,se);
axes(handles.axes6); cla; imshow(foc1);
figure
subplot(2,2,1),imshow(RGB),xlabel('(a)');
subplot(2,2,2),imshow(Igray),xlabel('(b)');
subplot(2,2,[3 4]),imshow(u),xlabel('(c)');
%subplot(2,2,4),imshow(foc1),xlabel('(d)');




function B = BoundMirrorExpand(A)

[m,n] = size(A);
yi = 2:m+1;
xi = 2:n+1;
B = zeros(m+2,n+2);
B(yi,xi) = A;
B([1 m+2],[1 n+2]) = B([3 m],[3 n]);  % mirror corners
B([1 m+2],xi) = B([3 m],xi);          % mirror left and right boundary
B(yi,[1 n+2]) = B(yi,[3 n]);          % mirror top and bottom boundary
function B = BoundMirrorShrink(A)

[m,n] = size(A);
yi = 2:m-1;
xi = 2:n-1;
B = A(yi,xi);
function B = BoundMirrorEnsure(A)
[m,n] = size(A);

if (m<10 || n<10) 
    error('either the number of rows or columns is smaller than 2');
end
yi = 2:m-1;
xi = 2:n-1;
B = A;
B([1 m],[1 n]) = B([3 m-2],[3 n-2]);  % mirror corners
B([1 m],xi) = B([3 m-2],xi);          % mirror left and right boundary
B(yi,[1 n]) = B(yi,[3 n-2]);          % mirror top and bottom boundary

function [u,v] = GVF(Igray, mu, ITER)


[m,n] = size(Igray);
fmin  = min(Igray(:));
fmax  = max(Igray(:));
f = (Igray-fmin)/(fmax-fmin);  % Normalize f to the range [0,1]

f = BoundMirrorExpand(f);  % Take care of boundary condition
[fx,fy] = gradient(f);     % Calculate the gradient of the edge map
u = fx; v = fy;            % Initialize GVF to the gradient
SqrMagf = fx.*fx + fy.*fy; % Squared magnitude of the gradient field

% Iteratively solve for the GVF u,v
for i=1:ITER,
  u = BoundMirrorEnsure(u);
  v = BoundMirrorEnsure(v);
  u = u + mu*4*del2(u) - SqrMagf.*(u-fx);
  v = v + mu*4*del2(v) - SqrMagf.*(v-fy);
  fprintf(1, '%3d', i);
  if (rem(i,20) == 0)
     fprintf(1, '\n');
  end 
end
fprintf(1, '\n');

u = BoundMirrorShrink(u);
v = BoundMirrorShrink(v);


%% Hausdorff Fractal Dimension Calculation Funcation
function [ D ] = hausDim(I )
    % Pad the image with background pixels so that its dimensions are a power of 2.
    maxDim = max(size(I));
    newDimSize = 2^ceil(log2(maxDim));
    rowPad = newDimSize - size(I, 1);
    colPad = newDimSize - size(I, 2);
    I = padarray(I, [rowPad, colPad], 'post');

    boxCounts = zeros(1, ceil(log2(maxDim)));
    resolutions = zeros(1, ceil(log2(maxDim)));
    %Set the box size ?e? to the size of the image
    boxSize = size(I, 1);
    boxesPerDim = 1;
    idx = 0;
    while boxSize >= 1
        boxCount = 0;
        
        for boxRow = 1:boxesPerDim
            for boxCol = 1:boxesPerDim
                minRow = (boxRow - 1) * boxSize + 1;
                maxRow = boxRow * boxSize;
                minCol = (boxCol - 1) * boxSize + 1;
                maxCol = boxCol * boxSize;
                
                objFound = false;
                for row = minRow:maxRow
                    for col = minCol:maxCol
                        if I(row, col)
                            boxCount = boxCount + 1;
                            objFound = true; % Break from nested loop.
                        end;
                        
                        if objFound
                            break; % Break from nested loop.
                        end;
                    end;
                    
                    if objFound
                        break; % Break from nested loop.
                    end;
                end;
            end;
    %Compute N (e), which corresponds to the number of boxes of size 'e'  
        
        idx = idx + 1;
        boxCounts(idx) = boxCount;
        resolutions(idx) = 1 / boxSize;
        
        boxesPerDim = boxesPerDim * 2;
        boxSize = boxSize / 2;
    end;
    %Compute the points log (N (e)) x log (1/e) and use the least squares method to fit a line to the points.
    
    D = polyfit(log(resolutions), log(boxCounts), 1);
    %The returned Haussdorf fractal dimension D is the slope of the line
    D = D(1);
    
    end
    
function pushbutton2_Callback(hObject, eventdata, handles)

close gcf;



