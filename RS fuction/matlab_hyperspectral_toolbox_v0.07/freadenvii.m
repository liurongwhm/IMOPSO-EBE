function [image_2d,image_3d,p,t] = freadenvii(fname);
% freadenvit - read ENVI image (V. Guissard, Apr 29, 2004)
%
%     Reads an image of ENVI standard type 
%    to a [(lines x samples) x bands] MATLAB array
% and transposes all bands (Allan A. Nielsen)
%
% Syntax
%
% image = freadenvit(fname);
% [image,p] = freadenvit(fname);
% [image,p,t] = freadenvit(fname);
%
% INPUT :
%
% fname string giving the full pathname of the ENVI image to read
%
% OUTPUT :
%
% image lines by samples by bands array containing the ENVI image values organised in
%    (lines x samples) x bands.
% p  1 by 3 vector that contains (1) the number of samples, (2) the
%       number of liness and (3) the number of bands of the opened image.
%
% t  string describing the image data type in MATLAB convention.
%
% NOTE : freadenvi needs the corresponding image header file generated
%   automatically by ENVI. The ENVI header file must have the same name
%   as the ENVI image file + the '.hdr' extension.
% Modified (very moderately) including tranposing all frames by
% Allan Aasbjerg Nielsen
% aa@imm.dtu.dk
% Dec 2008
% Parameters initialization
elements = {'samples ' 'lines   ' 'bands   ' 'data type ' 'interleave '};
d = {'uint8' 'int16' 'int32' 'float32' 'float64' 'uint16' 'uint32' 'int64' 'uint64'};
%d = {'*uint8' 'int16' 'int32' 'float32' 'float64' 'uint16' 'uint32' 'int64' 'uint64'};
% Check user input
if ~ischar(fname)
    error('fname should be a char string');
end
% Open ENVI header file to retreive s, l, b & d variables
rfid = fopen(strcat(fname,'.hdr'),'r');

% Check if the header file is correctely opened
if rfid == -1
    error('Input header file does not exist');
end;
% Read ENVI image header file and get p(1) : nb samples,
% p(2) : nb lines, p(3) : nb bands and t : data type
while 1
    tline = fgetl(rfid);
    if ~ischar(tline), break, end
    [first,second] = strtok(tline,'=');
    
    switch first
        case elements(1)
            [f,s] = strtok(second);
            p(1) = str2num(s);
        case elements(2)
            [f,s] = strtok(second);
            p(2) = str2num(s);
        case elements(3)
            [f,s] = strtok(second);
            p(3) = str2num(s);
        case elements(4)
            [f,s] = strtok(second);
            t = str2num(s);
            switch t
                case 1
                    t = d(1);
                case 2
                    t = d(2);
                case 3
                    t = d(3);
                case 4
                    t = d(4);
                case 5
                    t = d(5);
                case 12
                    t = d(6);
                case 13
                    t = d(7);
                case 14
                    t = d(8);
                case 15
                    t = d(9);
                otherwise
                    error('Unknown image data type');
            end
        case elements(5) % input must be BSQ here
            [f,s] = strtok(second);
            if s(1,1)~=' bsq', error('input must be bsq'); end
    end
end
fclose(rfid); 
t = t{1,1};
% Open the ENVI image and store it in the 'image' MATLAB array
%disp(['Input is (',(num2str(p(2))),' lines) x (',(num2str(p(1))),' samples) x (',(num2str(p(3))),' bands)',' type ', (t), (' image ...')]);
disp([('Opening '),(num2str(p(1))),('cols x '),(num2str(p(2))),('lines x '),(num2str(p(3))),('bands')]);
disp([('of type '), (t), (' image...')]);
fid = fopen(fname,'r');
if fid==-1, error(strcat(fname,' not found')); end

image_2d=fread(fid,[p(1)*p(2),p(3)],t);% original freadenvi reshape
fclose(fid);

image_3d = zeros(p(2),p(1),p(3));
for ii=1:p(3)
    image_3d(:,:,ii) = (reshape(image_2d(:,ii),p(1),p(2)))';
end
% clear image1;
