function writeTiffStack(img,fname)
%img is a matrix of size xyczt with image data
%the function should write this to a multipage tiff that is readable by
%Fiji/ImageJ

%account for either 16-bit or 8-bit images
if isa(img,'uint16')
    nbits = 16;
elseif isa(img,'uint8')
    nbits = 8;
else
    error('only uint16 and uint8 images currently handled');
end
maxstr = ['max=',num2str(intmax(class(img))),'.0'];


fiji_descr = ['ImageJ=' newline ...%do we need to specify ImageJ version?
            'images=' num2str(size(img,3)*size(img,4)*size(img,5)) newline... 
            'channels=' num2str(size(img,3)) newline...
            'slices=' num2str(size(img,4)) newline...
            'frames=' num2str(size(img,5)) newline... 
            'hyperstack=true' newline...
            'mode=grayscale' newline...  
            'loop=false' newline...  
            'min=0.0' newline...      
            maxstr];  % max value
        
% fiji_descr = ["ImageJ=" newline...
%     "hyperstack=true" newline...
%     "images=" num2str(channels*z*t) newline...
%     "channels=" num2str(channels) newline...
%     "slices=" num2str(z) newline...
%     "frames=" num2str(t)];
            
t = Tiff(fname,'w');
tagstruct.ImageLength = size(img,1);
tagstruct.ImageWidth = size(img,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = nbits; %depends on image type
tagstruct.SamplesPerPixel = 1;
tagstruct.Compression = Tiff.Compression.None; %would lossless compression be better?
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
tagstruct.ImageDescription = fiji_descr;

for frame = 1:size(img,5)
    for slice = 1:size(img,4)
        for channel = 1:size(img,3)
            t.setTag(tagstruct)
            t.write(img(:,:,channel,slice,frame));
            t.writeDirectory(); % saves a new page in the tiff file
        end
    end
end
t.close()
