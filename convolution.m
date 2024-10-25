function [h]=convolution(f,g)
    [n, m, c] = size(f);
    [x,y]=size(g);
    h=double(zeros(n,m,c));
    padx=floor(x/2);
    pady=floor(y/2);
    paddedf=padarray(f,[padx,pady]);

    for channel=1:c
        imgChannel = paddedf(:,:,channel);
        tmpChannel=double(zeros(n,m));
        for i=2:n
            for j=2:m
                tmp=double(imgChannel(i:i+x-1,j:j+y-1));
                tmpChannel(i,j)=sum(sum(tmp.*double(g)));
            end
        end
        h(:,:,channel)=tmpChannel;
    end
    
end

% img = imread("img\shore.jpg");
% g=[-1 -1 -1; -1 8 -1; -1 -1 -1];
% h=(convolution(im2double(img),g));
% imshow(h)
% hm=uint8(conv2(im2double(img), double(g)));
% figure,imshow(hm)