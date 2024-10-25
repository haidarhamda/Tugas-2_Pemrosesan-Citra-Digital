function bright=freqFiltering(filename, factor)
    img = imread(filename);

    img_fft = fftn(img, size(img));
    img_fftshift = fftshift(img_fft);
    
    [x, y, z] = size(img_fftshift);
    for i=1:z
        img_fftshift(x/2 + 1, y/2 + 1, i) = img_fftshift(x/2 + 1, y/2 + 1, i) * factor;
    end
    
    filtered_fft = ifftshift(img_fftshift);
    bright = real(ifftn(filtered_fft, size(img))) / 255;
end