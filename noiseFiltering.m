function [noise, result]=noiseFiltering(img, noise_density, mean, std_dev, noiseMethod, method, Q, alpha)
    [x, y, z] = size(img);

    switch noiseMethod
        case 'saltAndPepper'
            pixels = x * y;
            salt = round(noise_density * pixels / 2);
            pepper = round(noise_density * pixels / 2);
        
            salt_coords = randperm(pixels, salt);
            pepper_coords = randperm(pixels, pepper);
        
            noise = img;
            for i = 1:z
                noise(salt_coords + (i - 1)*pixels) = 1;
                noise(pepper_coords + (i - 1)*pixels) = 0;
            end
        case 'gauss'
            noise = img + (std_dev * randn(x, y, z) + mean);
    end

    result = real(filter.noise_filter(noise, method, Q, alpha));
end