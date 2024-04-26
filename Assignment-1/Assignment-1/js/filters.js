"use strict";

const Filters = {};

////////////////////////////////////////////////////////////////////////////////
// General utility functions
////////////////////////////////////////////////////////////////////////////////

// Hardcoded Pi value
// const pi = 3.14159265359;
const pi = Math.PI;

// Constrain val to the range [min, max]
function clamp(val, min, max) {
    /* Shorthand for:
    * if (val < min) {
    *   return min;
    * } else if (val > max) {
    *   return max;
    * } else {
    *   return val;
    * }
    */
    return val < min ? min : val > max ? max : val;
}

// Extract vertex coordinates from a URL string
function stringToCoords(vertsString) {
    const centers = [];
    const coordStrings = vertsString.split("x");
    for (let i = 0; i < coordStrings.length; i++) {
        const coords = coordStrings[i].split("y");
        const x = parseInt(coords[0]);
        const y = parseInt(coords[1]);
        if (!isNaN(x) && !isNaN(y)) {
            centers.push({ x: x, y: y });
        }
    }

    return centers;
}

// Blend scalar start with scalar end. Note that for image blending,
// end would be the upper layer, and start would be the background
function blend(start, end, alpha) {
    return start * (1 - alpha) + end * alpha;
}

// ----------- STUDENT CODE BEGIN ------------
// ----------- Our reference solution uses 72 lines of code.
// ----------- STUDENT CODE END ------------

////////////////////////////////////////////////////////////////////////////////
// Filters
////////////////////////////////////////////////////////////////////////////////

// You've already implemented this in A0! Feel free to copy your code into here
Filters.fillFilter = function(image, color) {
    for (var x = 0; x < image.width; x++) {
        for (var y = 0; y < image.height; y++) {
          // uncomment this line to enable this function
          image.setPixel(x, y, color);
        }
      }
      return image;
};

// You've already implemented this in A0! Feel free to copy your code into here
Filters.brushFilter = function(image, radius, color, vertsString) {
    // centers is an array of (x, y) coordinates that each defines a circle center
    const centers = stringToCoords(vertsString);

    // draw a filled circle centered at every location in centers[].
    // radius and color are specified in function arguments.
    // ----------- STUDENT CODE BEGIN ------------
     // centers is an array of (x, y) coordinates that each defines a circle center


    // draw a filled circle centered at every location in centers[].
    // radius and color are specified in function arguments.
    // ----------- STUDENT CODE BEGIN ------------
    for (var i = 0; i < centers.length; i++) {
        var center = centers[i]
        
        for (var x = center.x - radius; x <= center.x + radius; x++) {
        for (var y = center.y - radius; y <= center.y + radius; y++) {
            var dx = x - center.x
            var dy = y - center.y 

            var d = Math.sqrt((dx * dx) + (dy * dy))

            if (d < radius) {
            image.setPixel(x, y, color)
            }
        }
        }
    }
    // ----------- STUDENT CODE END ------------
    //Gui.alertOnce ('brushFilter is not implemented yet');

    return image;
    // ----------- STUDENT CODE END ------------
};

// You've already implemented this in A0! Feel free to copy your code into here
Filters.softBrushFilter = function(image, radius, color, alpha_at_center, vertsString) {
    // centers is an array of (x, y) coordinates that each defines a circle center
    const centers = stringToCoords(vertsString);

    // draw a filled circle with opacity equals to alpha_at_center at the center of each circle
    // the opacity decreases linearly along the radius and becomes zero at the edge of the circle
    // radius and color are specified in function arguments.
    // ----------- STUDENT CODE BEGIN ------------
    for (var i = 0; i < centers.length; i++) {
        var center = centers[i]
          
        for (var x = center.x - radius; x <= center.x + radius; x++) {
          for (var y = center.y - radius; y <= center.y + radius; y++) {
            var dx = x - center.x
            var dy = y - center.y 
    
            var distance = Math.sqrt((dx * dx) + (dy * dy))
    
            if (distance < radius) {
              var Height = radius * 2
              var Width = radius * 2
              color.a = alpha_at_center - distance / (Math.sqrt((Width / 2) * (Width / 2) + (Height / 2) * (Height / 2)))
              image.setPixel(x, y, color)
            }
          }
        }
      }
    // ----------- STUDENT CODE END ------------

    return image;
};

// Ratio is a value in the domain [-1, 1]. When ratio is < 0, linearly blend the image
// with black. When ratio is > 0, linearly blend the image with white. At the extremes
// of -1 and 1, the image should be completely black and completely white, respectively.
Filters.brightnessFilter = function(image, ratio) {
    let alpha, dirLuminance;
    if (ratio < 0.0) {
        alpha = 1 + ratio;
        dirLuminance = 0; // blend with black
    } else {
        alpha = 1 - ratio;
        dirLuminance = 1; // blend with white
    }

    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            pixel.data[0] = alpha * pixel.data[0] + (1 - alpha) * dirLuminance;
            pixel.data[1] = alpha * pixel.data[1] + (1 - alpha) * dirLuminance;
            pixel.data[2] = alpha * pixel.data[2] + (1 - alpha) * dirLuminance;

            image.setPixel(x, y, pixel);
        }
    }

    return image;
};

// Reference at this:
//      https://en.wikipedia.org/wiki/Image_editing#Contrast_change_and_brightening
// value = (value - 0.5) * (tan ((contrast + 1) * PI/4) ) + 0.5;
// Note that ratio is in the domain [-1, 1]
Filters.contrastFilter = function(image, ratio) {
    // ----------- STUDENT CODE BEGIN ------------
    function ajustContrast(value) {
        let ajustedValue = value
        if (ratio < 0.0)  {
            ajustedValue = value * ( 1.0 + ratio)
        } else {
            ajustedValue = value + ((1 - value) * ratio);
        }
        ajustedValue = (ajustedValue - 0.5) * Math.tan((ratio + 1) * Math.PI/4) + 0.5

        return ajustedValue
    }

    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y)

            pixel.data[0] = ajustContrast(pixel.data[0])
            pixel.data[1] = ajustContrast(pixel.data[1])
            pixel.data[2] = ajustContrast(pixel.data[2])

            image.setPixel(x, y, pixel)
        }
    }
    
    // ----------- STUDENT CODE END ------------
    //Gui.alertOnce ('contrastFilter is not implemented yet');
    return image;
};

// Note that the argument here is log(gamma)
Filters.gammaFilter = function(image, logOfGamma) {
    const gamma = Math.exp(logOfGamma);
    // ----------- STUDENT CODE BEGIN ------------
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y)

            pixel.data[0] = Math.pow(pixel.data[0], gamma)
            pixel.data[1] = Math.pow(pixel.data[1], gamma)
            pixel.data[2] = Math.pow(pixel.data[2], gamma)
            
            image.setPixel(x, y, pixel)

        }
    }
    // ----------- STUDENT CODE END ------------
    //Gui.alertOnce ('gammaFilter is not implemented yet');
    return image;
};

/*
* The image should be perfectly clear up to innerRadius, perfectly dark
* (black) at outerRadius and beyond, and smoothly increase darkness in the
* circular ring in between. Both are specified as multiples of half the length
* of the image diagonal (so 1.0 is the distance from the image center to the
* corner).
*
* Note that the vignette should still form a perfect circle!
*/
Filters.vignetteFilter = function(image, innerR, outerR) {
    // Let's ensure that innerR is at least 0.1 smaller than outerR
    innerR = clamp(innerR, 0, outerR - 0.1);
    // ----------- STUDENT CODE BEGIN ------------
        
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('vignetteFilter is not implemented yet');
    return image;
};

/*
* You will want to build a normalized CDF of the L channel in the image.
*/
Filters.histogramEqualizationFilter = function(image) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 33 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('histogramEqualizationFilter is not implemented yet');
    return image;
};

// Set each pixel in the image to its luminance
Filters.grayscaleFilter = function(image) {
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);
            const luminance = 0.2126 * pixel.data[0] + 0.7152 * pixel.data[1] + 0.0722 * pixel.data[2];
            pixel.data[0] = luminance;
            pixel.data[1] = luminance;
            pixel.data[2] = luminance;

            image.setPixel(x, y, pixel);
        }
    }

    return image;
};

// Adjust each channel in each pixel by a fraction of its distance from the average
// value of the pixel (luminance).
// See: http://www.graficaobscura.com/interp/index.html
Filters.saturationFilter = function(image, ratio) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 13 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('saturationFilter is not implemented yet');
    return image;
};

// Apply the Von Kries method: convert the image from RGB to LMS, divide by
// the LMS coordinates of the white point color, and convert back to RGB.
Filters.whiteBalanceFilter = function(image, white) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 23 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('whiteBalanceFilter is not implemented yet');
    return image;
};

// This is similar to the histogram filter, except here you should take the
// the CDF of the L channel in one image and
// map it to another
//
Filters.histogramMatchFilter = function(image, refImg) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 58 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('histogramMatchFilter is not implemented yet');
    return image;
};

// Convolve the image with a gaussian filter.
// NB: Implement this as a seperable gaussian filter
Filters.gaussianFilter = function(image, sigma) {
    // note: this function needs to work in a new copy of the image
    //       to avoid overwriting original pixels values needed later
    // create a new image with the same size as the input image
    let newImg = image.createImg(image.width, image.height);
    // the filter window will be [-winR, winR] for a total diameter of roughly Math.round(3*sigma)*2+1;
    const winR = Math.round(sigma * 3);
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 58 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('gaussianFilter is not implemented yet');
    return newImg;
};

/*
* First the image with the edge kernel and then add the result back onto the
* original image.
*/
Filters.sharpenFilter = function(image) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 33 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('sharpenFilter is not implemented yet');
    return image;
};

/*
* Convolve the image with the edge kernel from class. You might want to define
* a convolution utility that convolves an image with some arbitrary input kernel
*
* For this filter, we recommend inverting pixel values to enhance edge visualization
*/
Filters.edgeFilter = function(image) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 57 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('edgeFilter is not implemented yet');
    return image;
};

// Set a pixel to the median value in its local neighbor hood. You might want to
// apply this seperately to each channel.
Filters.medianFilter = function(image, winR) {
    // winR: the window will be  [-winR, winR];
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 36 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('medianFilter is not implemented yet');
    return image;
};

// Apply a bilateral filter to the image. You will likely want to reference
// precept slides, lecture slides, and the assignments/examples page for help.
Filters.bilateralFilter = function(image, sigmaR, sigmaS) {
    // reference: https://en.wikipedia.org/wiki/Bilateral_filter
    // we first compute window size and preprocess sigmaR
    const winR = Math.round((sigmaR + sigmaS) * 1.5);
    sigmaR = sigmaR * (Math.sqrt(2) * winR);

    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 53 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('bilateralFilter is not implemented yet');
    return image;
};

// Conver the image to binary
Filters.quantizeFilter = function(image) {
    // convert to grayscale
    image = Filters.grayscaleFilter(image);

    // use center color
    for (let i = 0; i < image.height; i++) {
        for (let j = 0; j < image.width; j++) {
            const pixel = image.getPixel(j, i);
            for (let c = 0; c < 3; c++) {
                pixel.data[c] = Math.round(pixel.data[c]);
            }
            pixel.clamp();
            image.setPixel(j, i, pixel);
        }
    }
    return image;
};

// To apply random dithering, first convert the image to grayscale, then apply
// random noise, and finally quantize
Filters.randomFilter = function(image) {
    // convert to grayscale
    image = Filters.grayscaleFilter(image);

    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 12 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('randomFilter is not implemented yet');
    return image;
};

// Apply the Floyd-Steinberg dither with error diffusion
Filters.floydFilter = function(image) {
    // convert to grayscale
    image = Filters.grayscaleFilter(image);

    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 27 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('floydFilter is not implemented yet');
    return image;
};

// Apply ordered dithering to the image. We recommend using the pattern from the
// examples page and precept slides.
Filters.orderedFilter = function(image) {
    // convert to gray scale
    image = Filters.grayscaleFilter(image);

    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 31 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('orderedFilter is not implemented yet');
    return image;
};

// Implement bilinear and Gaussian sampling (in addition to the basic point sampling).
// This operation doesn't appear on GUI and should be used as a utility function.
// Call this function from filters that require sampling (e.g. scale, rotate)
Filters.samplePixel = function(image, x, y, mode) {
    if (mode === "bilinear") {
        // ----------- STUDENT CODE BEGIN ------------
        // ----------- Our reference solution uses 21 lines of code.
        // ----------- STUDENT CODE END ------------
        Gui.alertOnce ('bilinear sampling is not implemented yet');
    } else if (mode === "gaussian") {
        // ----------- STUDENT CODE BEGIN ------------
        // ----------- Our reference solution uses 38 lines of code.
        // ----------- STUDENT CODE END ------------
        Gui.alertOnce ('gaussian sampling is not implemented yet');
    } else {
        // point sampling
        y = Math.max(0, Math.min(Math.round(y), image.height - 1));
        x = Math.max(0, Math.min(Math.round(x), image.width - 1));
        return image.getPixel(x, y);
    }
};

// Translate the image by some x, y and using a requested method of sampling/resampling
Filters.translateFilter = function(image, x, y, sampleMode) {
    // Note: set pixels outside the image to RGBA(0,0,0,0)
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 21 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('translateFilter is not implemented yet');
    return image;
};

// Scale the image by some ratio and using a requested method of sampling/resampling
Filters.scaleFilter = function(image, ratio, sampleMode) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 19 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('scaleFilter is not implemented yet');
    return image;
};

// Rotate the image by some angle and using a requested method of sampling/resampling
Filters.rotateFilter = function(image, radians, sampleMode) {
    // Note: set pixels outside the image to RGBA(0,0,0,0)
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 29 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('rotateFilter is not implemented yet');
    return image;
};

// Swirl the filter about its center. The rotation of the swirl should be in linear increase
// along the radial axis up to radians
Filters.swirlFilter = function(image, radians, sampleMode) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 26 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('swirlFilter is not implemented yet');
    return image;
};

// Set alpha from luminance
Filters.getAlphaFilter = function(backgroundImg, foregroundImg) {
    for (let i = 0; i < backgroundImg.height; i++) {
        for (let j = 0; j < backgroundImg.width; j++) {
            const pixelBg = backgroundImg.getPixel(j, i);
            const pixelFg = foregroundImg.getPixel(j, i);
            const luminance =
            0.2126 * pixelFg.data[0] + 0.7152 * pixelFg.data[1] + 0.0722 * pixelFg.data[2];
            pixelBg.a = luminance;
            backgroundImg.setPixel(j, i, pixelBg);
        }
    }

    return backgroundImg;
};

// Composites the foreground image over the background image, using the alpha
// channel of the foreground image to blend two images.
Filters.compositeFilter = function(backgroundImg, foregroundImg) {
    // Assume the input images are of the same sizes.
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 14 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('compositeFilter is not implemented yet');
    return backgroundImg;
};

// Morph two images according to a set of correspondance lines
Filters.morphFilter = function(initialImg, finalImg, alpha, sampleMode, linesFile) {
    const lines = Parser.parseJson("images/" + linesFile);

    // The provided linesFile represents lines in a flipped x, y coordinate system
    //  (i.e. x for vertical direction, y for horizontal direction).
    // Therefore we first fix the flipped x, y coordinates here.
    for (let i = 0; i < lines.initial.length; i++) {
        [lines.initial[i].x0, lines.initial[i].y0] = [lines.initial[i].y0, lines.initial[i].x0];
        [lines.initial[i].x1, lines.initial[i].y1] = [lines.initial[i].y1, lines.initial[i].x1];
        [lines.final[i].x0, lines.final[i].y0] = [lines.final[i].y0, lines.final[i].x0];
        [lines.final[i].x1, lines.final[i].y1] = [lines.final[i].y1, lines.final[i].x1];
    }

    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 114 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('morphFilter is not implemented yet');
    return image;
};

// Use k-means to extract a pallete from an image
Filters.paletteFilter = function(image, colorNum) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 89 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('paletteFilter is not implemented yet');
    return image;
};

// Read the following paper and implement your own "painter":
//      http://mrl.nyu.edu/publications/painterly98/hertzmann-siggraph98.pdf
Filters.paintFilter = function(image, value) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 59 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('paintFilter is not implemented yet');
    return image;
};

/*
* Read this paper for background on eXtended Difference-of-Gaussians:
*      http://www.cs.princeton.edu/courses/archive/spring19/cos426/papers/Winnemoeller12.pdf
* Read this paper for an approach that develops a flow field based on a bilateral filter
*      http://www.cs.princeton.edu/courses/archive/spring19/cos426/papers/Kang09.pdf
*/
Filters.xDoGFilter = function(image, value) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 70 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('xDoGFilter is not implemented yet');
    return image;
};

// You can use this filter to do whatever you want, for example
// trying out some new idea or implementing something for the
// art contest.
// Currently the 'value' argument will be 1 or whatever else you set
// it to in the URL. You could use this value to switch between
// a bunch of different versions of your code if you want to
// code up a bunch of different things for the art contest.
Filters.customFilter = function(image, value) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 0 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('customFilter is not implemented yet');
    return image;
};
