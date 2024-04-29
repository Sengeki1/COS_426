  ### Constrast

To apply contrast filter its really simple. ```contrastFilter(image, ratio)``` changes the contrast of an image by interpolating between a constant gray image ```(ratio = -1)``` with the average luminance and the original image ```(ratio = 0)```. Interpolation reduces contrast, extrapolation boosts contrast, and negative factors generate inverted images.

Noting that ratio is in the domain ```[-1, 1]``` we only have to change the color channels by ajusting the contrast which we use this formula as reference.

```js
  if (ratio < 0.0)  {
      ajustedValue = value * ( 1.0 + ratio)
  } else {
      ajustedValue = value + ((1 - value) * ratio)
  }
  ajustedValue = (ajustedValue - 0.5) * Math.tan((ratio + 1) * Math.PI/4) + 0.5
```

### Vignette

The vignetting effect typically involves smoothly darkening the pixels of an image towards the edges, while maintaining clarity in the center. The inner and outer radii define the boundaries within which this effect occurs.

For each pixel in the image, calculate its distance from the center of the image. This distance can be computed using the Euclidean distance formula.

$$r = \sqrt{(x - x_c)^2 + (y - y_c)^2}$$

Where $(x_c, y_c)$ are the coordinates of the center of the image. 

Using a Quadratic function we can smoothly transitions from 1 (no effect) at the center to 0 (fully darkened) at the outer radius. 

$$f(r) = 1 - (\frac{r - innerR}{outerR - innerR})^2$$

This function smoothly decreases from ```1``` to ```0``` as ```r``` increases from the inner radius to the outer radius.
We then normalize it to ensure that points are applied properly within the circular region

### Histogram Equalization

<img src="https://github.com/Sengeki1/JS_Paint_COS426/assets/106749775/bd0ff0f2-7c8b-4d5c-b476-d5a17dc0d0f1" alt="Logo">

Here, ```δ``` is the Kronecker delta function, which is ```1``` if $I(x, y) = k$ and ```0``` otherwise. This formula computes the histogram H(k) by counting the frequency of each intensity level k in the image.

* $H(k)$ represents the histogram value for intensity level $k$
* $N$ and $M$ are the dimensions of the image.
* $I(x, y)$ represents the intensity value of the image at position $(x, y)$
* $δ(I(x, y) - k)$ is the Kronecker delta function, which equals 1 if $I(x, y)$ = k and 0 otherwise.

First we iterate over the pixels then we retrive the intensity value $I(x, y)$ or Lightness of the Image of pixel at coordinate $(x, y)$
The line ```histogram[lightness]++```; increments the count in the histogram array for the specific lightness value encountered at a pixel. Let me clarify why we're using this line:

In the histogram array, each index represents a particular lightness value, ranging from ```0``` to ```100```. For example, ```histogram[0]``` represents the count of pixels with a lightness value of ```0```, ```histogram[1]``` represents the count of pixels with a lightness value of ```1```, and so on, up to ```histogram[100]``` for a lightness value of ```100```.

When we encounter a pixel with a particular lightness value, we increment the count in the histogram array at the index corresponding to that lightness value. For example, if we encounter a pixel with a ```lightness value of 50```, we increment ```histogram[50]``` by ```1``` to indicate that we've found ```another pixel with a lightness value of 50```.

After computing the histogram, the next step in histogram equalization is to compute the ```Cumulative Distribution Function (CDF)``` based on the histogram. The CDF represents the cumulative probability distribution of the lightness values in the image.

First, we compute the histogram $H(k)$, which counts the frequency of each intensity level $k$ in the image. Diferent from grayscale which its levels range from ```[0, 255]```, hsl range from ```0%``` to ```100%```.For example, let's say our histogram is:
* $H (0) = 10$
* $H (1) = 20$
* $H (2) = 15$
* ...
* $H (100) = 5$
  
Next, we normalize the histogram so that its values represent probabilities. We divide each histogram bin count by the total number of pixels in the image $Npixels$. Suppose our image has $Npixels = 1000$ pixels. The normalized histogram would be:
* $Normalized Histogram (0) = 10/1000 = 0.01$
* $Normalized Histogram (1) = 20/1000 = 0.02$
* $Normalized Histogram (2) = 15/1000 = 0.015$
* ...
* $Normalized Histogram (100) = 5/1000 = 0.005$

Finally, we compute the cumulative sum of the normalized histogram values up to each intensity level $k$. This gives us the cumulative distribution function.
* $CDF (0) = Normalized Histogram(0) = 0.01$
* $CDF (1) = Normalized Histogram (0) + Normalized Histogram (1) = 0.01 + 0.02 = 0.03$
* $CDF (2) = Normalized Histogram (0) + Normalized Histogram (1) + Normalized Histogram (2) = 0.01 + 0.02 + 0.0015 = 0.045$
* ...
* $CDF (100) = 1$

Lastly we map the original intensity values of the image to their corresponding values in the CDF. This mapping redistributes the intensity values to achieve a more uniform distribution. 
Let $f(x)$ be the histogram equalization function that maps the original intensity value $x$ to their corresponding values in the CDF $CDF(x)$. This function can be fefine as:

$$f(x) = round(\frac{CDF(x) - min(CDF)}{max(CDF)-min(CDF)} * (L - 1))$$

Once the histogram equalization function is defined, you can update each pixel's intensity value in the image using this function

### Saturation

Firstly we iterate over the image and get the original pixels of that image so that we can perform the grayscale filter. We interpolate it with the saturation of the original image. 
Next we perform the interpolation by using the formula from the provided reference which is:

Alpha being ```alpha = ratio + 1```

$$(1 - alpha) * Sgrayscale + alpha * Soriginal$$

### White Balance

The first step in the Von Kries method is to convert the image from the RGB color space to the LMS color space. This conversion involves matrix multiplication with a transformation matrix. 

After converting to the LMS color space, you need to divide each pixel's LMS values by the LMS coordinates of the white point. Finally, the adjusted LMS values need to be converted back to the RGB color space to get the final white-balanced image.

Reference: <https://medium.com/@KuldeepDileep/chromatic-adaptation-with-matlab-code-9af2aaf9096a>

### Edges

To implement an edge detection filter using convolution first we set its kernel then we normalize it so that it doesn't change the overall brightness of the image. Then we invert each pixel in the image so that the image is clearer for visualization once we determine it's edges. 

$$pixel = 1 - pixel$$

Then we set the image as grayscaled so that we can find areas in the image where the intensity is higher which will allow us to work better with the image in calculating its edges.

Before convolving the image with the kernel, we create a new image and set only the edge pixels to white, which will allow us to work better with the image in calculating edges.

### Bilateral

A bilateral filter is a non-linear, edge-preserving, and noise-reducing smoothing filter for images. It replaces the intensity of each pixel with a weighted average of intensity values from nearby pixels. This weight can be based on a Gaussian distribution.

For the pixel $I(i,j)$ located at $(i,j)$, the bilateral filter weight of pixel $I(k,l)$, $I$ being the intensity value of the pixel, is given by:

$$w(i,j,k,l) = exp(- \frac{(i - k)^2 + (j - l)^2}{2σ_d^2} - \frac{||I(i,j) - I(k,l)||^2}{2σ_r^2})$$

After calculating the weights, normalize them: 

$$I_D(i,j) = \frac{Σ_{k, l} I(k,l)w(i,j,k,l)}{Σ_{k, l}w(i,j,k,l)}$$

where $I_D$ is the denoised intensity of pixel $(i, j)$

### FlyodFilter

The algorithm achieves dithering using error diffusion, meaning it pushes (adds) the residual quantization error of a pixel onto its neighboring pixels, to be dealt with later. It spreads the debt out according to the distribution.

#### Pseudocode

```jai
for each y from top to bottom do
    for each x from left to right do
        oldpixel := pixels[x][y]
        newpixel := find_closest_palette_color(oldpixel)
        pixels[x][y] := newpixel
        quant_error := oldpixel - newpixel
        pixels[x + 1][y    ] := pixels[x + 1][y    ] + quant_error × 7 / 16
        pixels[x - 1][y + 1] := pixels[x - 1][y + 1] + quant_error × 3 / 16
        pixels[x    ][y + 1] := pixels[x    ][y + 1] + quant_error × 5 / 16
        pixels[x + 1][y + 1] := pixels[x + 1][y + 1] + quant_error × 1 / 16
```

```jai
  find_closest_palette_color(oldpixel) = round(oldpixel / 255)
```

For more information: <https://en.wikipedia.org/wiki/Floyd%E2%80%93Steinberg_dithering>
