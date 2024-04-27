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
