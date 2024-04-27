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
