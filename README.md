# Paint

A simple paint program with the purpose to ease you into both JavaScript programming and computational photography as 
comfortably.
The implementation done on this project wore made on the ```filter.js``` file, to create features.

## Assignment Features

### Brush

Since the program already captures the click and computes the mouse position in image-coordinates the only thing left us to do is to fill-out the pixels in a circle of radius ```r``` that is centered at some point ```(x, y)```.
Given the mouse position in image-coordinates we already accomplish the first step which is to get the pixel in the screen.

To make a solid circle with a specified radius and color first we have to iterate for each point clicked in the image.
```js
var centers = stringToCoords(vertsString);
  for (var i = 0; i < centers.length; i++) {
    var center = centers[i]
```

After iterating over each center points we compute the bouding box of a circle given the radius that was provided in the function argument.
You may ask yourself, what is a bouding box? 
A bouding box is the box (square) in which a object is inscribed, in this case a circle.

![image](https://github.com/Sengeki1/JS_Paint_COS426/assets/106749775/59acaa35-c38b-4839-b44b-28cd8ed14e31)

Notice that the radius of the circle is exactly half the length of a side of the square.
This allows us to compute the points in the circumference in a more efficient way.
To calculate the points in the circumference assuming a standard coordinate system (```y - values``` increase in the up direction and ```x - values``` increase in the right direction)
we iterate over each of the corners of the box and well-done you made up a solid square.

Now to make up a circle we need to check the distance between $(x_c, y_c)$ and $(x_p, y_p)$ which is given by the Pythagorean theorem as
$$d = \sqrt{(x_p - x_c)^2 + (y_p - y_c)^2}$$
The point  $(x_p, y_p)$ is inside the circle if $d < r$, on the circle if $d = r$, and outside the circle if $d > r$.
```js
  var dx = x - center.x
  var dy = y - center.y 

  var d = Math.sqrt((dx * dx) + (dy * dy))

  if (d < radius) {
    image.setPixel(x, y, color)
  }
```

### Soft Brush

Now to make a soft brush, we draw a colored circle with opacity decreasing from center to the edge. We already finished the first step which was to create a brush now we only have to add opacity.
To add opacity is simple. We already know that a pixel is the combination of a color system (RGB) and alpha value. Now to change that alpha value we only have to set that color channel into the value passed in the function arguments
```js
  color.a = alpha_at_center
  image.setPixel(x, y, color)
```
Now to decrease it linearly along the radius and becomes zero at the edge of the circle, mathematically we calculate the width of the image (number of pixels from left to right) and the height of the image (number of pixels from top to bottom), ```x``` be the horizontal position of the pixel (starting from 0 at the left), and ```y``` be the vertical position of the pixel (starting from 0 at the top).

First, we find the coordinates of the center of the circle $(\frac{W}{2}, \frac{H}{2})$. Calculated the distance of each pixel we only have to scale the distance to the ```range[0,1]``` and subtract it from 1 to get the alpha value.
$$alpha(x, y) = 1 -\frac{distance}{\sqrt{(\frac{W}{2})^2 + (\frac{H}{2})^2}}$$
```js
  var Height = radius * 2
  var Width = radius * 2
  color.a = alpha_at_center - distance / (Math.sqrt((Width / 2) * (Width / 2) + (Height / 2) * (Height / 2)))
  image.setPixel(x, y, color)
```

### Gaussian Blur
