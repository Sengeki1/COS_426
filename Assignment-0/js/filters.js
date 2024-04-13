"use strict";

var Filters = Filters || {};

////////////////////////////////////////////////////////////////////////////////
// General utility functions
////////////////////////////////////////////////////////////////////////////////

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
  return ((val < min) ? min : ((val > max) ? max : val));
}

// extract vertex coordinates from a URL string
function stringToCoords( vertsString ) {
  var centers = [];
  var coordStrings = vertsString.split("x");
  var coordsSoFar = 0;
  for (var i = 0; i < coordStrings.length; i++) {
    var coords = coordStrings[i].split("y");
    var x = parseInt(coords[0]);
    var y = parseInt(coords[1]);
    if (!isNaN(x) && !isNaN(y)) {
      centers.push({x: x, y: y})
    }
  }

  return centers;
}

////////////////////////////////////////////////////////////////////////////////
// Filters
////////////////////////////////////////////////////////////////////////////////

// Fill the entire image with color
Filters.fillFilter = function( image, color ) {
  for (var x = 0; x < image.width; x++) {
    for (var y = 0; y < image.height; y++) {
      // uncomment this line to enable this function
      image.setPixel(x, y, color);
    }
  }
  return image;
};

// At each center, draw a solid circle with the specified radius and color
Filters.brushFilter = function( image, radius, color, vertsString ) {
  // centers is an array of (x, y) coordinates that each defines a circle center
  var centers = stringToCoords(vertsString);

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
};

/*
 * At each center, draw a soft circle with the specified radius and color.
 * Pixel opacity should linearly decrease with the radius from alpha_at_center to 0.
 */
Filters.softBrushFilter = function( image, radius, color, alpha_at_center, vertsString ) {
  // centers is an array of (x, y) coordinates that each defines a circle center
  var centers = stringToCoords(vertsString);

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
  //Gui.alertOnce ('softBrushFilter is not implemented yet');

  return image;
};

Filters.customFilter = function( image, value ) {
  // You can use this filter to do whatever you want, for example
  // trying out some new idea or implementing something for the
  // art contest.
  // Currently the 'value' argument will be 1 or whatever else you set
  // it to in the URL. You could use this value to switch between
  // a bunch of different versions of your code if you want to
  // code up a bunch of different things for the art contest.
  // ----------- STUDENT CODE BEGIN ------------
  // Gaussian Kernel
  var kernel = [[1, 2, 1],
                [2, 4, 2],
                [1, 2, 1]]

  var sum = 0
  kernel.forEach(row => {
    row.forEach(num => {
      sum += num
    })
  })
  kernel = kernel.map(row => row.map(num => num / sum)); // Normalizing the kernel

  var outputImage = image.copyImg() // Create a copy of the input image for output

  for (var x = 0; x < image.width; x++) {
    for (var y = 0; y < image.height; y++) {

      var accumulator = new Pixel(0, 0, 0)
      
      // Compute the bounds of the neighborhood
      var minX = Math.max(0, x - 1)
      var maxX = Math.min(image.width - 1, x + 1)
      var minY = Math.max(0, y - 1)
      var maxY = Math.min(image.height - 1, y + 1)

      for (var i = minX; i <= maxX; i++) {
        for (var j = minY; j <= maxY; j++) {

          var neighborhoodPixel = image.getPixel(i, j)
          var weightedPixel = neighborhoodPixel.multipliedBy(kernel[i - minX][j - minY])
          accumulator = accumulator.plus(weightedPixel)
        }
      }
      outputImage.setPixel(x, y, accumulator)
    }
  }
  // ----------- STUDENT CODE END ------------
  //Gui.alertOnce ('customFilter is not implemented yet');

  return outputImage;
};
