var Raytracer = Raytracer || {};
// ID
Raytracer.objectID = 0;
Raytracer.lightID = 0;

// material types - how to light the surface
Raytracer.BASICMATERIAL = 1;
Raytracer.PHONGMATERIAL = 2;
Raytracer.LAMBERTMATERIAL = 3;

// special materials
Raytracer.NONE = 0;
Raytracer.CHECKERBOARD = 1;
Raytracer.MYSPECIAL = 2;

// reflect types - how to bounce rays
Raytracer.NONEREFLECT = 1;
Raytracer.MIRRORREFLECT = 2;
Raytracer.GLASSREFLECT = 3;

Raytracer.curMaterial = {};

// for animation
Raytracer.frame = 0;
Raytracer.needsToDraw = true;

Raytracer.mouseDown = false;
Raytracer.lastMouseX = null;
Raytracer.lastMouseY = null;
//
Raytracer.camera = [0.15, -0.05, 0];

Raytracer.handleMouseDown = function(event) {
  Raytracer.mouseDown = true;
  Raytracer.lastMouseX = event.clientX;
  Raytracer.lastMouseY = event.clientY;
};

Raytracer.handleMouseUp = function(event) {
  Raytracer.mouseDown = false;
};

Raytracer.handleZoom = function(delta) {
  mat4.translate(Raytracer.RotationMatrix, [0.0, 0.0, 0.5 * delta]);
  Raytracer.needsToDraw = true;
};

Raytracer.handleMouseMove = function(event) {
  var newX = event.clientX;
  var newY = event.clientY;
  var deltaX = newX - Raytracer.lastMouseX;
  var deltaY = newY - Raytracer.lastMouseY;
  var moved = deltaX != 0 || deltaY != 0;

  if (!Raytracer.mouseDown || !moved) {
    return;
  }

  var degToRad = function(degrees) {
    return (degrees * Math.PI) / 180;
  };
  var newRotationMatrix = mat4.create();
  mat4.identity(newRotationMatrix);
  mat4.rotate(newRotationMatrix, degToRad(deltaX / 10), [0, 1, 0]);
  mat4.rotate(newRotationMatrix, degToRad(deltaY / 10), [1, 0, 0]);
  mat4.multiply(newRotationMatrix, Raytracer.RotationMatrix, Raytracer.RotationMatrix);

  Raytracer.lastMouseX = newX;
  Raytracer.lastMouseY = newY;
  Raytracer.needsToDraw = true;
};

Raytracer.initShader = function(program, shaderType, src, debug) {
  var shader = this.gl.createShader(shaderType);
  this.gl.shaderSource(shader, src);
  this.gl.compileShader(shader);

  // check compile status and report error
  var ok = this.gl.getShaderParameter(shader, this.gl.COMPILE_STATUS);

  if (debug || !ok) {
    var log = this.gl.getShaderInfoLog(shader);
    var msg = debug ? "Debug status of " : "Compile error in ";
    msg += "shader type " + shaderType + ": " + log;
    alert(msg);
    console.log(msg);
  }
  this.gl.attachShader(program, shader);
  return shader;
};

Raytracer.init = function(height, width, debug) {
  canvas = document.getElementById("canvas");
  this.gl = canvas.getContext("experimental-webgl", { preserveDrawingBuffer: true });
  canvas.width = width;
  canvas.height = height;

  this.gl.viewportWidth = canvas.width;
  this.gl.viewportHeight = canvas.height;

  this.gl.viewport(0, 0, this.gl.drawingBufferWidth, this.gl.drawingBufferHeight);

  this.gl.clearColor(0.0, 0.0, 0.0, 1.0);
  this.gl.clear(this.gl.COLOR_BUFFER_BIT);

  var fSrcBase = Parser.parseTxt("shaders/fragmentShader.frag");
  var vSrc = Parser.parseTxt("shaders/vertexShader.vert");
  var fSrc = fSrcBase + Scene.getIntersectFunction();

  this.program = this.gl.createProgram();

  var compileStartTime = performance.now();

  this.initShader(this.program, this.gl.VERTEX_SHADER, vSrc, debug);
  this.initShader(this.program, this.gl.FRAGMENT_SHADER, fSrc, debug);

  var compileTime = Math.round(performance.now() - compileStartTime);
  // console.log('shader compilation completed in ' + compileTime + ' ms.');

  this.gl.linkProgram(this.program);
  this.gl.useProgram(this.program);

  this.gl.uniform1f(this.gl.getUniformLocation(this.program, "width"), width);
  this.gl.uniform1f(this.gl.getUniformLocation(this.program, "height"), height);

  positionLocation = this.gl.getAttribLocation(this.program, "a_position");
  this.gl.enableVertexAttribArray(positionLocation);

  var bufferGeom = new Float32Array([
    -1.0,
    -1.0,
    1.0,
    -1.0,
    -1.0,
    1.0,
    -1.0,
    1.0,
    1.0,
    -1.0,
    1.0,
    1.0,
  ]);
  var buffer = this.gl.createBuffer();
  buffer.itemSize = 2;
  buffer.numItems = 6;

  this.gl.bindBuffer(this.gl.ARRAY_BUFFER, buffer);
  this.gl.bufferData(this.gl.ARRAY_BUFFER, bufferGeom, this.gl.STATIC_DRAW);
  this.gl.vertexAttribPointer(positionLocation, 2, this.gl.FLOAT, false, 0, 0);

  this.gl.uniform1i(this.gl.getUniformLocation(this.program, "frame"), this.frame);

  Raytracer.RotationMatrix = mat4.create();
  mat4.identity(Raytracer.RotationMatrix);

  canvas.onmousedown = Raytracer.handleMouseDown;
  document.onmouseup = Raytracer.handleMouseUp;
  document.onmousemove = Raytracer.handleMouseMove;
};

Raytracer.setCamera = function(cameraAngle) {
  //rotation matrix
  var newRotationMatrix = mat4.create();
  mat4.identity(newRotationMatrix);
  mat4.rotate(newRotationMatrix, cameraAngle[0] * Math.PI, [1, 0, 0]);
  mat4.rotate(newRotationMatrix, cameraAngle[1] * Math.PI, [0, 1, 0]);
  mat4.rotate(newRotationMatrix, cameraAngle[2] * Math.PI, [0, 0, 1]);
  console.log(cameraAngle);
  mat4.multiply(newRotationMatrix, Raytracer.RotationMatrix, Raytracer.RotationMatrix);
};
Raytracer.addLight = function(px, py, pz, cr, cg, cb, intensity, attenuate) {
  var lightID = "lights[" + this.lightID + "].";
  this.setUniform("3f", lightID + "position", px, py, pz);
  this.setUniform("3f", lightID + "color", cr, cg, cb);
  this.setUniform("1f", lightID + "intensity", intensity);
  this.setUniform("1f", lightID + "attenuate", attenuate);
  this.lightID++;
};

Raytracer.setUniform = function(varType, varName, v0, v1, v2) {
  var unifName = "uniform" + varType;
  v0 = v0 || 0.0;
  v1 = v1 || 0.0;
  v2 = v2 || 0.0;
  // v1 and v2 may be undefined because we only need 1 argument, but this is ok
  this.gl[unifName](this.gl.getUniformLocation(this.program, varName), v0, v1, v2);
};

Raytracer.render = function(animated) {
  this.frame++;
  if (animated) {
    this.setUniform("1i", "frame", this.frame);
  }
  //rotation matrix

  this.setUniform("Matrix4fv", "uMVMatrix", false, this.RotationMatrix);

  if (this.needsToDraw || animated) {
    this.gl.drawArrays(this.gl.TRIANGLES, 0, 6);
    this.needsToDraw = false;
  }
};
