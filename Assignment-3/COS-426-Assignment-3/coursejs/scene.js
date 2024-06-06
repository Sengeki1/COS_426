var Scene = Scene || {};

Scene.defaultMaterial = {
  type: "LAMBERTMATERIAL",
  color: [0.5, 0.5, 0.5],
  shininess: 10000,
  specular: [0.0, 0.0, 0.0],
  reflectType: "NONEREFLECT",
  reflectivity: 0,
  refractionRatio: 0.8,
  special: "NONE",
};

Scene.setUniforms = function() {
  var sceneData = Parser.parseJson("scenes/" + Scene.sceneName + ".json");

  if (sceneData.camera !== undefined) {
    Raytracer.setCamera(sceneData.camera);
  }

  for (var i = 0; i < sceneData.lights.length; i++) {
    var l = sceneData.lights[i];
    Raytracer.addLight(
      l.pos[0],
      l.pos[1],
      l.pos[2],
      l.color[0],
      l.color[1],
      l.color[2],
      l.intensity,
      l.attenuate
    );
  }
  Raytracer.setUniform("1i", "numLights", sceneData.lights.length);
  Raytracer.setUniform("3f", "camera", 0.0, 0.0, -25.0);
};

var arrayToVec3String = function(array) {
  return "vec3(" + array[0] + ", " + array[1] + ", " + array[2] + ")";
};
var fixPrecision = function(a, precision) {
  var base = Math.pow(10, precision);
  return Math.round(a * base) / base;
};
// construct the glsl code for rayIntersectScene with the scene hardcoded in
// (this is a performance/compatibility optimization)
Scene.getIntersectFunction = function() {
  var MAX_TRIANGLE = 50;

  var sceneData = Parser.parseJson("scenes/" + Scene.sceneName + ".json");
  var count_triangle = 0;

  var matStr = "";

  var funcStr = "\n\n";
  funcStr +=
    "float rayIntersectScene( Ray ray, out Material out_mat, out Intersection out_intersect ) {\n" +
    "Intersection intersect;\n" +
    "float cur_dist;\n" +
    "float dist = INFINITY;\n";

  for (var i = 0; i < sceneData.objects.length; i++) {
    var obj = sceneData.objects[i];
    var m = obj.material;
    for (var prop in Scene.defaultMaterial) {
      if (m[prop] === undefined) {
        m[prop] = Scene.defaultMaterial[prop];
      }
    }
    matStr +=
      "Material material" +
      i +
      " = Material(" +
      Raytracer[m.type] +
      ", " +
      arrayToVec3String(m.color) +
      ", " +
      m.shininess.toFixed(2) +
      ", " +
      arrayToVec3String(m.specular) +
      ", " +
      Raytracer[m.reflectType] +
      ", " +
      m.reflectivity.toFixed(2) +
      ", " +
      m.refractionRatio.toFixed(2) +
      ", " +
      Raytracer[m.special] +
      ");\n";

    switch (obj.type) {
      case "sphere":
        funcStr +=
          "cur_dist = findIntersectionWithSphere(ray, " +
          arrayToVec3String(obj.center) +
          ", " +
          obj.radius.toFixed(2) +
          ", " +
          "intersect);\n";

        break;

      case "plane":
        funcStr +=
          "cur_dist = findIntersectionWithPlane(ray, " +
          arrayToVec3String(obj.normal) +
          ", " +
          obj.dist.toFixed(2) +
          ", " +
          "intersect);\n";
        break;

      case "cylinder":
        funcStr +=
          "cur_dist = findIntersectionWithCylinder(ray, " +
          arrayToVec3String(obj.bottomCenter) +
          ", " +
          arrayToVec3String(obj.topCenter) +
          ", " +
          obj.radius.toFixed(2) +
          ", " +
          "intersect);\n";
        break;

      case "cone":
        funcStr +=
          "cur_dist = findIntersectionWithCone(ray, " +
          arrayToVec3String(obj.bottomCenter) +
          ", " +
          arrayToVec3String(obj.topCenter) +
          ", " +
          obj.radius.toFixed(2) +
          ", " +
          "intersect);\n";
        break;

      case "box":
        funcStr +=
          "cur_dist = findIntersectionWithBox(ray, " +
          arrayToVec3String(obj.minCorner) +
          ", " +
          arrayToVec3String(obj.maxCorner) +
          ", " +
          "intersect);\n";
        break;

      case "triangle":
        funcStr +=
          "cur_dist = findIntersectionWithTriangle(ray, " +
          arrayToVec3String(obj.t1) +
          ", " +
          arrayToVec3String(obj.t2) +
          ", " +
          arrayToVec3String(obj.t3) +
          ", " +
          "intersect);\n";
        break;

      case "mesh":
        var manager = new THREE.LoadingManager();
        var loader = new OBJLoader(manager);
        var container = loader.parse(Parser.parseTxt(obj.objfile));

        var vertices = container[0];
        var faces = container[1];

        for (var k = 0; k < vertices.length; k++) {
          var v = vertices[k];

          if (obj.scale !== undefined) {
            v.x = v.x * obj.scale;
            v.y = v.y * obj.scale;
            v.z = v.z * obj.scale;
          }
          if (obj.offset !== undefined) {
            v.x = v.x + obj.offset[0];
            v.y = v.y + obj.offset[1];
            v.z = v.z + obj.offset[2];
          }

          vertices[k] = [fixPrecision(v.x, 2), fixPrecision(v.y, 2), fixPrecision(v.z, 2)];
        }
        for (var k = 0; k < faces.length; k++) {
          var cur_face = faces[k];
          for (var v = 1; v < cur_face.length - 1; v++) {
            if (++count_triangle > MAX_TRIANGLE) {
              break;
            }
            // triangulate to 0, v, v+1
            funcStr +=
              "cur_dist = findIntersectionWithTriangle(ray, " +
              arrayToVec3String(vertices[cur_face[0]]) +
              ", " +
              arrayToVec3String(vertices[cur_face[v]]) +
              ", " +
              arrayToVec3String(vertices[cur_face[v + 1]]) +
              ", " +
              "intersect);\n";
            funcStr +=
              "if ( chooseCloserIntersection(cur_dist, dist, intersect, out_intersect) ) {\n" +
              "out_mat = material" +
              i +
              ";\n" +
              "}\n";
          }
        }
        break;

      default:
        break;
    }
    funcStr +=
      "if ( chooseCloserIntersection(cur_dist, dist, intersect, out_intersect) ) {\n" +
      "out_mat = material" +
      i +
      ";\n" +
      "}\n";
  }
  funcStr += "return dist;\n } \n";

  return matStr + "\n" + funcStr;
};
