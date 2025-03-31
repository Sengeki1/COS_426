// set the precision of the float values (necessary if using float)
#ifdef GL_FRAGMENT_PRECISION_HIGH
precision highp float;
#else
precision mediump float;
#endif
precision mediump int;

// flag for using soft shadows
#define SOFT_SHADOWS 1

// define number of soft shadow samples to take
#define SOFT_SAMPLING 3

// define constant parameters
// EPS is for the precision issue
#define INFINITY 1.0e+12
#define EPS 1.0e-3

// define maximum recursion depth for rays
#define MAX_RECURSION 8

// define constants for scene setting
#define MAX_LIGHTS 10

// define texture types
#define NONE 0
#define CHECKERBOARD 1
#define MYSPECIAL 2

// define material types
#define BASICMATERIAL 1
#define PHONGMATERIAL 2
#define LAMBERTMATERIAL 3

// define reflect types - how to bounce rays
#define NONEREFLECT 1
#define MIRRORREFLECT 2
#define GLASSREFLECT 3

struct Shape {
  int shapeType;
  vec3 v1;
  vec3 v2;
  float rad;
};

struct Material {
  int materialType;
  vec3 color;
  float shininess;
  vec3 specular;

  int materialReflectType;
  float reflectivity;
  float refractionRatio;
  int special;
};

struct Object {
  Shape shape;
  Material material;
};

struct Light {
  vec3 position;
  vec3 color;
  float intensity;
  float attenuate;
};

struct Ray {
  vec3 origin;
  vec3 direction;
};

struct Intersection {
  vec3 position;
  vec3 normal;
};

// uniform
uniform mat4 uMVMatrix;
uniform int frame;
uniform float height;
uniform float width;
uniform vec3 camera;
uniform int numObjects;
uniform int numLights;
uniform Light lights[MAX_LIGHTS];
uniform vec3 objectNorm;

// varying
varying vec2 v_position;

// find then position some distance along a ray
vec3 rayGetOffset(Ray ray, float dist) {
  return ray.origin + (dist * ray.direction);
}

// if a newly found intersection is closer than the best found so far, record
// the new intersection and return true; otherwise leave the best as it was and
// return false.
bool chooseCloserIntersection(float dist, inout float best_dist,
                              inout Intersection intersect,
                              inout Intersection best_intersect) {
  if (best_dist <= dist)
    return false;
  best_dist = dist;
  best_intersect.position = intersect.position;
  best_intersect.normal = intersect.normal;
  return true;
}

// put any general convenience functions you want up here
// ----------- STUDENT CODE BEGIN ------------

// Function to compute the dot product between two vectors
float dotProduct(vec3 a, vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Function to compute the length of a vector
float length(vec3 v) {
    return sqrt(dotProduct(v, v));
}

// Function to normalize a vector
vec3 normalize(vec3 v) {
    float len = length(v);
    return vec3(v.x / len, v.y / len, v.z / len);
}

// Function to reflect a vector around a normal
vec3 reflect(vec3 incident, vec3 normal) {
    return incident - 2.0 * dot(incident, normal) * normal;
}

// Function to mix two colors
vec3 mixColors(vec3 color1, vec3 color2, float ratio) {
    return color1 * ratio + color2 * (1.0 - ratio);
}

// ----------- STUDENT CODE END ------------


// forward declaration
float rayIntersectScene(Ray ray, out Material out_mat,
                        out Intersection out_intersect);

// Plane
// this function can be used for plane, triangle, and box
float findIntersectionWithPlane(Ray ray, vec3 norm, float dist,
                                out Intersection intersect) {
  float a = dot(ray.direction, norm);
  float b = dot(ray.origin, norm) - dist;

  if (a < EPS && a > -EPS)
    return INFINITY;

  float len = -b / a;
  if (len < EPS)
    return INFINITY;

  intersect.position = rayGetOffset(ray, len);
  intersect.normal = norm;
  return len;
}

// Triangle
float findIntersectionWithTriangle(Ray ray, vec3 t1, vec3 t2, vec3 t3,
                                   out Intersection intersect) {
    vec3 edge1 = t2 - t1;
    vec3 edge2 = t3 - t1;
    vec3 pvec = cross(ray.direction, edge2);
    float det = dot(edge1, pvec);

    // Ray and triangle are parallel if det is close to 0
    if (abs(det) < EPS) {
        return INFINITY;
    }

    float invDet = 1.0 / det;
    vec3 tvec = ray.origin - t1;
    float u = dot(tvec, pvec) * invDet;

    // u is outside the range [0, 1], so no intersection
    if (u < 0.0 || u > 1.0) {
        return INFINITY;
    }

    vec3 qvec = cross(tvec, edge1);
    float v = dot(ray.direction, qvec) * invDet;

    // v is outside the range [0, 1], so no intersection
    if (v < 0.0 || u + v > 1.0) {
        return INFINITY;
    }

    float t = dot(edge2, qvec) * invDet;

    if (t > EPS) {
        intersect.position = rayGetOffset(ray, t);
        intersect.normal = normalize(cross(edge1, edge2));
        return t;
    }

    return INFINITY;
}

// Sphere
float findIntersectionWithSphere(Ray ray, vec3 center, float radius,
                                 out Intersection intersect) {
    vec3 oc = ray.origin - center;
    float a = dot(ray.direction, ray.direction);
    float b = 2.0 * dot(oc, ray.direction);
    float c = dot(oc, oc) - radius * radius;
    float discriminant = b * b - 4.0 * a * c;

    // Se o discriminante for negativo, não há interseção
    if (discriminant < 0.0) {
        return INFINITY;
    }

    float sqrtDiscriminant = sqrt(discriminant);
    float t1 = (-b - sqrtDiscriminant) / (2.0 * a);
    float t2 = (-b + sqrtDiscriminant) / (2.0 * a);

    // Verificar se as interseções estão na direção do raio
    if (t1 > EPS || t2 > EPS) {
        float t = min(t1, t2);
        intersect.position = rayGetOffset(ray, t);
        intersect.normal = normalize(intersect.position - center);
        return t;
    }

    return INFINITY;
}

// Box
float findIntersectionWithBox(Ray ray, vec3 pmin, vec3 pmax,
                              out Intersection out_intersect) {
    vec3 tmin = (pmin - ray.origin) / ray.direction;
    vec3 tmax = (pmax - ray.origin) / ray.direction;

    vec3 t1 = min(tmin, tmax);
    vec3 t2 = max(tmin, tmax);

    float tmin_max = max(max(t1.x, t1.y), t1.z);
    float tmax_min = min(min(t2.x, t2.y), t2.z);

    // Se o intervalo mínimo máximo for maior que o intervalo máximo mínimo, não há interseção
    if (tmin_max > tmax_min) {
        return INFINITY;
    }

    // Verificar se o intervalo mínimo máximo é maior que zero
    float t = (tmin_max > 0.0) ? tmin_max : tmax_min;

    if (t < EPS) {
        return INFINITY;
    }

    vec3 intersection_point = ray.origin + t * ray.direction;

    // Verificar se o ponto de interseção está dentro da caixa
    if (intersection_point.x < pmin.x - EPS || intersection_point.x > pmax.x + EPS ||
        intersection_point.y < pmin.y - EPS || intersection_point.y > pmax.y + EPS ||
        intersection_point.z < pmin.z - EPS || intersection_point.z > pmax.z + EPS) {
        return INFINITY;
    }

    out_intersect.position = intersection_point;

    // Calcular a normal da face mais próxima
    vec3 d = normalize(ray.direction);
    vec3 normal = -sign(d) * step(t1, t2);
    out_intersect.normal = normal;

    return t;
}

// Cylinder
float getIntersectOpenCylinder(Ray ray, vec3 center, vec3 axis, float len,
                               float rad, out Intersection intersect) {
    // Transformar o raio para o sistema de coordenadas do cilindro
    vec3 oc = ray.origin - center;
    vec3 oc_axis = cross(oc, axis);
    vec3 d_axis = cross(ray.direction, axis);

    float a = dot(d_axis, d_axis);
    float b = 2.0 * dot(d_axis, oc_axis);
    float c = dot(oc_axis, oc_axis) - rad * rad;

    // Resolver a equação quadrática
    float discriminant = b * b - 4.0 * a * c;

    if (discriminant < 0.0) {
        return INFINITY; // Não há interseção
    }

    float sqrt_discriminant = sqrt(discriminant);
    float t1 = (-b - sqrt_discriminant) / (2.0 * a);
    float t2 = (-b + sqrt_discriminant) / (2.0 * a);

    // Verificar se a interseção está dentro do cilindro
    vec3 p1 = ray.origin + t1 * ray.direction;
    vec3 p2 = ray.origin + t2 * ray.direction;

    float t = INFINITY;
    vec3 intersection_point;

    if (t1 >= 0.0 && t1 <= len && p1.y >= center.y && p1.y <= center.y + len) {
        t = t1;
        intersection_point = p1;
    } else if (t2 >= 0.0 && t2 <= len && p2.y >= center.y && p2.y <= center.y + len) {
        t = t2;
        intersection_point = p2;
    } else {
        return INFINITY; // Não há interseção dentro do cilindro
    }

    intersect.position = intersection_point;

    // Calcular a normal no ponto de interseção
    vec3 normal = normalize(intersection_point - center - axis * dot(intersection_point - center, axis));
    intersect.normal = normal;

    return t;
}

// Disc
float getIntersectDisc(Ray ray, vec3 center, vec3 norm, float rad,
                       out Intersection intersect) {
    // Verificar se o raio é paralelo ao plano do disco
    float denom = dot(norm, ray.direction);
    if (abs(denom) < EPS) {
        return INFINITY; // Raio é paralelo ao plano do disco
    }

    // Calcular o parâmetro de distância do ponto de interseção ao longo do raio
    float t = dot(center - ray.origin, norm) / denom;

    // Verificar se a interseção está dentro do disco
    vec3 intersection_point = ray.origin + t * ray.direction;
    vec3 v = intersection_point - center;
    float dist_squared = dot(v, v);
    if (dist_squared > rad * rad) {
        return INFINITY; // Interseção está fora do disco
    }

    intersect.position = intersection_point;
    intersect.normal = norm;

    return t;
}


float findIntersectionWithCylinder(Ray ray, vec3 center, vec3 apex,
                                   float radius,
                                   out Intersection out_intersect) {
  vec3 axis = apex - center;
  float len = length(axis);
  axis = normalize(axis);

  Intersection intersect;
  float best_dist = INFINITY;
  float dist;

  // -- infinite cylinder
  dist = getIntersectOpenCylinder(ray, center, axis, len, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);

  // -- two caps
  dist = getIntersectDisc(ray, center, -axis, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);
  dist = getIntersectDisc(ray, apex, axis, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);
  return best_dist;
}

// Cone
float getIntersectOpenCone(Ray ray, vec3 apex, vec3 axis, float len,
                           float radius, out Intersection intersect) {
    // Vetor do ponto de origem do raio ao vértice do cone
    vec3 O = ray.origin - apex;

    // Vetor de direção do raio
    vec3 D = ray.direction;

    // Computar o ângulo ao quadrado
    float cos2_theta = radius * radius / (radius * radius + len * len);

    // Computar o vetor direção do cone
    vec3 axis_dir = normalize(axis);

    // Computar o vetor D' = D - (D dot axis_dir) * axis_dir
    vec3 D_prime = D - dot(D, axis_dir) * axis_dir;

    // Computar o vetor O' = O - (O dot axis_dir) * axis_dir
    vec3 O_prime = O - dot(O, axis_dir) * axis_dir;

    // Calcular as partes quadráticas da equação do cone
    float A = dot(D_prime, D_prime) - cos2_theta * dot(D, axis_dir) * dot(D, axis_dir);
    float B = 2.0 * (dot(D_prime, O_prime) - cos2_theta * dot(D, axis_dir) * dot(O, axis_dir));
    float C = dot(O_prime, O_prime) - cos2_theta * dot(O, axis_dir) * dot(O, axis_dir);

    // Resolver a equação quadrática
    float discriminant = B * B - 4.0 * A * C;
    if (discriminant < 0.0) {
        return INFINITY; // Não há interseção com o cone
    }

    float t1 = (-B - sqrt(discriminant)) / (2.0 * A);
    float t2 = (-B + sqrt(discriminant)) / (2.0 * A);

    // Verificar se o ponto de interseção está dentro do cone
    vec3 intersection_point = ray.origin + t1 * ray.direction;
    vec3 AO = intersection_point - apex;
    float dot_product = dot(AO, axis_dir);
    if (dot_product < 0.0 || dot_product > len) {
        return INFINITY; // Ponto de interseção está fora do cone
    }

    intersect.position = intersection_point;
    intersect.normal = normalize(AO - dot_product * axis_dir);

    return t1;
}

float findIntersectionWithCone(Ray ray, vec3 center, vec3 apex, float radius,
                               out Intersection out_intersect) {
  vec3 axis = center - apex;
  float len = length(axis);
  axis = normalize(axis);

  // -- infinite cone
  Intersection intersect;
  float best_dist = INFINITY;
  float dist;

  // -- infinite cone
  dist = getIntersectOpenCone(ray, apex, axis, len, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);

  // -- caps
  dist = getIntersectDisc(ray, center, axis, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);

  return best_dist;
}

vec3 calculateSpecialDiffuseColor(Material mat, vec3 posIntersection, vec3 normalVector) {
    if (mat.special == CHECKERBOARD) {
        // Calcular uma textura de tabuleiro de xadrez
        int check = int(mod(floor(posIntersection.x) + floor(posIntersection.z), 2.0));
        if (check == 0) {
            return vec3(0.0); // Cor preta
        } else {
            return vec3(1.0); // Cor branca
        }
    } else if (mat.special == MYSPECIAL) {
        // Calcular a cor especial MYSPECIAL
        return vec3(0.5, 0.5, 0.0); // Amarelo
    }

    // Se não for um material especial, retorne a cor do material
    return mat.color;
}

vec3 calculateDiffuseColor(Material mat, vec3 posIntersection,
                           vec3 normalVector) {
  // Special colors
  if (mat.special != NONE) {
    return calculateSpecialDiffuseColor(mat, posIntersection, normalVector);
  }
  return vec3(mat.color);
}

// check if position pos in in shadow with respect to a particular light.
// lightVec is the vector from that position to that light -- it is not
// normalized, so its length is the distance from the position to the light
bool pointInShadow(vec3 pos, vec3 lightVec) {
    // Direção do raio para a luz
    Ray shadowRay;
    shadowRay.origin = pos;
    shadowRay.direction = normalize(lightVec);

    // Trace o raio para verificar se ele atinge algum objeto no caminho para a luz
    float shadowDist = rayIntersectScene(shadowRay, out Material tempMat, out Intersection tempIntersect);
    
    // Se a distância até a interseção for menor que a distância até a luz, a posição está na sombra
    return shadowDist < length(lightVec);
}

// use random sampling to compute a ratio that represents the
// fractional contribution of the light to the position pos.
// lightVec is the vector from that position to that light -- it is not
// normalized, so its length is the distance from the position to the light
float softShadowRatio(vec3 pos, vec3 lightVec) {
    float totalSamples = 0.0;
    float shadowedSamples = 0.0;
    
    // Realize amostragem aleatória para calcular a razão de sombra suave
    for (int i = 0; i < SOFT_SAMPLING; i++) {
        // Gere uma posição de amostra aleatória dentro de uma área pequena em torno da posição
        vec3 samplePos = pos + (normalize(lightVec) * (rand(float(i)) - 0.5));
        
        // Verifique se a posição de amostra está na sombra
        if (pointInShadow(samplePos, lightVec)) {
            shadowedSamples += 1.0;
        }
        totalSamples += 1.0;
    }
    
    // Calcule a razão de sombra suave como a fração de amostras na sombra
    return shadowedSamples / totalSamples;
}

vec3 getLightContribution(Light light, Material mat, vec3 posIntersection,
                          vec3 normalVector, vec3 eyeVector, bool phongOnly,
                          vec3 diffuseColor) {
  vec3 lightVector = light.position - posIntersection;


  float ratio = 1.0; // default to 1.0 for hard shadows
  if (SOFT_SHADOWS == 1) {
    // if using soft shadows, call softShadowRatio to determine
    // fractional light contribution
    ratio = softShadowRatio(posIntersection, lightVector);
  }
  else {
    // check if point is in shadow with light vector
    if (pointInShadow(posIntersection, lightVector)) {
      return vec3(0.0, 0.0, 0.0);
    }
  }

  // Slight optimization for soft shadows
  if (ratio < EPS) {
    return vec3(0.0, 0.0, 0.0);
  }


  // normalize the light vector for the computations below
  float distToLight = length(lightVector);
  lightVector /= distToLight;

  if (mat.materialType == PHONGMATERIAL ||
      mat.materialType == LAMBERTMATERIAL) {
    vec3 contribution = vec3(0.0, 0.0, 0.0);

    // get light attenuation
    float attenuation = light.attenuate * distToLight;
    float diffuseIntensity =
        max(0.0, dot(normalVector, lightVector)) * light.intensity;

    // glass and mirror objects have specular highlights but no diffuse lighting
    if (!phongOnly) {
      contribution +=
          diffuseColor * diffuseIntensity * light.color / attenuation;
    }

    if (mat.materialType == PHONGMATERIAL) {
      // Start with just black by default (i.e. no Phong term contribution)
      vec3 phongTerm = vec3(0.0, 0.0, 0.0);
      // ----------- STUDENT CODE BEGIN ------------
      vec3 h = normalize(lightVector + eyeVector);
      float specularIntensity = pow(max(0.0, dot(normalVector, h)), mat.shininess);
      phongTerm = mat.specular * specularIntensity * light.color / attenuation;
      // ----------- Our reference solution uses 4 lines of code.
      // ----------- STUDENT CODE END ------------
      contribution += phongTerm;
    }

    return ratio * contribution;
  } else {
    return ratio * diffuseColor;
  }
}

vec3 calculateColor(Material mat, vec3 posIntersection, vec3 normalVector,
                    vec3 eyeVector, bool phongOnly) {
  // The diffuse color of the material at the point of intersection
  // Needed to compute the color when accounting for the lights in the scene
  vec3 diffuseColor = calculateDiffuseColor(mat, posIntersection, normalVector);

  // color defaults to black when there are no lights
  vec3 outputColor = vec3(0.0, 0.0, 0.0);

  // Loop over the MAX_LIGHTS different lights, taking care not to exceed
  // numLights (GLSL restriction), and accumulate each light's contribution
  // to the point of intersection in the scene.
  // ----------- STUDENT CODE BEGIN ------------
  for (int i = 0; i < numLights; i++) {
    // Get the contribution of the current light to the point of intersection
    vec3 lightContribution = getLightContribution(lights[i], mat, posIntersection, normalVector, eyeVector, phongOnly, diffuseColor);
    // Accumulate the contribution to the output color
    outputColor += lightContribution;
  }
  // ----------- Our reference solution uses 9 lines of code.
  // Return diffuseColor by default, so you can see something for now.
  return diffuseColor;
  // ----------- STUDENT CODE END ------------
}

// find reflection or refraction direction (depending on material type)
vec3 calcReflectionVector(Material material, vec3 direction, vec3 normalVector,
                          bool isInsideObj) {
  if (material.materialReflectType == MIRRORREFLECT) {
    return reflect(direction, normalVector);
  }
  // If it's not mirror, then it is a refractive material like glass.
  // Compute the refraction direction.
  // See lecture 13 slide (lighting) on Snell's law.
  // The eta below is eta_i/eta_r.
  // ----------- STUDENT CODE BEGIN ------------
  float eta =
      (isInsideObj) ? 1.0 / material.refractionRatio : material.refractionRatio;
  vec3 refractDir = refract(direction, normalVector, eta);
  // If the refraction is total internal reflection, return reflection
  if (refractDir == vec3(0.0, 0.0, 0.0)) {
    return reflect(direction, normalVector);
  }
  return refractDir;
  // ----------- STUDENT CODE END ------------
}

vec3 traceRay(Ray ray) {
  // Accumulate the final color from tracing this ray into resColor.
  vec3 resColor = vec3(0.0, 0.0, 0.0);

  // Accumulate a weight from tracing this ray through different materials
  // based on their BRDFs. Initially all 1.0s (i.e. scales the initial ray's
  // RGB color by 1.0 across all color channels). This captures the BRDFs
  // of the materials intersected by the ray's journey through the scene.
  vec3 resWeight = vec3(1.0, 1.0, 1.0);

  // Flag for whether the ray is currently inside of an object.
  bool isInsideObj = false;

  // Iteratively trace the ray through the scene up to MAX_RECURSION bounces.
  for (int depth = 0; depth < MAX_RECURSION; depth++) {
    // Fire the ray into the scene and find an intersection, if one exists.
    //
    // To do so, trace the ray using the rayIntersectScene function, which
    // also accepts a Material struct and an Intersection struct to store
    // information about the point of intersection. The function returns
    // a distance of how far the ray travelled before it intersected an object.
    //
    // Then, check whether or not the ray actually intersected with the scene.
    // A ray does not intersect the scene if it intersects at a distance
    // "equal to zero" or far beyond the bounds of the scene. If so, break
    // the loop and do not trace the ray any further.
    // (Hint: You should probably use EPS and INFINITY.)
    // ----------- STUDENT CODE BEGIN ------------
    float intersectionDistance = rayIntersectScene(ray, hitMaterial, intersect);
    if (intersectionDistance == INFINITY || intersectionDistance <= EPS) {
      break;
    }
    // ----------- Our reference solution uses 4 lines of code.
    // ----------- STUDENT CODE END ------------

    // Compute the vector from the ray towards the intersection.
    vec3 posIntersection = intersect.position;
    vec3 normalVector    = intersect.normal;

    vec3 eyeVector = normalize(ray.origin - posIntersection);

    // Determine whether we are inside an object using the dot product
    // with the intersection's normal vector
    if (dot(eyeVector, normalVector) < 0.0) {
        normalVector = -normalVector;
        isInsideObj = true;
    } else {
        isInsideObj = false;
    }

    // Material is reflective if it is either mirror or glass in this assignment
    bool reflective = (hitMaterial.materialReflectType == MIRRORREFLECT ||
                       hitMaterial.materialReflectType == GLASSREFLECT);

    // Compute the color at the intersection point based on its material
    // and the lighting in the scene
    vec3 outputColor = calculateColor(hitMaterial, posIntersection,
      normalVector, eyeVector, reflective);

    // A material has a reflection type (as seen above) and a reflectivity
    // attribute. A reflectivity "equal to zero" indicates that the material
    // is neither reflective nor refractive.

    // If a material is neither reflective nor refractive...
    // (1) Scale the output color by the current weight and add it into
    //     the accumulated color.
    // (2) Then break the for loop (i.e. do not trace the ray any further).
    // ----------- STUDENT CODE BEGIN ------------
    if (hitMaterial.reflectivity == 0.0) {
        resColor += resWeight * outputColor;
        break;
    }
    // ----------- Our reference solution uses 4 lines of code.
    // ----------- STUDENT CODE END ------------
    // If the material is reflective or refractive...
    // (1) Use calcReflectionVector to compute the direction of the next
    //     bounce of this ray.
    // (2) Update the ray object with the next starting position and
    //     direction to prepare for the next bounce. You should modify the
    //     ray's origin and direction attributes. Be sure to normalize the
    //     direction vector.
    // (3) Scale the output color by the current weight and add it into
    //     the accumulated color.
    // (4) Update the current weight using the material's reflectivity
    //     so that it is the appropriate weight for the next ray's color.
    // ----------- STUDENT CODE BEGIN ------------
    ray.origin = posIntersection;
    ray.direction = calcReflectionVector(hitMaterial, ray.direction, normalVector, isInsideObj);
    resColor += resWeight * outputColor;
    resWeight *= hitMaterial.reflectivity;
    // ----------- Our reference solution uses 8 lines of code.
    // ----------- STUDENT CODE END ------------
  }

  return resColor;
}

void main() {
  float cameraFOV = 0.8;
  vec3 direction = vec3(v_position.x * cameraFOV * width / height,
                        v_position.y * cameraFOV, 1.0);

  Ray ray;
  ray.origin = vec3(uMVMatrix * vec4(camera, 1.0));
  ray.direction = normalize(vec3(uMVMatrix * vec4(direction, 0.0)));

  // trace the ray for this pixel
  vec3 res = traceRay(ray);

  // paint the resulting color into this pixel
  gl_FragColor = vec4(res.x, res.y, res.z, 1.0);
}
