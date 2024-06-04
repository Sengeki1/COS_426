var Filters = Filters || {};

// Space for your helper functions
// ----------- STUDENT CODE BEGIN ------------
function createRotationMatrixX(angle) {
  const cos = Math.cos(angle);
  const sin = Math.sin(angle);
  return [
    [1, 0, 0],
    [0, cos, -sin],
    [0, sin, cos]
  ];
}

function createRotationMatrixY(angle) {
  const cos = Math.cos(angle);
  const sin = Math.sin(angle);
  return [
    [cos, 0, sin],
    [0, 1, 0],
    [-sin, 0, cos]
  ];
}

function createRotationMatrixZ(angle) {
  const cos = Math.cos(angle);
  const sin = Math.sin(angle);
  return [
    [cos, -sin, 0],
    [sin, cos, 0],
    [0, 0, 1]
  ];
}

function cotangentWeightCalculation(v1, v2, mesh) {
  const neighbors1 = mesh.verticesOnVertex(v1)
  const neighbors2 = mesh.verticesOnVertex(v2)

  const commonNeighbors = []
  for (let neighbor of neighbors1) {
    if (neighbors2.includes(neighbor)) commonNeighbors.push(neighbor)
  }

  if (commonNeighbors.length !== 2) {
    return 0;
  }

  const [v3, v4] = commonNeighbors;
  const cotAlpha = cotangent(v1.position, v2.position, v3.position);
  const cotBeta = cotangent(v1.position, v2.position, v4.position);

  return cotAlpha + cotBeta
}

function cotangent(p1, p2, p3) {
  const u = new THREE.Vector3().subVectors(p1, p3);
  const v = new THREE.Vector3().subVectors(p2, p3);
  const dot = u.dot(v)
  const cross = new THREE.Vector3().crossVectors(u, v)
  const crossLength = cross.length()
  if (crossLength === 0) {
    return 0;
  }
  const cotangent = dot / Math.abs(crossLength)
  return cotangent;
}

function distanceTo(vertex1, vertex2) {
  const dx = vertex2.x - vertex1.x;
  const dy = vertex2.y - vertex1.y;
  const dz = vertex2.z - vertex1.z;

  // Calculate the Euclidean distance
  const distance = Math.sqrt(dx * dx + dy * dy + dz * dz);

  return distance;
}

// Function to compute principal curvatures for a vertex using Euler operations
function computePrincipalCurvatures(vertex, mesh) {
  const edgeLengths = [];
  const dihedralAngles = [];

  // Iterate over the halfedges incident to the vertex
  let startHalfedge = vertex.halfedge;
  let halfedge = startHalfedge;
  const visitedHalfedges = new Set();
  
  while (halfedge && !visitedHalfedges.has(halfedge)) {
    visitedHalfedges.add(halfedge);

    // Access the start and end vertices of the edge
    const v1 = halfedge.vertex.position;
    const v2 = halfedge.opposite.vertex.position;

    // Compute the length of the edge
    const edgeLength = distanceTo(v1, v2);
    edgeLengths.push(edgeLength);

    // Compute the dihedral angle for adjacent faces
    const dihedralAngle = computeDihedralAngle(halfedge, mesh);
    dihedralAngles.push(dihedralAngle);

    // Move to the next halfedge incident to the vertex
    halfedge = halfedge.opposite ? halfedge.opposite.next : null;

    // If we return to the starting halfedge, we should stop
    if (halfedge === startHalfedge) {
      break;
    }
  }

  // Estimate principal curvatures using edge lengths and dihedral angles
  const principalCurvatures = estimateCurvature(edgeLengths, dihedralAngles);

  // Return the computed principal curvatures
  return principalCurvatures;
}

// Function to compute the dihedral angle between adjacent faces incident to an edge
function computeDihedralAngle(edge, mesh) {
  // Find the two adjacent faces sharing the edge
  const face1 = edge.face;
  const face2 = edge.opposite.face;

  // Assuming the mesh is triangular, get the normal vectors of the faces
  const normal1 = mesh.calculateFaceNormal(face1);
  const normal2 = mesh.calculateFaceNormal(face2);

  // Compute the angle between the normal vectors using dot product
  const dotProduct = normal1.dot(normal2);
  const angle = Math.acos(Math.min(Math.max(dotProduct, -1), 1)); // Ensure the value is within [-1, 1] for safe acos

  return angle;
}

// Function to estimate principal curvatures using edge lengths and dihedral angles
function estimateCurvature(edgeLengths, dihedralAngles) {
  const meanEdgeLength = edgeLengths.reduce((acc, val) => acc + val, 0) / edgeLengths.length;
  const meanDihedralAngle = dihedralAngles.reduce((acc, val) => acc + val, 0) / dihedralAngles.length;

  // Estimate curvature based on linear and quadratic relationships
  const linearTerm = 1 / meanEdgeLength; // Linear term based on inverse of mean edge length
  const quadraticTerm = Math.pow(meanDihedralAngle, 2); // Quadratic term based on squared mean dihedral angle

  // Combine linear and quadratic terms to estimate curvature
  const curvature = linearTerm + quadraticTerm;

  return curvature;
}

function mapCurvatureToColor(curvature) {
  const lowCurvatureColor = [0, 255, 255]; // Blue
  const highCurvatureColor = [255, 0, 255]; // Red

  const minCurvature = -1; // Minimum possible curvature value
  const maxCurvature = 1; // Maximum possible curvature value

  // Map the curvature value to a color within the defined range
  const mappedColor = [];

  for (let i = 0; i < 3; i++) {
    const lowColorComponent = lowCurvatureColor[i];
    const highColorComponent = highCurvatureColor[i];
    
    // Linearly interpolate between low and high color components based on curvature value
    const colorComponent = lowColorComponent + (curvature - minCurvature) * (highColorComponent - lowColorComponent) / (maxCurvature - minCurvature);
    
    // Ensure the color component is within the valid range [0, 255]
    mappedColor[i] = Math.min(255, Math.max(0, Math.round(colorComponent)));
  }

  return mappedColor;
}

function getCommonFace(v1, v2, mesh) {
  let face1 = mesh.facesOnVertex(v1);
  let face2 = mesh.facesOnVertex(v2);

  for (let face of face2) {
    if (face1.includes(face)) return face;
  }
}
// ----------- STUDENT CODE END ------------

// Translate all selected vertices in the mesh by the given x,y,z offsets.
Filters.translation = function(mesh, x, y, z) {
  const t = new THREE.Vector3(x, y, z);

  const verts = mesh.getModifiableVertices();

  const n_vertices = verts.length;
  for (let i = 0; i < n_vertices; ++i) {
    verts[i].position.add(t);
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Given x,y,z, the desired rotation around each axis, in radians,
// apply this rotation to all selected vertices in the mesh.
Filters.rotation = function(mesh, x, y, z) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  function applyMatrixToVertex(vertex, matrix) {
    let x = vertex.x
    let y = vertex.y
    let z = vertex.z
    return new THREE.Vector3(
      x * matrix[0][0] + y * matrix[0][1] + z * matrix[0][2],
      x * matrix[1][0] + y * matrix[1][1] + z * matrix[1][2],
      x * matrix[2][0] + y * matrix[2][1] + z * matrix[2][2]
    )
  }
  const rotationMatrixX = createRotationMatrixX(x);
  const rotationMatrixY = createRotationMatrixY(y);
  const rotationMatrixZ = createRotationMatrixZ(z);
  for (let i = 0; i < verts.length; i++) {
    let vertex = verts[i].position.clone()

    vertex = applyMatrixToVertex(vertex, rotationMatrixX);
    vertex = applyMatrixToVertex(vertex, rotationMatrixY);
    vertex = applyMatrixToVertex(vertex, rotationMatrixZ);

    verts[i].position.copy(vertex)
  }
  // ----------- STUDENT CODE END ------------
  //Gui.alertOnce("Rotation is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Uniformly scale the position of all selected vertices in the mesh
// by the provided scale factor s
Filters.scale = function(mesh, s) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  const t = new THREE.Vector3(s, s, s);
  for (let i = 0; i < verts.length; i++) {
    verts[i].position.multiply(t)
  }
  // ----------- STUDENT CODE END ------------
  //Gui.alertOnce("Scaling is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// estimate the per-vertex gaussian vurvature of the mesh at each vertex.
// set that vertex's color to some value based on its curvature value.
// (the precise mapping of curvature to color is left to you)
Filters.curvature = function(mesh) {
  // ----------- STUDENT CODE BEGIN ------------
  const vertices = mesh.getModifiableVertices();
  for (let vertex of vertices) {
    const principalCurvatures = computePrincipalCurvatures(vertex, mesh)
    const gaussianCurvature = principalCurvatures[0] * principalCurvatures[1]
    const vertexColor = mapCurvatureToColor(gaussianCurvature)
    vertex.color = vertexColor
  }
  // ----------- STUDENT CODE END ------------
  //Gui.alertOnce("Curvature is not implemented yet");
};

// Apply a random offset to each selected vertex in the direction of its normal
// scale the random offset by the provided factor and by
// the average length of edges at that vertex
Filters.noise = function(mesh, factor) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  for (let i = 0; i < verts.length; i++) {
    const avgEdgeLength = mesh.averageEdgeLength(verts[i])
    let randomValue = Math.random() * 2 - 1 // Generate a random value in the range [-1, 1)
    let offset = verts[i].normal.clone().multiplyScalar(randomValue * avgEdgeLength * factor)
    verts[i].position.add(offset)
  }
  mesh.updateVertexNormals()
  // ----------- STUDENT CODE END ------------
  //Gui.alertOnce("Noise is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Smooth the mesh using the specified weighting scheme.
// In the standard case, this is done using uniform Laplacian smoothing,
// by moving each vertex towards the average position of its neighbors.
//
// Arguments:
//  - mesh: the mesh to smooth
//  - iter: the number of iterations of smoothing to apply
//  - delta: a scaling factor for the amount of smoothing
//  - curvFlow: a bool. if true, use cotangent weights instead of uniform (requires triangular mesh)
//  - scaleDep: a bool. if true, scale offsets differently for each vertex (see spec.)
//  - implicit: a bool. if true, perform implicit smoothing (see spec.)
//
// Note that the reference solution calls a giant utility function so the line
// count is not terribly representative of the true solution
//
// For implicit, you will want to compute I - M*L*delta, where I is the identity
// matrix, M is a diagonal "mass" matrix, and L is a Laplacian matrix. Then
// you will want to call math.lup() on your result in order to decompose the
// matrix. Finally, call math.lusolve() to compute the X,Y, and Z positions of
// vertices. Note that the decomposition step allows for fast solving multiple
// times. It would be possible to replace a few of these steps with simple matrix
// inversion; however, matrix inversion is far slower and less numerically stable
//
Filters.smooth = function(mesh, iter, delta, curvFlow, scaleDep, implicit) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  const n = verts.length
  for (let i = 0; i < iter; i++) {
    if (implicit) {
      // Implicit smoothing setup
      const M = new Array(n).fill(0).map(() => 1);  // Assuming uniform mass
      const L = new Array(n).fill().map(() => new Array(n).fill(0));
      const I = new Array(n).fill().map((_, idx) => {
        const row = new Array(n).fill(0);
        row[idx] = 1;
        return row;
      });

      for (let j = 0; j < n; j++) {
        let neighbors = mesh.verticesOnVertex(verts[j]);
        let weightSum = 0;
        
        for (let neighbor of neighbors) {
          const k = verts.indexOf(neighbor);
          L[j][k] = -1;
          weightSum += 1;
        }
        L[j][j] = weightSum;

        const A = new Array(n).fill().map(() => new Array(n).fill(0));
        for (let k = 0; k < n; k++) {
          A[j][k] = I[j][k] - delta * M[j] * L[j][k];
        }
      }
      
      // Step 3: Decompose A using LU decomposition
      const { L: Lmatrix, U: Umatrix, P: Pmatrix } = math.lup(A);

      // Step 4: Solve for new positions using LU decomposition
      for (let axis of ['x', 'y', 'z']) {
        const b = verts.map(v => v.position[axis]);
        const Pb = math.multiply(Pmatrix, b);
        const y = math.lsolve(Lmatrix, Pb);
        const newPos = math.usolve(Umatrix, y);

        for (let j = 0; j < n; j++) {
          verts[j].position[axis] = newPos[j];
        }
      }
    } else {
      // Explicit smoothing
      for (let j = 0; j < verts.length; j++) {
        let neighbors = mesh.verticesOnVertex(verts[j]);
        let averagePosition = new THREE.Vector3(0, 0, 0);
        let weight = 0;
        let areaSum = 0;

        if (scaleDep) {
          let adjacent = mesh.facesOnVertex(verts[j]);
          for (let face of adjacent) {
            let area = mesh.calculateFaceArea(face);
            areaSum += area;
          }
        }

        if (!curvFlow) {
          for (let neighbor of neighbors) {
            averagePosition.add(neighbor.position);
          }
          weight = neighbors.length;
        } else {
          for (let neighbor of neighbors) {
            let w = cotangentWeightCalculation(verts[j], neighbor, mesh);
            averagePosition.addScaledVector(neighbor.position, w);
            weight += w;
          }
        }

        averagePosition.divideScalar(weight);
        let direction = new THREE.Vector3(
          averagePosition.x - verts[j].position.x,
          averagePosition.y - verts[j].position.y,
          averagePosition.z - verts[j].position.z
        );

        if (scaleDep) {
          let adjacent = mesh.facesOnVertex(verts[j]);
          const A = areaSum / adjacent.length;
          direction.multiplyScalar(delta * (areaSum / A));
        } else {
          direction.multiplyScalar(delta);
        }

        verts[j].position.copy(verts[j].position.clone().add(direction));
      }
    }
  }
  // ----------- STUDENT CODE END ------------
  //Gui.alertOnce("Smooth is not implemented yet");
  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Sharpen the mesh by moving selected vertices away from the average position
// of their neighbors (i.e. Laplacian smoothing in the negative direction)
Filters.sharpen = function(mesh, iter, delta) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  for (let i = 0; i < iter; i++) {
    for (let j = 0; j < verts.length; j++) {
      let neighbors = mesh.verticesOnVertex(verts[j])

      let averagePosition = new THREE.Vector3(0, 0, 0)
      for (let neighbor of neighbors) {
        averagePosition.add(neighbor.position)
      }
      averagePosition.divideScalar(neighbors.length)
      let direction = new THREE.Vector3(
        verts[j].position.x - averagePosition.x,
        verts[j].position.y - averagePosition.y,
        verts[j].position.z - averagePosition.z
      ).multiplyScalar(delta)

      verts[j].position.copy(verts[j].position.clone().add(direction))
    }
  }
  // ----------- STUDENT CODE END ------------
  //Gui.alertOnce("Sharpen is not implemented yet");
  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Move every selected vertex along its normal direction
// Scale the amount by the provided factor and average edge length at that vertex
Filters.inflate = function(mesh, factor) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  for (let i = 0; i < verts.length; i++) {     
    verts[i].position.add(new THREE.Vector3(verts[i].normal.x * factor,verts[i].normal.y * factor,verts[i].normal.z * factor))
  }
  mesh.updateVertexNormals()
  // ----------- STUDENT CODE END ------------
  //Gui.alertOnce("Inflate is not implemented yet");
  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// rotate selected vertices around the Y axis by an amount
// proportional to its Y value times the scale factor.
Filters.twist = function(mesh, factor) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  for (let i = 0; i < verts.length; i++) {
    let vertex = verts[i].position
    const matrix = createRotationMatrixY(vertex.y * factor)

    verts[i].position = new THREE.Vector3(
        vertex.x * matrix[0][0] + vertex.y * matrix[0][1] + vertex.z * matrix[0][2],
        vertex.x * matrix[1][0] + vertex.y * matrix[1][1] + vertex.z * matrix[1][2],
        vertex.x * matrix[2][0] + vertex.y * matrix[2][1] + vertex.z * matrix[2][2]
      )
  }
  // ----------- STUDENT CODE END ------------
  //Gui.alertOnce("Twist is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// warp a mesh using a nonlinear mapping of your choice
Filters.wacky = function(mesh, factor) {
  // ----------- STUDENT CODE BEGIN ------------
  const verts = mesh.getModifiableVertices()
  for (let i = 0; i < verts.length; i++) {
    verts[i].position.multiply(new THREE.Vector3(Math.cos(factor), Math.sin(factor), 1))
  }
  // ----------- STUDENT CODE END ------------
  //Gui.alertOnce("Wacky is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Convert the selected faces from arbitrary polygons into all triangles
Filters.triangulate = function(mesh) {
  const faces = mesh.getModifiableFaces();

  // ----------- STUDENT CODE BEGIN ------------
  faces.forEach(face => {
    if (face.halfedge && mesh.edgesOnFace(face).length > 3) {
      console.log(face)
      mesh.triangulateFace(face);
    }
  });

  // Compute principal curvatures for each vertex
  const vertices = mesh.getModifiableVertices();
  vertices.forEach(vertex => {
    // Compute principal curvatures for the vertex
    const principalCurvatures = computePrincipalCurvatures(vertex, mesh);

    // Store the computed principal curvatures in the vertex
    vertex.curvature = principalCurvatures;
  });
  // ----------- STUDENT CODE END ------------
  //Gui.alertOnce("triangulate is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// wrapper for splitEdgeMakeVert in mesh.js
Filters.splitEdge = function(mesh) {
  const verts = mesh.getSelectedVertices();

  if (verts.length === 2) {
    mesh.splitEdgeMakeVert(verts[0], verts[1], 0.5);
  } else {
    console.log("ERROR: to use split edge, select exactly 2 adjacent vertices");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// wrapper for joinEdgeKillVert in mesh.js
Filters.joinEdges = function(mesh) {
  const verts = mesh.getSelectedVertices();

  if (verts.length === 3) {
    let v0 = verts[0],
      v1 = verts[1],
      v2 = verts[2];

    const he01 = mesh.edgeBetweenVertices(v0, v1);
    const he12 = mesh.edgeBetweenVertices(v1, v2);

    if (he01) {
      if (he12) {
        mesh.joinEdgeKillVert(verts[0], verts[1], verts[2]);
      } else {
        mesh.joinEdgeKillVert(verts[1], verts[0], verts[2]);
      }
    } else {
      if (he12) {
        mesh.joinEdgeKillVert(verts[0], verts[2], verts[1]);
      } else {
        console.log(
          "ERROR: to use join edge, select exactly 3 vertices such that one only has edges to the other two"
        );
      }
    }
  } else {
    console.log("ERROR: to use join edge, select exactly 3 vertices");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// wrapper for splitFaceMakeEdge in mesh.js
Filters.splitFace = function(mesh) {
  const verts = mesh.getSelectedVertices();
  const faces = mesh.getModifiableFaces();

  if (verts.length === 2 && faces.length === 1) {
    mesh.splitFaceMakeEdge(faces[0], verts[0], verts[1]);
  } else {
    console.log("ERROR: to use split face, select exactly 1 face and 2 nonadjacent vertices on it");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// wrapper for joinFaceKillEdge in mesh.js
Filters.joinFaces = function(mesh) {
  const verts = mesh.getSelectedVertices();
  const faces = mesh.getModifiableFaces();

  if (verts.length === 2 && faces.length === 2) {
    mesh.joinFaceKillEdge(faces[0], faces[1], verts[0], verts[1]);
  } else {
    console.log(
      "ERROR: to use split face, select exactly 2 adjacent faces the 2 vertices between them"
    );
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// extrude the selected faces from the mesh in the direction of their normal
// vector, scaled by the provided factor.
// See the spec for more detail.
Filters.extrude = function(mesh, factor) {
  const _faces = mesh.getModifiableFaces();

  // ----------- STUDENT CODE BEGIN ------------
  let faces = [..._faces];
  for(let face of faces)
    {
      let halfedges = mesh.edgesOnFace(face);
      console.log(halfedges)
      let new_verts = [];
      
      for(let i=0; i < halfedges.length; i++)
      {
        let new_vert = mesh.splitEdgeMakeVert(halfedges[i].vertex,
        halfedges[i].opposite.vertex, 0);
        
        let adj_face = new_vert.halfedge.opposite.face;
        
        mesh.splitFaceMakeEdge(adj_face, new_vert.halfedge.vertex, 
        new_vert.halfedge.opposite.next.vertex);
        new_verts.push(new_vert);
      }
      
      for(let i = 0; i < new_verts.length; i++)
      {
        mesh.splitFaceMakeEdge(face, new_verts[i], new_verts[(i + 1) % new_verts.length]);
        
        mesh.joinFaceKillEdge(
          new_verts[i].halfedge.opposite.next.next.opposite.face,
          new_verts[i].halfedge.opposite.face,
          new_verts[i].halfedge.opposite.next.vertex,
          new_verts[i].halfedge.vertex
          );
      }
        
      console.log(faces.length);
      for(let new_vert of new_verts)
      {
          new_vert.position.addScaledVector(face.normal.normalize(), factor);
      }
    }
  // ----------- STUDENT CODE END ------------
  //Gui.alertOnce("Extrude is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Truncate the selected vertices of the mesh by "snipping off" corners
// and replacing them with faces. factor specifies the size of the truncation.
// See the spec for more detail.
Filters.truncate = function(mesh, factor) {
  const vertis = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  const meshCopy = new Mesh();
  meshCopy.copy(mesh);

  // Get the vertices list from the copied mesh
  const verts = meshCopy.getModifiableVertices();
  for (let i = 0; i < vertis.length; i++) {
    const vertex = verts[i];
    const neighbors = meshCopy.verticesOnVertex(vertex);
    if (neighbors.length < 3) {
      continue;
    }

    // Calculate new vertices along each edge
    const newVertices = [];
    for (let j = 0; j < neighbors.length; j++) {
      const neighbor = neighbors[j];
      const newVertex = meshCopy.splitEdgeMakeVert(vertex, neighbor, factor);
      if (newVertex === false) {
        return;
      }
      newVertices.push(newVertex);
    }
    // Create new faces to truncate the corner
    const numNeighbors = neighbors.length;
    for (let j = 0; j < numNeighbors - 1; j++) {
      const v1 = newVertices[j];
      const v2 = newVertices[j + 1];

      // Create a new face with vertices v1, v2, and the original vertex
      const face = getCommonFace(v1, v2, meshCopy)
      const newFace = meshCopy.splitFaceMakeEdge(face, v1, v2, undefined, false);
      if (!newFace) {
        return;
      }
    }
  }
  mesh.clear();
  mesh.copy(meshCopy);

  for (let face of mesh.getModifiableFaces()) {
    mesh.calculateFaceNormal(face)
  }
  // ----------- STUDENT CODE END ------------
  //Gui.alertOnce("Truncate is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Apply the bevel operation to the mesh, scaling the degree of bevelling by factor
Filters.bevel = function ( mesh, factor ) {

    var verts = mesh.getModifiableVertices();

    // ----------- STUDENT CODE BEGIN ------------
    // Store original vertices and their truncated counterparts
    const truncatedVertices = [];

    // Truncate each vertex along the edge
    for (const vertex of verts) {
      const newPos = new THREE.Vector3().lerpVectors(vertex, vertex+1, factor);
      console.log(newPos)
      const truncatedVertex = mesh.splitEdgeMakeVert(vertex, newPos);
      truncatedVertices.push(truncatedVertex);
    }
    // Create new faces to bevel the edge
    const face = verts[0].halfedge.face;
    console.log(face)
    const newFace = mesh.splitFaceMakeEdge(face, truncatedVertices[0], truncatedVertices[1]);
    if (!newFace) {
      console.error(`Failed to split face with vertices ${truncatedVertices[0]}, ${truncatedVertices[1]}`);
      return;
    }
    // ----------- STUDENT CODE END ------------
    //Gui.alertOnce ('Bevel is not implemented yet');

    mesh.calculateFacesArea();
    mesh.updateNormals();
};

// Split the longest edges in the mesh into shorter edges.
// factor is a float in [0,1]. it tells the proportion
// of the total number of edges in the mesh that should be split.
Filters.splitLong = function(mesh, factor) {
  // ----------- STUDENT CODE BEGIN ------------
  const faces = mesh.getModifiableFaces()
  let edges = []
  for (let face of faces) {
    let edge = mesh.edgesOnFace(face)
    edges.push(edge) 
  }
  const totalEdges = edges.length
  let splits = 0

  while (splits < factor * totalEdges) {
    // Find the longest edge
    let longestEdge = null;
    let maxLengthSquared = 0;

    for (const edge of edges) {
        const lengthSquared = edge.lengthSquared();
        if (lengthSquared > maxLengthSquared) {
            maxLengthSquared = lengthSquared;
            longestEdge = edge;
        }
    }

    // Split the longest edge
    const midpoint = longestEdge.midpoint();
    const v = mesh.splitEdgeMakeVert(longestEdge.vertex1, longestEdge.vertex2, 0.5);

    // Increment split count
    splits++;

    // Update edges
    edges.push(longestEdge.halfedge.opposite.edge);
    edges.push(longestEdge.halfedge.edge);
  }
  // ----------- STUDENT CODE END ------------
  //Gui.alertOnce("Split Long Edges is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Triangulate a mesh, and apply triangular subdivision to its faces.
// Repeat for the specified number of levels.
Filters.triSubdiv = function(mesh, levels) {
  Filters.triangulate(mesh);

  for (let l = 0; l < levels; l++) {
    const faces = mesh.getModifiableFaces();
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 43 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce("Triangle subdivide is not implemented yet");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Triangulate the mesh and apply loop subdivision to the faces
// repeat for the specified number of levels.
Filters.loop = function(mesh, levels) {
  Filters.triangulate(mesh);

  for (let l = 0; l < levels; l++) {
    const faces = mesh.getModifiableFaces();
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 123 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce("Triangle subdivide is not implemented yet");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Requires a quad mesh. Apply quad subdivision to the faces of the mesh.
// Repeat for the specified number of levels.
Filters.quadSubdiv = function(mesh, levels) {
  for (let l = 0; l < levels; l++) {
    const faces = mesh.getModifiableFaces();
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 55 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce("Quad subdivide is not implemented yet");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Apply catmull clark subdivision to the faces of the mesh.
// Repeat for the specified number of levels.
Filters.catmullClark = function(mesh, levels) {
  for (let l = 0; l < levels; l++) {
    const faces = mesh.faces;
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 102 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce("Catmull-Clark subdivide is not implemented yet");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// ================= internal functions =======================

// internal function for selecting faces in the form of a loop
Filters.procFace = function(mesh, f) {
  const faceFlags = new Array(mesh.faces.length);
  for (let i = 0; i < mesh.faces.length; i++) {
    faceFlags[i] = 0;
  }
  let sum = f.area;
  const start_he = f.halfedge.opposite.next;
  let curr_he = start_he;
  do {
    if (faceFlags[curr_he.face.id] > 0) {
      break;
    }
    sum += curr_he.face.area;
    curr_he.face.selected = true;
    faceFlags[curr_he.face.id]++;
    const last_he = curr_he;
    curr_he = curr_he.opposite.next;
    if (curr_he.face == f) {
      curr_he = last_he.next.opposite.next;
    }
  } while (true);
};

Filters.parseSelected = function(sel) {
  if (sel === undefined || sel.replace === undefined) {
    return [];
  }
  if (typeof sel === "number") {
    return [sel];
  }
  // sel = sel.replace(/[\(\)]/g,'');
  sel = sel.split(",");
  const parsedSel = [];
  for (let i = 0; i < sel.length; i++) {
    const idx = parseInt(sel[i]);
    if (!isNaN(idx)) {
      parsedSel.push(idx);
    }
  }
  return parsedSel;
};

// internal filter for updating selection
Filters.selection = function(mesh, vertIdxs, faceIdxs) {
  mesh.setSelectedVertices(Filters.parseSelected(vertIdxs));
  mesh.setSelectedFaces(Filters.parseSelected(faceIdxs));
};

// internal filter for setting display settings
Filters.displaySettings = function(
  mesh,
  showLabels,
  showHalfedge,
  shading,
  showVN,
  showFN,
  showGrid,
  showVertDots,
  showAxes,
  showVC,
  meshColor
) {
  Main.displaySettings.showIdLabels = showLabels;
  Main.displaySettings.wireframe = showHalfedge;
  Main.displaySettings.shading = shading;
  Main.displaySettings.showVN = showVN;
  Main.displaySettings.showFN = showFN;
  Main.displaySettings.showGrid = showGrid;
  Main.displaySettings.showVertDots = showVertDots;

  Main.displaySettings.showAxes = showAxes;
  Main.displaySettings.showVC = showVC;
  // Main.displaySettings.meshColor = meshColor;

  // Main.refreshDisplaySettings();
};
