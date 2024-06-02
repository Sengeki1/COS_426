// In this file you will implement traversal and analysis for your assignment.
// Make sure to familiarize yourself with the utility functions in meshUtils.js
// they might be useful for the second part of your assignment!

////////////////////////////////////////////////////////////////////////////////
// Traversal
////////////////////////////////////////////////////////////////////////////////

// Return all vertices on face f
Mesh.prototype.verticesOnFace = function(f) {
  const vertices = [];
  let he = f.halfedge;
  const first = he;
  while (true) {
    vertices.push(he.vertex);
    he = he.next;
    if (he === first) {
      break;
    }
  }
  return vertices;
};

// Return all halfedges on face f
Mesh.prototype.edgesOnFace = function(f) {
  const halfedges = [];

  // ----------- STUDENT CODE BEGIN ------------
  let startEdge = f.halfedge;
  let currentEdge = startEdge;
  do {
    halfedges.push(currentEdge);
    currentEdge = currentEdge.next;
  } while (currentEdge.next !== startEdge);
  // ----------- STUDENT CODE END ------------

  return halfedges;
};

// Return all faces adjacent to input face, not
// including input face.
Mesh.prototype.facesOnFace = function(f) {
  const faces = [];

  // ----------- STUDENT CODE BEGIN ------------
  let startHalfEdge = f.halfedge;
  let currentHalfEdge = startHalfEdge;

  // Traverse around the vertex v to collect all neighboring vertices
  do {
    // Get the vertex pointed to by the current half-edge
    faces.push(currentHalfEdge.face);

    // Move to the next half-edge around the vertex
    currentHalfEdge = currentHalfEdge.next;
  } while (currentHalfEdge !== startHalfEdge);
  // ----------- STUDENT CODE END ------------

  return faces;
};

// Return one-ring neighbors of input vertex, not
// including the input vertex itself
Mesh.prototype.verticesOnVertex = function(v) {
  const vertices = [];

  // ----------- STUDENT CODE BEGIN ------------
  // Start from any half-edge that originates from vertex v
  let startHalfEdge = v.halfedge;
  let currentHalfEdge = startHalfEdge;

  // Traverse around the vertex v to collect all neighboring vertices
  do {
    // Get the vertex pointed to by the current half-edge
    vertices.push(currentHalfEdge.vertex);

    // Move to the next half-edge around the vertex
    currentHalfEdge = currentHalfEdge.next;
  } while (currentHalfEdge !== startHalfEdge);
  // ----------- STUDENT CODE END ------------
  return vertices;
};

// Return all halfedges that point away from v
Mesh.prototype.edgesOnVertex = function(v) {
  const halfedges = [];

  // ----------- STUDENT CODE BEGIN ------------
  let startHalfEdge = v.halfedge;
  let currentHalfEdge = startHalfEdge;

  // Traverse around the vertex v to collect all outgoing halfedges
  do {
    halfedges.push(currentHalfEdge.halfedge);
    currentHalfEdge = currentHalfEdge.next;
  } while (currentHalfEdge !== startHalfEdge);
  // ----------- STUDENT CODE END ------------

  return halfedges;
};

// Return all faces that include v as a vertex.
Mesh.prototype.facesOnVertex = function(v) {
  const faces = [];

  // ----------- STUDENT CODE BEGIN ------------
  for (const face of this.faces) {
    let vetFace = this.verticesOnFace(face)
    if (vetFace.includes(v)) {
      faces.push(face)
    }
  }
  // ----------- STUDENT CODE END ------------

  return faces;
};

// Return the vertices that form the endpoints of a given edge
Mesh.prototype.verticesOnEdge = function(e) {
  const vertices = [];

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 2 lines of code.
  // ----------- STUDENT CODE END ------------

  return vertices;
};

// Return the faces that include a given edge
Mesh.prototype.facesOnEdge = function(e) {
  const faces = [];
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 2 lines of code.
  // ----------- STUDENT CODE END ------------
  return faces;
};

// Return the edge pointing from v1 to v2
Mesh.prototype.edgeBetweenVertices = function(v1, v2) {
  let out_he = undefined;
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 11 lines of code.
  // ----------- STUDENT CODE END ------------
  return out_he;
};

////////////////////////////////////////////////////////////////////////////////
// Analysis
////////////////////////////////////////////////////////////////////////////////

// Return the surface area of a provided face f.
Mesh.prototype.calculateFaceArea = function(f) {
  let area = 0.0;
  // ----------- STUDENT CODE BEGIN ------------
  const v0 = f.halfedge.vertex.position;
  const v1 = f.halfedge.next.vertex.position;
  const v2 = f.halfedge.next.next.vertex.position;
  
  // Use the cross product of two edges to calculate the area
  const edge1 = new THREE.Vector3().subVectors(v1, v0);
  const edge2 = new THREE.Vector3().subVectors(v2, v0);
  const crossProduct = new THREE.Vector3().crossVectors(edge1, edge2);

  // Area of the triangle is half the magnitude of the cross product
  area = crossProduct.length() * 0.5;
  // ----------- STUDENT CODE END ------------
  f.area = area;
  return area;
};

// Update the area attributes of all faces in the mesh
Mesh.prototype.calculateFacesArea = function() {
  for (let i = 0; i < this.faces.length; ++i) {
    this.calculateFaceArea(this.faces[i]);
  }
};

// Calculate the vertex normal at a given vertex,
// using the face normals of bordering faces, weighted by face area
Mesh.prototype.calculateVertexNormal = function(v) {
  var v_normal = new THREE.Vector3(0, 0, 0);
  // ----------- STUDENT CODE BEGIN ------------
  function normalizeVector(vector) {
    var length = Math.sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z)
    if (length !== 0) {
        vector.x /= length
        vector.y /= length
        vector.z /= length
    }
    return vector
  }

  var adjacent = this.facesOnVertex(v) // Get all the faces adjacent to the vertex
  for (let face of adjacent) {
    var faceNormal = this.calculateFaceNormal(face) // calculate the normal of the given face
    var weightedNormal = faceNormal.multiplyScalar(face.area)
    v_normal.add(weightedNormal)
  }
  normalizeVector(v_normal)
  // ----------- STUDENT CODE END ------------
  v.normal = v_normal;
  return v_normal;
};

// update the vertex normals of every vertex in the mesh
Mesh.prototype.updateVertexNormals = function() {
  for (let i = 0; i < this.vertices.length; ++i) {
    this.calculateVertexNormal(this.vertices[i]);
  }
};

// compute the average length of edges touching v
Mesh.prototype.averageEdgeLength = function(v) {
  let avg = 0.0;

  // ----------- STUDENT CODE BEGIN ------------
  let totalLength = 0.0
  let count = 0

  let startEdge = v.halfedge
  let currentEdge = startEdge

  while (currentEdge) {
    const vertexPosition = currentEdge.vertex.position
    const nextVertexPosition = currentEdge.next ? currentEdge.next.vertex.position : null

    if (nextVertexPosition) {
      // Calculate Euclidean distance manually
      const dx = vertexPosition.x - nextVertexPosition.x
      const dy = vertexPosition.y - nextVertexPosition.y
      const dz = vertexPosition.z - nextVertexPosition.z
      const edgeLength = Math.sqrt(dx * dx + dy * dy + dz * dz)
      totalLength += edgeLength
      count++
    }

    currentEdge = currentEdge.next
    
    if (currentEdge === startEdge) {
      break; // Break the loop if we've reached the starting edge again
    }
  }

  avg = count > 0 ? totalLength / count : 0.0
  // ----------- STUDENT CODE END ------------
  return avg;
};

////////////////////////////////////////////////////////////////////////////////
// Topology
////////////////////////////////////////////////////////////////////////////////

// Given a face in the shape of an arbitrary polygon,
// split that face so it consists only of several triangular faces. 
Mesh.prototype.triangulateFace = function(f) {
  // ----------- STUDENT CODE BEGIN ------------
  const halfedges = this.edgesOnFace(f);

  // Check if the face is already triangular
  if (halfedges.length === 3) {
    return; 
  }
  const startVertex = halfedges[0].vertex;

  // Iterate over the remaining vertices of the face to create triangles
  for (let i = 1; i < halfedges.length - 1; i++) {
    const v2 = halfedges[i].vertex;
    const v3 = halfedges[i + 1].vertex;

    // Create a new triangle face using the start vertex and the current pair of vertices
    const newFace = this.addFace();
    this.addHalfEdge(startVertex, v2, newFace);
    this.addHalfEdge(v2, v3, newFace);
    this.addHalfEdge(v3, startVertex, newFace);
  }

  // Remove the original face from the mesh
  this.removeFace(f);
  // ----------- STUDENT CODE END ------------
};
