import * as THREE from 'three';

// A solid drawn as triangles shows its facets wherever it meets another --
// the crease where two faceted tubes cross zig-zags however finely they are
// cut, and close up those creases are most of what shows. An impostor draws
// the solid exactly instead: a box around each instance is rasterised, and
// every fragment casts its view ray at the true surface, discarding misses
// and shading the hit with the hit's own normal and depth. Creases then
// fall wherever the depth test puts them, as smooth as the solids are.
//
// A bar is a capsule: an open cylinder between two spheres of its radius,
// so a chain of cylinders with a sphere at every joint is one smooth tube
// through its bends, and several chains meeting at a sphere are welded
// there by construction.

const CUBE = new THREE.BoxGeometry(2, 2, 2);

export function shell() {
  return CUBE.clone();
}

// A sphere instance maps the cube onto the cube around a ball of the
// instance's scale.
export function sphereMatrix(centre, radius) {
  return new THREE.Matrix4().compose(
    centre,
    new THREE.Quaternion(),
    new THREE.Vector3(radius, radius, radius)
  );
}

// A cylinder instance maps the cube onto the slab around the segment:
// its y axis runs from one end to the other and the other two span the
// radius.
export function cylinderMatrix(from, to, radius) {
  const along = new THREE.Vector3().subVectors(to, from);

  return new THREE.Matrix4().compose(
    new THREE.Vector3().addVectors(from, to).multiplyScalar(0.5),
    new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 1, 0), along.clone().normalize()),
    new THREE.Vector3(radius, 0.5 * along.length(), radius)
  );
}

export function sphereMaterial(parameters) {
  return impostorMaterial(parameters, `
    float surfaceAlong = dot(surfaceRay, impostorCentre);
    float surfaceReach = surfaceAlong * surfaceAlong - dot(impostorCentre, impostorCentre) + impostorRadius * impostorRadius;
    if (surfaceReach < 0.0) discard;
    vec3 surfaceHit = (surfaceAlong - sqrt(surfaceReach)) * surfaceRay;
    vec3 surfaceNormal = normalize(surfaceHit - impostorCentre);
  `);
}

// The ray meets the infinite cylinder where its distance from the axis is
// the radius; a hit beyond either end belongs to the sphere there.
export function cylinderMaterial(parameters) {
  return impostorMaterial(parameters, `
    vec3 surfaceAxis = impostorHigh - impostorLow;
    float surfaceLength = length(surfaceAxis);
    vec3 surfaceUnit = surfaceAxis / surfaceLength;
    vec3 surfaceRayAcross = surfaceRay - dot(surfaceRay, surfaceUnit) * surfaceUnit;
    vec3 surfaceOriginAcross = -impostorLow + dot(impostorLow, surfaceUnit) * surfaceUnit;
    float surfaceQuadratic = dot(surfaceRayAcross, surfaceRayAcross);
    float surfaceLinear = dot(surfaceOriginAcross, surfaceRayAcross);
    float surfaceConstant = dot(surfaceOriginAcross, surfaceOriginAcross) - impostorRadius * impostorRadius;
    float surfaceReach = surfaceLinear * surfaceLinear - surfaceQuadratic * surfaceConstant;
    if (surfaceReach < 0.0) discard;
    vec3 surfaceHit = ((-surfaceLinear - sqrt(surfaceReach)) / surfaceQuadratic) * surfaceRay;
    float surfaceStation = dot(surfaceHit - impostorLow, surfaceUnit);
    if (surfaceStation < 0.0 || surfaceStation > surfaceLength) discard;
    vec3 surfaceNormal = normalize(surfaceHit - impostorLow - surfaceStation * surfaceUnit);
  `);
}

// The instance's own frame is carried to view space in the vertex stage:
// its centre, the two ends of its y axis, and its radius. Everything
// after the hit is the stock material's, lit from the hit rather than from
// the box.
function impostorMaterial(parameters, intersection) {
  const material = new THREE.MeshPhysicalMaterial(parameters);
  material.customProgramCacheKey = () => intersection;

  material.onBeforeCompile = shader => {
    shader.vertexShader = shader.vertexShader
      .replace('void main() {', `
        varying vec3 impostorCentre;
        varying vec3 impostorLow;
        varying vec3 impostorHigh;
        varying float impostorRadius;
        void main() {
      `)
      .replace('#include <project_vertex>', `
        #include <project_vertex>
        mat4 impostorToView = modelViewMatrix * instanceMatrix;
        impostorCentre = (impostorToView * vec4(0.0, 0.0, 0.0, 1.0)).xyz;
        impostorLow = (impostorToView * vec4(0.0, -1.0, 0.0, 1.0)).xyz;
        impostorHigh = (impostorToView * vec4(0.0, 1.0, 0.0, 1.0)).xyz;
        impostorRadius = length((impostorToView * vec4(1.0, 0.0, 0.0, 0.0)).xyz);
      `);

    shader.fragmentShader = shader.fragmentShader
      .replace('void main() {', `
        varying vec3 impostorCentre;
        varying vec3 impostorLow;
        varying vec3 impostorHigh;
        varying float impostorRadius;
        uniform mat4 projectionMatrix;
        void main() {
          vec3 surfaceRay = normalize(-vViewPosition);
          ${intersection}
          vec4 surfaceClip = projectionMatrix * vec4(surfaceHit, 1.0);
          gl_FragDepth = 0.5 * surfaceClip.z / surfaceClip.w + 0.5;
      `)
      .replace('#include <normal_fragment_begin>', `
        #include <normal_fragment_begin>
        normal = surfaceNormal;
        nonPerturbedNormal = surfaceNormal;
      `)
      .replace('#include <lights_fragment_begin>', THREE.ShaderChunk.lights_fragment_begin
        .replace('vec3 geometryPosition = - vViewPosition;', 'vec3 geometryPosition = surfaceHit;')
        .replace('normalize( vViewPosition )', 'normalize( - surfaceHit )')
      );
  };

  return material;
}
