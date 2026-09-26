import * as THREE from 'three';
import { RoomEnvironment } from 'three/addons/environments/RoomEnvironment.js';

export const BACKGROUND = 0x0e0f11;

export const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
renderer.setSize(innerWidth, innerHeight);
renderer.toneMapping = THREE.ACESFilmicToneMapping;
renderer.toneMappingExposure = 1.15;
document.body.appendChild(renderer.domElement);

export const scene = new THREE.Scene();
scene.background = new THREE.Color(BACKGROUND);

export const camera = new THREE.PerspectiveCamera(38, innerWidth / innerHeight, 0.02, 100);
camera.position.set(0, 0.85, 3.9);
camera.lookAt(0, 0, 0);

const environment = new THREE.PMREMGenerator(renderer);
scene.environment = environment.fromScene(new RoomEnvironment(), 0.04).texture;

const key = new THREE.DirectionalLight(0xffffff, 2.2);
key.position.set(3, 4, 2);
scene.add(key);
scene.add(new THREE.AmbientLight(0xffffff, 0.25));

// Everything a drawing puts in the scene hangs here, so clearing one
// drawing away cannot disturb the lights or the environment.
export const group = new THREE.Group();
scene.add(group);

// The figure turns, not the camera: dragging rolls it arcball-style about
// whatever axis the drag implies, a flick hands the idle spin that axis
// and speed, and left alone the spin eases back to a slow drift about the
// last axis it was given. Letting go of a figure that is held still --
// a click, or a drag that comes to rest before the release -- stops it,
// and it stays stopped until the next flick. The wheel walks the camera
// in and out -- all the way inside.
const IDLE = 0.06;
const FASTEST = 0.9;
const EASE = 4;
const HELD = 0.1;
const RADIUS = { nearest: 0.05, farthest: 8 };

const spin = { axis: new THREE.Vector3(0, 1, 0), speed: IDLE, resting: IDLE };
const swing = new THREE.Quaternion();
const dragAxis = new THREE.Vector3();

let dragging = false;
let previous = null;

renderer.domElement.addEventListener('pointerdown', event => {
  dragging = true;
  previous = { x: event.clientX, y: event.clientY, time: performance.now() };
  renderer.domElement.setPointerCapture(event.pointerId);
});

addEventListener('pointerup', () => {
  if (!dragging) return;

  dragging = false;

  const held = (performance.now() - previous.time) / 1000 > HELD;
  spin.resting = held ? 0 : IDLE;

  if (held) spin.speed = 0;
});

renderer.domElement.addEventListener('pointermove', event => {
  if (!dragging) return;

  const time = performance.now();
  const acrossX = event.clientX - previous.x;
  const acrossY = event.clientY - previous.y;
  const seconds = Math.max((time - previous.time) / 1000, 1e-3);
  previous = { x: event.clientX, y: event.clientY, time };

  const span = Math.hypot(acrossX, acrossY);

  if (span < 1e-6) return;

  // A drag turns the figure about the screen axis perpendicular to it.
  dragAxis.set(acrossY, acrossX, 0).normalize().applyQuaternion(camera.quaternion);
  const angle = span / (0.35 * Math.min(innerWidth, innerHeight));

  swing.setFromAxisAngle(dragAxis, angle);
  group.quaternion.premultiply(swing);

  spin.axis.copy(dragAxis);
  spin.speed = Math.min(FASTEST, 0.65 * spin.speed + 0.35 * (angle / seconds));
});

renderer.domElement.addEventListener('wheel', event => {
  event.preventDefault();

  camera.position.setLength(Math.min(
    RADIUS.farthest,
    Math.max(RADIUS.nearest, camera.position.length() * Math.exp(event.deltaY * 0.0012))
  ));
}, { passive: false });

// The drawn structure is normalised to unit radius, so one framing fits
// every snapshot: back off far enough that a unit ball fills the narrower
// field of view.
function frame() {
  camera.aspect = innerWidth / innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(innerWidth, innerHeight);
}

function fit() {
  const vertical = THREE.MathUtils.degToRad(camera.fov) / 2;
  const horizontal = Math.atan(Math.tan(vertical) * camera.aspect);

  camera.position.setLength(1.12 / Math.sin(Math.min(vertical, horizontal)));
}

addEventListener('resize', frame);
frame();
fit();

let previousTime = performance.now();

renderer.setAnimationLoop(() => {
  const time = performance.now();
  const seconds = Math.min((time - previousTime) / 1000, 0.1);
  previousTime = time;

  if (!dragging) {
    spin.speed += (spin.resting - spin.speed) * Math.min(1, seconds / EASE);
    swing.setFromAxisAngle(spin.axis, spin.speed * seconds);
    group.quaternion.premultiply(swing);
  }

  // Fog reads as depth only while it straddles the structure, so its band
  // follows the camera rather than sitting at fixed distances -- clamped
  // near, so it keeps working with the camera inside the figure.
  if (scene.fog) {
    const distance = camera.position.length();
    scene.fog.near = Math.max(distance - 1.05, 0.05);
    scene.fog.far = distance + 1.7;
  }

  renderer.render(scene, camera);
});
