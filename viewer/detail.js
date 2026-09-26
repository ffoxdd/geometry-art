// Every ball in a drawing is one instance of one geometry, so whatever it
// costs in triangles is paid once per instance. Detail is chosen to hold
// the total roughly fixed instead: an icosahedron subdivided as far as the
// budget allows, which is detail falling with the square root of the count.
// A subdivided icosahedron is the right ball to spend those triangles on --
// its faces are near enough equal, so nothing is wasted pinching poles.

const TRIANGLE_BUDGET = 2e6;
const ICOSAHEDRON_FACES = 20;
const FINEST = 3;
const SEGMENTS_PER_ICOSAHEDRON_EDGE = 6;

export function ballDetail(count) {
  const affordable = TRIANGLE_BUDGET / (ICOSAHEDRON_FACES * Math.max(count, 1));

  return Math.max(0, Math.min(FINEST, Math.floor(Math.log(affordable) / Math.log(4))));
}

// A tube that ends in a ball should be as round as the ball, or the ball
// shows through its flats: a great circle crosses about six of the
// icosahedron's edges, each subdivided detail + 1 times.
export function ringSegments(detail) {
  return SEGMENTS_PER_ICOSAHEDRON_EDGE * (detail + 1);
}
