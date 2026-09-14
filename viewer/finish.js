// How a surface answers light, shared by every drawing that offers the
// choice: a part's colour says what it is, its finish says what it is
// made of.

export const CHOICES = { satin: 'satin', matte: 'matte', metal: 'metal' };

const FINISHES = {
  satin: { metalness: 0.15, roughness: 0.55 },
  matte: { metalness: 0, roughness: 0.8 },
  metal: { metalness: 1, roughness: 0.24 },
};

export function of(name) {
  return FINISHES[name] ?? FINISHES.satin;
}
