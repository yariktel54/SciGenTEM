// js/bonds.js
// Robust bond guessing similar in spirit to JSmol/3Dmol:
//  - covalent radii (from RDKit) + tolerance
//  - optional PBC minimum-image if atoms._cell_vectors is provided
//  - valence/degree limits to avoid garbage (e.g. O–H–O triangles)

import { RDKit } from './rdkit_wrap.js';
import { norm_sym, fallback_symbol_to_Z, fallback_rcovalent } from './ptable_fallback.js';

// Keep for compatibility if somewhere used
export function element_symbol_to_Z(sym) {
  const n = norm_sym(sym);
  let Z = 0;
  try { Z = RDKit.ptable_getAtomicNumber(n) || 0; } catch { Z = 0; }
  if (!Z) Z = fallback_symbol_to_Z(n) || 0;
  return Z;
}

export function get_covalent_radius(Z) {
  let r = 0;
  try { r = RDKit.ptable_getRcovalent(Z); } catch { r = 0; }
  if (!Number.isFinite(r) || r <= 0) r = fallback_rcovalent(Z) || 0.77;
  return r;
}

// Heuristic max degree (not formal valence, but works well for "guess bonds")
function max_degree(Z) {
  // common organic / bio / water
  if (Z === 1) return 1;   // H
  if (Z === 5) return 3;   // B
  if (Z === 6) return 4;   // C
  if (Z === 7) return 3;   // N
  if (Z === 8) return 2;   // O
  if (Z === 9) return 1;   // F
  if (Z === 14) return 4;  // Si
  if (Z === 15) return 5;  // P
  if (Z === 16) return 6;  // S
  if (Z === 17) return 1;  // Cl
  if (Z === 35) return 1;  // Br
  if (Z === 53) return 1;  // I

  // alkali metals
  if (Z === 3 || Z === 11 || Z === 19 || Z === 37 || Z === 55 || Z === 87) return 1;
  // alkaline earth
  if (Z === 4 || Z === 12 || Z === 20 || Z === 38 || Z === 56 || Z === 88) return 2;

  // transition metals: allow more neighbours
  if ((Z >= 21 && Z <= 30) || (Z >= 39 && Z <= 48) || (Z >= 72 && Z <= 80)) return 8;
  // lanthanides / actinides
  if ((Z >= 57 && Z <= 71) || (Z >= 89 && Z <= 103)) return 10;

  // default
  return 4;
}

function dot(a, b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
function cross(a, b) {
  return [
    a[1]*b[2] - a[2]*b[1],
    a[2]*b[0] - a[0]*b[2],
    a[0]*b[1] - a[1]*b[0]
  ];
}

function pair_max_abs(Zi, Zj, global_cap) {
  // Extra conservative caps for light elements to avoid counting H-bonds / nonbonded contacts as covalent.
  if (Zi === 1 || Zj === 1) return Math.min(global_cap, 1.25); // H–X
  if ((Zi === 8 && Zj === 8) || (Zi === 7 && Zj === 8) || (Zi === 6 && Zj === 8) || (Zi === 6 && Zj === 7)) {
    return Math.min(global_cap, 2.0);
  }
  return global_cap;
}

function build_pbc_tools(cell_vectors) {
  if (!cell_vectors || cell_vectors.length !== 3) return null;
  const a = cell_vectors[0], b = cell_vectors[1], c = cell_vectors[2];
  if (!a || !b || !c) return null;
  if (![...a, ...b, ...c].every(Number.isFinite)) return null;

  const bxc = cross(b, c);
  const cxa = cross(c, a);
  const axb = cross(a, b);
  const det = dot(a, bxc);

  if (!Number.isFinite(det) || Math.abs(det) < 1e-10) return null;

  // invM rows: (b×c)/det, (c×a)/det, (a×b)/det
  const inv = [
    [bxc[0]/det, bxc[1]/det, bxc[2]/det],
    [cxa[0]/det, cxa[1]/det, cxa[2]/det],
    [axb[0]/det, axb[1]/det, axb[2]/det],
  ];

  function cart_to_frac(v) {
    return [
      inv[0][0]*v[0] + inv[0][1]*v[1] + inv[0][2]*v[2],
      inv[1][0]*v[0] + inv[1][1]*v[1] + inv[1][2]*v[2],
      inv[2][0]*v[0] + inv[2][1]*v[1] + inv[2][2]*v[2],
    ];
  }

  function frac_to_cart(f) {
    return [
      f[0]*a[0] + f[1]*b[0] + f[2]*c[0],
      f[0]*a[1] + f[1]*b[1] + f[2]*c[1],
      f[0]*a[2] + f[1]*b[2] + f[2]*c[2],
    ];
  }

  function min_image_vec(dr) {
    const f = cart_to_frac(dr);
    // wrap to [-0.5, 0.5] in fractional
    f[0] -= Math.round(f[0]);
    f[1] -= Math.round(f[1]);
    f[2] -= Math.round(f[2]);
    return frac_to_cart(f);
  }

  return { min_image_vec };
}

export function guess_bonds_by_distance(atoms, scale_factor = 1.25, max_abs = 2.3) {
  if (!Array.isArray(atoms) || atoms.length < 2) return [];

  const usePBC = !!atoms._pbc && Array.isArray(atoms._cell_vectors);
  const pbc = usePBC ? build_pbc_tools(atoms._cell_vectors) : null;

  const n = atoms.length;

  const deg = new Array(n).fill(0);
  const maxDeg = new Array(n);
  const rad = new Array(n);

  for (let i = 0; i < n; i++) {
    const Zi = atoms[i]?.Z ?? 6;
    maxDeg[i] = max_degree(Zi);
    rad[i] = get_covalent_radius(Zi);
  }

  // collect candidates
  const cand = [];
  for (let i = 0; i < n; i++) {
    const ai = atoms[i];
    const Zi = ai?.Z ?? 6;

    for (let j = i + 1; j < n; j++) {
      const aj = atoms[j];
      const Zj = aj?.Z ?? 6;

      // never H-H
      if (Zi === 1 && Zj === 1) continue;

      let dx = (aj.x - ai.x);
      let dy = (aj.y - ai.y);
      let dz = (aj.z - ai.z);

      if (![dx, dy, dz].every(Number.isFinite)) continue;

      if (pbc) {
        const v = pbc.min_image_vec([dx, dy, dz]);
        dx = v[0]; dy = v[1]; dz = v[2];
      }

      const d2 = dx*dx + dy*dy + dz*dz;
      if (!Number.isFinite(d2) || d2 < 1e-6) continue;

      const d = Math.sqrt(d2);
      const rsum = rad[i] + rad[j];

      // threshold: scaled radii + small padding, with sane cap
      const pad = (Zi === 1 || Zj === 1) ? 0.15 : 0.2;
      let thr = scale_factor * rsum + pad;
      thr = Math.min(thr, pair_max_abs(Zi, Zj, max_abs));
      if (thr < 0.6) thr = 0.6;

      if (d <= thr) {
        // score: how "tight" relative to radii sum
        const score = d / Math.max(1e-6, rsum);
        cand.push({ i, j, d, score });
      }
    }
  }

  // sort best (shortest/most plausible) first
  cand.sort((a, b) => (a.score - b.score) || (a.d - b.d));

  // greedy accept with degree limits
  const bonds = [];
  const seen = new Set();

  for (const c of cand) {
    const i = c.i, j = c.j;

    // degree limits
    if (deg[i] >= maxDeg[i]) continue;
    if (deg[j] >= maxDeg[j]) continue;

    // extra hard cap for H: only 1 bond ever
    if ((atoms[i].Z === 1 && deg[i] >= 1) || (atoms[j].Z === 1 && deg[j] >= 1)) continue;

    const key = i < j ? (i + ',' + j) : (j + ',' + i);
    if (seen.has(key)) continue;

    seen.add(key);
    bonds.push([i, j, 1]);
    deg[i]++; deg[j]++;
  }

  return bonds;
}
