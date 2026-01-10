// js/phases.js
// System builder: detects input format (SMILES / CIF / CSV) and returns atoms + bonds.
// NOTE: "liquid/solid" modes were removed. We always build what the last input represents.
//
// Safe-refactor wave A:
//  - phases.js is a thin dispatcher
//  - CIF tiling math moved to lattice.js
//  - per-format building moved to builders/

import { build_from_smiles } from './builders/build_smiles.js';
import { build_from_cif } from './builders/build_cif.js';
import { build_from_csv } from './builders/build_csv.js';
import { is_dev_mode, normalize_build_result, validate_system } from './system_contract.js';

export function parse_system_string(s) {
    s = (s || '').trim();
    const low = s.toLowerCase();

    // ---------- helpers ----------
    const parse_size = (val) => {
        val = (val || '').trim().toLowerCase();
        if (val === 'auto') return 'auto';
        const parts = val.split(/[x×]/);
        if (parts.length !== 3) return [1, 1, 1];
        const nx = parseInt(parts[0], 10);
        const ny = parseInt(parts[1], 10);
        const nz = parseInt(parts[2], 10);
        if ([nx, ny, nz].some(Number.isNaN)) return [1, 1, 1];
        return [nx, ny, nz];
    };

    const parse_kv = (body) => {
        const items = body.replace(/,/g, ';').split(';').map(x => x.trim()).filter(Boolean);
        const kv = {};
        for (const it of items) {
            if (it.includes('=')) {
                const [k, v] = it.split('=', 2);
                kv[k.trim().toLowerCase()] = (v ?? '').trim();
            }
        }
        return kv;
    };

    const looks_like_csv_text = (txt) => {
        // Strict CSV uses section header "ATOMS" and optional "BONDS".
        return /(^|\n)\s*ATOMS\s*(\n|$)/i.test(txt);
    };

    const looks_like_cif_text = (txt) => {
        // CIF usually contains data_ and/or atom_site/cell fields.
        if (/(^|\n)\s*data_/i.test(txt)) return true;
        if (/(^|\n)\s*_atom_site_/i.test(txt)) return true;
        if (/(^|\n)\s*_cell_length_/i.test(txt)) return true;
        return false;
    };

    // -------- format registry: add new formats by appending handlers --------
    // Each handler: { name, detect(str, low, helpers) -> spec|null }

    const helpers = { parse_size, parse_kv, looks_like_csv_text, looks_like_cif_text };

    const detect_cif = (str, lowStr, h) => {
        // CIF with optional args ("...cif;size=2x2x1")
        // Works for normal URLs/paths AND for blob urls with fragments: "blob:...#file.cif;size=...".
        if (lowStr.includes('.cif;')) {
            const [path, rest] = str.split(';', 2);
            let size = [1, 1, 1];
            const kv = h.parse_kv(rest);
            if ('size' in kv) size = h.parse_size(kv.size);
            return { mode: 'cif', path: path.trim(), size };
        }

        // Plain .cif (also matches "...#file.cif")
        if (lowStr.endsWith('.cif') || /\.cif([?#].*)?$/i.test(str)) {
            return { mode: 'cif', path: str, size: [1, 1, 1] };
        }

        // Raw pasted CIF (robustness)
        if (h.looks_like_cif_text(str)) return { mode: 'cif', path: str, size: [1, 1, 1] };
        return null;
    };

    const detect_csv = (str, lowStr, h) => {
        // Plain .csv (also matches "...#file.csv")
        if (lowStr.endsWith('.csv') || /\.csv([?#].*)?$/i.test(str)) {
            return { mode: 'csv', path: str };
        }

        // Raw pasted strict CSV (robustness)
        if (h.looks_like_csv_text(str)) return { mode: 'csv', path: str };
        return null;
    };

    const detect_smiles = (str) => {
        // Default fallback.
        return { mode: 'smiles', smiles: str };
    };

    const FORMAT_HANDLERS = [
        { name: 'cif', detect: detect_cif },
        { name: 'csv', detect: detect_csv },
        { name: 'smiles', detect: (str, lowStr) => detect_smiles(str, lowStr) },
    ];

    for (const hnd of FORMAT_HANDLERS) {
        const spec = hnd.detect(s, low, helpers);
        if (spec) return spec;
    }

    // should never happen due to SMILES fallback
    return { mode: 'smiles', smiles: s };
}

export async function build_system_from_input(input_str) {
    const spec = parse_system_string(input_str);

    // Builder registry: add new formats here (one line) and you're done.
    const BUILDERS = {
        cif: build_from_cif,
        csv: build_from_csv,
        smiles: build_from_smiles,
    };

    const build_fn = BUILDERS[spec.mode];
    if (!build_fn) throw new Error('Непідтримуваний режим');

    const warnings = [];
    let res = await build_fn(spec);

    // Merge builder-provided warnings (if any).
    try {
        const mw = res?.meta?.warnings;
        if (Array.isArray(mw)) warnings.push(...mw);
    } catch (_) {}


    const sys = normalize_build_result(res);
    validate_system(sys, `build:${spec.mode}`);

    // Lightweight warnings (no console spam). UI may choose to show them.
    if (spec.mode === 'cif') {
        if (Array.isArray(spec.size) && (spec.size[0] !== 1 || spec.size[1] !== 1 || spec.size[2] !== 1)) {
            warnings.push(`CIF tiling: ${spec.size[0]}x${spec.size[1]}x${spec.size[2]}`);
        }
        if (sys.bonds === null) {
            warnings.push('CIF: зв’язки не визначені (bonds=null)');
        } else if (Array.isArray(sys.bonds) && sys.bonds.length === 0) {
            warnings.push('CIF: зв’язки не знайдено (спробуйте збільшити size)');
        }
    }
    if (spec.mode === 'csv') {
        if (Array.isArray(sys.bonds) && sys.bonds.length === 0) {
            warnings.push('CSV: bonds порожні (BONDS секція відсутня або порожня)');
        }
    }

    // Provide metadata without changing the public return shape.
    const out = [sys.atoms, sys.bonds, sys.title];
    out.meta = { mode: spec.mode, warnings: Array.from(new Set(warnings)).slice(0, 8) };

    // Optional DEV-only single-line debug summary.
    if (is_dev_mode()) {
        const nb = (sys.bonds === null) ? 'null' : String(sys.bonds.length);
        // keep it short; no spam
        console.debug(`[TEM] build ${spec.mode}: atoms=${sys.atoms.length}; bonds=${nb}`);
        if (warnings.length) console.debug(`[TEM] warnings: ${warnings.join(' | ')}`);
    }

    // Also expose last warnings for quick UI/debug access.
    try {
        globalThis.TEM_LAST_BUILD = { mode: spec.mode, warnings, atoms: sys.atoms.length, bonds: sys.bonds === null ? null : sys.bonds.length, title: sys.title };
    } catch (_) {}

    return out;
}
