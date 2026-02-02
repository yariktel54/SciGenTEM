// js/renderer.js
// FIXES:
// 1) прибрано «спам» console.log
// 2) виправлено виклики draw_atoms/draw_bonds у DoF-гілці (раніше параметри передавались у неправильному порядку)
// 3) draw_bonds спрощено (background_gray/scale були зайвими і не використовувались)
// 4) hide_front тепер реально працює (обрізає атоми/шари «перед» фокусом)

import { get_covalent_radius, guess_bonds_by_distance } from './bonds.js';

export function compute_scaled_coordinates(atoms, img_size, angstroms_per_pixel) {
    const xs = atoms.map(a => a.x);
    const ys = atoms.map(a => a.y);
    const meanX = xs.reduce((t, v) => t + v, 0) / Math.max(xs.length, 1);
    const meanY = ys.reduce((t, v) => t + v, 0) / Math.max(ys.length, 1);
    const scale = 1.0 / Math.max(angstroms_per_pixel, 1e-9);
    const cx = (img_size[1] >> 1), cy = (img_size[0] >> 1);
    const coords = xs.map((_, i) => [
        Math.round((xs[i] - meanX) * scale + cx),
        Math.round((ys[i] - meanY) * scale + cy)
    ]);
    return [coords, scale];
}

// «полотно» — масив Float32 (H×W), далі конвертуємо у 8-біт для Canvas
function newFloatImage(h, w, bg = 127) {
    const a = new Float32Array(h * w);
    a.fill(bg);
    return { a, h, w, bg };
}
function idx(img, x, y) { return y * img.w + x; }
function clamp255(v) { return v < 0 ? 0 : v > 255 ? 255 : v; }

function gaussianBlurFloat(img, sigma) {
    if (sigma <= 0) return img;
    const { h, w, a } = img;

    const out = new Float32Array(a);

    const ksize = Math.max(1, Math.floor(sigma * 3) * 2 + 1);
    const half = ksize >> 1;
    const kernel = new Float32Array(ksize);
    let sum = 0;
    for (let i = -half; i <= half; i++) {
        const v = Math.exp(-(i * i) / (2 * sigma * sigma));
        kernel[i + half] = v;
        sum += v;
    }
    for (let i = 0; i < ksize; i++) kernel[i] /= sum;

    // X
    for (let y = 0; y < h; y++) {
        for (let x = 0; x < w; x++) {
            let s = 0;
            for (let k = -half; k <= half; k++) {
                const xx = Math.min(w - 1, Math.max(0, x + k));
                s += out[y * w + xx] * kernel[k + half];
            }
            a[y * w + x] = s;
        }
    }

    // Y
    out.set(a);
    for (let y = 0; y < h; y++) {
        for (let x = 0; x < w; x++) {
            let s = 0;
            for (let k = -half; k <= half; k++) {
                const yy = Math.min(h - 1, Math.max(0, y + k));
                s += out[yy * w + x] * kernel[k + half];
            }
            a[y * w + x] = s;
        }
    }

    return img;
}

// opts: {compose, focal_z, dof_strength, hide_front}
export function draw_atoms(img, atoms, coords, scale, background_gray, opts = {}) {
    const {
        compose = 'sum',
        focal_z = 0.0,
        dof_strength = 0.0,
        hide_front = false
    } = opts;

    // (compose/dof_strength поки не використовуються напряму — залишено для сумісності)
    void compose; void dof_strength;

    const Zs = atoms.map(a => a.Z | 0);
    const radiiA = Zs.map(Z => get_covalent_radius(Z));

    for (let i = 0; i < atoms.length; i++) {
        const [x0, y0] = coords[i];

        if (hide_front && (atoms[i].z ?? 0) > focal_z + 1e-9) continue;

        // Реалістична "геометрична" ширина атома: залежить від ковалентного радіуса
        // atom_size_mul — загальний множник, atom_size_exp — підсилює різницю між H і важкими
        const atom_size_mul = opts.atom_size_mul ?? 1.2;     // 0.6..1.6 (тюнінг)
        const atom_size_exp = opts.atom_size_exp ?? 1.25;    // 1.0..1.5 (1.25 дає помітну різницю)
        const rA = radiiA[i];

        // без "мінімумів": лише числовий захист від нуля
        const sigma = Math.max(1e-6, Math.pow(rA, atom_size_exp) * scale * atom_size_mul);


        // Інтенсивність
        const atom_dark_mul = opts.atom_dark_mul ?? 2.5;   // 1.5..4.0
        const atom_dark_exp = opts.atom_dark_exp ?? 1.2;   // 1.1..1.6

        const baseZ = Math.pow(Zs[i], atom_dark_exp);
        const intensity = -(baseZ * atom_dark_mul);

        // «ядро»
        const core_sigma_rel = opts.core_sigma_rel ?? 0.18; // 0.12..0.22
        const core_rel = opts.core_rel ?? 0.6;             // 0.3..0.9
        const core_Z0 = opts.core_Z0 ?? 12;                // 8..20 (чим менше — тим сильніше "гасить" ядро для металів)

        const core_sigma = sigma * core_sigma_rel;

        // ядро не масштабуємо напряму як intensity, а робимо "контрольованим" і
        // трохи приглушуємо для великих Z, щоб важкі атоми не виглядали як "однакові круги з білим центром"
        const core_falloff = 1.0 / (1.0 + (Zs[i] / core_Z0) * (Zs[i] / core_Z0));
        const core_intensity = baseZ * core_rel * core_falloff;


        // --- Gaussian "тіло" атома ---
        const R = Math.ceil(3 * sigma);
        for (let dx = -R; dx <= R; dx++) {
            for (let dy = -R; dy <= R; dy++) {
                const x = x0 + dx;
                const y = y0 + dy;
                if (x < 0 || y < 0 || y >= img.h || x >= img.w) continue;
                const dist2 = dx * dx + dy * dy;
                const v = intensity * Math.exp(-dist2 / (2 * sigma * sigma));
                const p = (y * img.w + x);
                const current = img.a[p];
                const candidate = background_gray + v;
                img.a[p] = Math.min(current, candidate);
            }
        }

        // --- «ядро» (яскравий центр) ---
        const Rc = Math.ceil(3 * core_sigma);
        for (let dx = -Rc; dx <= Rc; dx++) {
            for (let dy = -Rc; dy <= Rc; dy++) {
                const x = x0 + dx;
                const y = y0 + dy;
                if (x < 0 || y < 0 || y >= img.h || x >= img.w) continue;
                const dist2 = dx * dx + dy * dy;
                const add = core_intensity * Math.exp(-dist2 / (2 * core_sigma * core_sigma));
                const p = (y * img.w + x);
                img.a[p] += add;
            }
        }
    }
}

function draw_gaussian_line(img, p1, p2, intensity = 30, halfwidth_px = 4.0) {
    const [x1, y1] = p1.map(Math.round), [x2, y2] = p2.map(Math.round);
    let dx = x2 - x1, dy = y2 - y1;
    const L = Math.hypot(dx, dy) | 0;
    if (L <= 0) return;

    dx /= L; dy /= L;
    const nx = -dy, ny = dx;
    const radius = Math.ceil(halfwidth_px);

    for (let i = 0; i <= L; i++) {
        const x = Math.round(x1 + dx * i), y = Math.round(y1 + dy * i);
        for (let sx = -radius; sx <= radius; sx++) {
            for (let sy = -radius; sy <= radius; sy++) {
                const xx = x + sx, yy = y + sy;
                if (xx < 0 || xx >= img.w || yy < 0 || yy >= img.h) continue;
                const d_perp = Math.abs(sx * nx + sy * ny);
                const s = d_perp / Math.max(halfwidth_px, 1e-6);
                if (s <= 1.0) {
                    const w = Math.cos(0.5 * Math.PI * s);
                    const val = intensity * (w * w);
                    if (val > 0.25) img.a[idx(img, xx, yy)] += val;
                }
            }
        }
    }
}

function draw_bond_glow_line(img, p1, p2, intensity = 50, sigma = 1.2) {
    const [x1, y1] = p1.map(Math.round), [x2, y2] = p2.map(Math.round);
    let dx = x2 - x1, dy = y2 - y1;
    const L = Math.hypot(dx, dy) | 0;
    if (L === 0) return;

    for (let i = 0; i <= L; i++) {
        const x = Math.round(x1 + (dx * i) / L), y = Math.round(y1 + (dy * i) / L);
        const R = Math.max(1, Math.floor(3 * sigma));
        for (let sx = -R; sx <= R; sx++) {
            for (let sy = -R; sy <= R; sy++) {
                const xx = x + sx, yy = y + sy;
                if (xx < 0 || xx >= img.w || yy < 0 || yy >= img.h) continue;
                const dist2 = sx * sx + sy * sy;
                const value = intensity * Math.exp(-dist2 / (2 * sigma * sigma));
                if (value > 0.5) img.a[idx(img, xx, yy)] += value;
            }
        }
    }
}

// opts: {wave_width_px, wave_amp}
export function draw_bonds(img, coords, bonds, atoms, opts = {}) {
    const { wave_width_px = 8.0, wave_amp = 1.0 } = opts;
    const halfw = Math.max(1.0, wave_width_px * 0.5);

    for (const b of bonds) {
        const [i, j, bt] = b;
        const p1 = coords[i], p2 = coords[j];
        if (!p1 || !p2) continue;

        const Z1 = atoms[i]?.Z ?? 0;
        const Z2 = atoms[j]?.Z ?? 0;
        const Zavg = 0.5 * (Z1 + Z2);

        const base_int = Math.min(255.0, Math.pow(Zavg, 0.9));
        const glow_int = Math.min(255.0, Math.pow(Zavg, 1.2) * 2.0);
        const bond_intensity = Math.min(255.0, base_int * wave_amp);
        const glow_intensity = Math.min(255.0, glow_int * wave_amp);

        // зсуви для 1x/2x/3x
        const direction = [(p2[0] - p1[0]), (p2[1] - p1[1])];
        const L = Math.hypot(direction[0], direction[1]) || 1;
        const nrm = [-direction[1] / L, direction[0] / L];
        const offset_base = Math.max(1.0, 0.8 * halfw);

        let offsets = [[0, 0]];
        if (bt === 2) offsets = [[nrm[0] * offset_base, nrm[1] * offset_base], [-nrm[0] * offset_base, -nrm[1] * offset_base]];
        else if (bt === 3) offsets = [[nrm[0] * 1.2 * offset_base, nrm[1] * 1.2 * offset_base], [0, 0], [-nrm[0] * 1.2 * offset_base, -nrm[1] * 1.2 * offset_base]];

        for (const off of offsets) {
            const q1 = [p1[0] + off[0], p1[1] + off[1]];
            const q2 = [p2[0] + off[0], p2[1] + off[1]];
            draw_gaussian_line(img, q1, q2, bond_intensity, halfw);

            // glow тільки для центральної лінії
            if (off[0] === 0 && off[1] === 0) {
                draw_bond_glow_line(img, p1, p2, glow_intensity * 0.85, Math.max(0.5, 0.6 * halfw));
            }
        }
    }
}

function nice_scale_length(angstroms_per_pixel, width_px, min_px = 80, max_px = 200) {
    if (width_px <= 0 || angstroms_per_pixel <= 0) return [0, 0];
    const target_px = Math.min(Math.max(width_px * 0.25, min_px), max_px);
    const target_A = target_px * angstroms_per_pixel;
    const k = target_A > 0 ? Math.floor(Math.log10(target_A)) : 0;

    let best = null;
    for (let shift = -3; shift <= 3; shift++) {
        const base = Math.pow(10, k + shift);
        for (const m of [1, 2, 5]) {
            const valA = m * base;
            const valPx = valA / angstroms_per_pixel;
            if (valPx >= min_px && valPx <= max_px) {
                if (!best || Math.abs(valPx - target_px) < Math.abs(best[0] - target_px)) best = [valPx, valA];
            }
        }
    }
    if (!best) {
        const valPx = max_px, valA = valPx * angstroms_per_pixel;
        return [Math.round(valPx), valA];
    }
    return [Math.round(best[0]), best[1]];
}

function draw_scale_bar(ctx, img, angstroms_per_pixel, { corner = 'bl', margin = 12, invert = false } = {}) {
    const h = img.h, w = img.w;
    const [bar_px, bar_A] = nice_scale_length(angstroms_per_pixel, w, 80, 200);
    if (bar_px <= 0) return;

    const col = invert ? 0 : 255;

    const thickness = Math.max(4, Math.min(18, (w / 80) | 0));
    const tick_w = Math.max(2, thickness - 2);
    const tick_h = (thickness * 1.8) | 0;
    const gap_txt = Math.max(8, (thickness >> 1) + 6);
    const pad_edge = Math.max(20, (margin * 1.6) | 0);

    const text = (Math.abs(bar_A - Math.round(bar_A)) < 1e-6)
        ? `${Math.round(bar_A)} A`
        : `${bar_A.toFixed(1)} A`;

    let x0, x1, y;
    const right = /r$/.test(corner);
    const top = /^t/.test(corner);

    if (right) {
        x1 = w - margin;
        x0 = x1 - bar_px;
    } else {
        x0 = margin;
        x1 = margin + bar_px;
    }

    y = top ? pad_edge : (h - pad_edge);

    ctx.fillStyle = `rgb(${col},${col},${col})`;

    ctx.fillRect(x0, y - (thickness >> 1), bar_px, thickness);
    ctx.fillRect(x0, y - (tick_h >> 1), tick_w, tick_h);
    ctx.fillRect(x1 - tick_w, y - (tick_h >> 1), tick_w, tick_h);

    ctx.font = `${Math.max(14, (thickness * 2) | 0)}px system-ui, Arial, sans-serif`;
    ctx.textBaseline = 'top';

    const tw = ctx.measureText(text).width | 0;
    const tx = right ? (x1 - tw) : x0;
    const ty = top
        ? (y + (thickness >> 1) + gap_txt)
        : (y - (thickness >> 1) - gap_txt - (thickness * 2));

    ctx.fillText(text, tx, Math.max(4, Math.min(h - 16, ty)));
}

export function render_image(
    atoms,
    {
        bonds = null,
        img_size = [400, 400],
        angstroms_per_pixel = 0.1,
        blur_sigma = 2,
        background_gray = 127,
        invert = false,
        noise_stddev = 2.0,
        contrast = 1.0,
        compose_mode = 'sum',
        draw_bonds_flag = true,
        bond_wave_width_px = 8.0,
        bond_wave_amplitude = 1.0,
        low_clip = null,
        high_clip = null,
        focal_z = 0.0,
        dof_strength = 0.0,
        hide_front = false,
        show_scale_bar = false,
        scale_bar_corner = 'bl',
        scale_bar_margin_px = 12,
        // додатково: контекст Canvas для малювання масштабної лінійки поверх
        canvasCtx = null
    } = {}
) {
    const [H, W] = img_size;
    const img = newFloatImage(H, W, background_gray);

    if (!atoms || !atoms.length) {
        // навіть якщо пусто — очистимо канву
        if (canvasCtx) {
            const imageData = canvasCtx.createImageData(W, H);
            for (let p = 0; p < imageData.data.length; p += 4) {
                imageData.data[p] = background_gray;
                imageData.data[p + 1] = background_gray;
                imageData.data[p + 2] = background_gray;
                imageData.data[p + 3] = 255;
            }
            canvasCtx.putImageData(imageData, 0, 0);
        }
        return new Uint8ClampedArray(H * W);
    }

    const [coords, scale] = compute_scaled_coordinates(atoms, img_size, angstroms_per_pixel);

    // Семантика bonds:
    //  • bonds === null або undefined → зв'язки невідомі → можна здогадувати по відстанях;
    //  • bonds — масив (навіть порожній [] ) → це явна топологія, НІЧОГО не "guess"-ити.
    let use_bonds = [];
    if (draw_bonds_flag) {
        use_bonds = (bonds == null) ? guess_bonds_by_distance(atoms) : bonds;
    }

    const dof_on = (dof_strength && dof_strength > 1e-6);

    if (dof_on) {
        const zvals = atoms.map(a => a.z ?? 0);
        let zmin = 0, zmax = 0;
        if (zvals.length) { zmin = Math.min(...zvals); zmax = Math.max(...zvals); }

        const K = (zmax > zmin) ? 12 : 1;
        const edges = Array.from({ length: K + 1 }, (_, i) => zmin + (zmax - zmin) * i / K);
        const centers = Array.from({ length: K }, (_, i) => 0.5 * (edges[i] + edges[i + 1]));

        const bond_z = (use_bonds.length)
            ? use_bonds.map(([i, j]) => 0.5 * ((atoms[i].z ?? 0) + (atoms[j].z ?? 0)))
            : null;

        for (let bi = 0; bi < K; bi++) {
            const z_lo = edges[bi], z_hi = edges[bi + 1], z_c = centers[bi];
            if (hide_front && z_c > focal_z + 1e-9) continue;

            // індекси атомів шару (глобальні)
            const idxs = [];
            for (let i = 0; i < atoms.length; i++) {
                const z = atoms[i].z ?? 0;
                if (z >= z_lo && (z < z_hi || (bi === K - 1 && z <= z_hi))) idxs.push(i);
            }

            if (!idxs.length && (!use_bonds.length || !bond_z)) continue;

            // мапа глобальний → локальний
            const indexMap = new Map();
            for (let li = 0; li < idxs.length; li++) indexMap.set(idxs[li], li);

            // зв’язки шару (локальні індекси)
            const bonds_sub = [];
            if (draw_bonds_flag && use_bonds.length && bond_z) {
                for (let k = 0; k < use_bonds.length; k++) {
                    const zb = bond_z[k];
                    if (!(zb >= z_lo && (zb < z_hi || (bi === K - 1 && zb <= z_hi)))) continue;
                    const [gi, gj, bt] = use_bonds[k];
                    const li = indexMap.get(gi);
                    const lj = indexMap.get(gj);
                    if (li === undefined || lj === undefined) continue;
                    bonds_sub.push([li, lj, bt]);
                }
            }

            if (!idxs.length && !bonds_sub.length) continue;

            const atoms_sub = idxs.map(i => atoms[i]);
            const coords_sub = idxs.map(i => coords[i]);

            const layer = newFloatImage(H, W, background_gray);
            draw_atoms(layer, atoms_sub, coords_sub, scale, background_gray, { compose: compose_mode, focal_z: 0, dof_strength: 0, hide_front: false });
            if (draw_bonds_flag && bonds_sub.length) {
                draw_bonds(layer, coords_sub, bonds_sub, atoms_sub, { wave_width_px: bond_wave_width_px, wave_amp: bond_wave_amplitude });
            }

            // delta (layer - bg) -> tmp
            const tmp = { a: new Float32Array(H * W), h: H, w: W, bg: 0 };
            for (let p = 0; p < tmp.a.length; p++) tmp.a[p] = layer.a[p] - background_gray;

            // DoF блюр у пікселях
            const sigma_px = Math.abs(z_c - focal_z) * dof_strength / Math.max(angstroms_per_pixel, 1e-9);
            if (sigma_px > 0) gaussianBlurFloat(tmp, sigma_px);

            for (let p = 0; p < tmp.a.length; p++) img.a[p] += tmp.a[p];
        }
    } else {
        draw_atoms(img, atoms, coords, scale, background_gray, { compose: compose_mode, focal_z, dof_strength: 0, hide_front });
        if (draw_bonds_flag && use_bonds.length) {
            draw_bonds(img, coords, use_bonds, atoms, { wave_width_px: bond_wave_width_px, wave_amp: bond_wave_amplitude });
        }
    }

    if (blur_sigma > 0) gaussianBlurFloat(img, blur_sigma);

    if (noise_stddev > 0) {
        const σ = noise_stddev;
        for (let i = 0; i < img.a.length; i++) {
            // Box–Muller
            const u1 = Math.random(), u2 = Math.random();
            const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
            img.a[i] += z * σ;
        }
    }

    if (contrast !== 1.0) {
        for (let i = 0; i < img.a.length; i++) {
            img.a[i] = (img.a[i] - background_gray) * contrast + background_gray;
        }
    }

    if (invert) {
        for (let i = 0; i < img.a.length; i++) {
            const d = background_gray - img.a[i];
            img.a[i] = background_gray + d;
        }
    }

    if (low_clip != null || high_clip != null) {
        const lo = (low_clip == null ? 0 : low_clip) | 0;
        const hi = (high_clip == null ? 255 : high_clip) | 0;
        const a = Math.min(lo, hi), b = Math.max(lo, hi);
        for (let i = 0; i < img.a.length; i++) {
            const v = img.a[i];
            img.a[i] = v < a ? a : v > b ? b : v;
        }
    }

    // до 8-біт
    const out = new Uint8ClampedArray(H * W);
    for (let i = 0; i < out.length; i++) out[i] = clamp255(img.a[i]);

    // Рисуємо на Canvas, якщо переданий контекст
    if (canvasCtx) {
        const imageData = canvasCtx.createImageData(W, H);
        // grayscale → RGBA
        for (let i = 0, p = 0; i < out.length; i++, p += 4) {
            const v = out[i];
            imageData.data[p] = v;
            imageData.data[p + 1] = v;
            imageData.data[p + 2] = v;
            imageData.data[p + 3] = 255;
        }
        canvasCtx.putImageData(imageData, 0, 0);
        if (show_scale_bar) {
            draw_scale_bar(canvasCtx, { h: H, w: W }, angstroms_per_pixel, {
                corner: scale_bar_corner,
                margin: scale_bar_margin_px,
                invert
            });
        }
    }

    return out;
}
