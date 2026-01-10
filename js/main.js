// js/main.js
// UI + state + rebuild/render pipeline (simplified)
//
// IMPORTANT:
//  - main.js MUST NOT construct ImageData from a grayscale buffer.
//    renderer.js already converts grayscale -> RGBA and draws into canvas.
//
// UI (current):
//  - CIF/CSV are loaded ONLY via file input / drag&drop.
//  - No "liquid/solid" modes.

import { build_system_from_input } from './phases.js';
import { render_image } from './renderer.js';

export async function interactive_em_image(smiles_text, _unused) {
    function el(id) { return document.getElementById(id); }

    // ---- DOM ----
    var cvs = el('canvas');
    if (!cvs) throw new Error('Canvas #canvas not found');
    var ctx = cvs.getContext('2d', { willReadFrequently: true });

    var lblTitle = el('lbl-title');

    // Controls
    var rngZoom = el('rng-zoom');
    var rngContrast = el('rng-contrast');
    var rngBlur = el('rng-blur');
    var rngNoise = el('rng-noise');
    var rngFocus = el('rng-focus');
    var rngDof = el('rng-dof');
    var rngBwidth = el('rng-bwidth');
    var rngBamp = el('rng-bamp');

    var cbInvert = el('cb-invert');
    var cbNoise = el('cb-noise');
    var cbBonds = el('cb-bonds');
    var cbHide = el('cb-hidefront');
    var cbScale = el('cb-scale');

    var tbTiles = el('tb-tiles');
    var tbSmiles = el('tb-smiles');
    var tbW = el('tb-w');
    var tbH = el('tb-h');

    var fileInput = el('file-input');

    // Buttons
    var btnExport = el('btn-export');

    // ---- State ----
    var atoms = [];
    var bonds = null; // null => unknown (guess allowed), Array => explicit (guess forbidden)
    var title = '';

    // last chosen source
    var active_source = 'smiles'; // 'smiles' | 'file'

    // currently loaded file (if any)
    var file_state = {
        kind: null,  // 'cif' | 'csv'
        name: '',
        url: null,
        isBlob: false
    };

    // ---- utils ----
    function trimStr(s) { return (s || '').trim(); }

    function clampInt(n, lo, hi, fallback) {
        var v = parseInt(n, 10);
        if (!Number.isFinite(v)) v = fallback;
        v = Math.min(hi, Math.max(lo, v | 0));
        return v;
    }

    function clearCanvas(bg) {
        ctx.save();
        ctx.fillStyle = 'rgb(' + bg + ',' + bg + ',' + bg + ')';
        ctx.fillRect(0, 0, cvs.width, cvs.height);
        ctx.restore();
    }

    function setTitle(t) {
        title = t || '';
        if (lblTitle) lblTitle.textContent = 'Simulated EM — ' + title;
    }

    function fileKindFromName(name) {
        var low = String(name || '').toLowerCase();
        // tolerate things like "file.cif;..." or "file.cif?x" in the #fragment
        if (low.indexOf('.cif') >= 0) return 'cif';
        if (low.indexOf('.csv') >= 0) return 'csv';
        return null;
    }

    function revokeFileUrl() {
        if (file_state.isBlob && file_state.url) {
            try { URL.revokeObjectURL(file_state.url); } catch (e) { }
        }
        file_state.kind = null;
        file_state.name = '';
        file_state.url = null;
        file_state.isBlob = false;
    }

    function parseTilesSize() {
        if (!tbTiles) return null;
        var t = trimStr(tbTiles.value);
        if (!t) return null;
        if (t.toLowerCase() === 'auto') return null;
        return t.replace(/×/g, 'x');
    }

    function applyCifSizeToSpec(spec, sizeStr) {
        if (!sizeStr) return spec;

        var s = String(spec || '').trim();
        if (!s) return s;

        var semi = s.indexOf(';');
        if (semi < 0) return s + ';size=' + sizeStr;

        var base = s.slice(0, semi).trim();
        var rest = s.slice(semi + 1).trim();
        var items = rest ? rest.split(';') : [];

        var kept = [];
        for (var i = 0; i < items.length; i++) {
            var it = trimStr(items[i]);
            if (!it) continue;
            var itLow = it.toLowerCase();
            if (itLow.indexOf('size=') === 0) continue;
            kept.push(it);
        }
        kept.push('size=' + sizeStr);

        return base + ';' + kept.join(';');
    }

    // ---- build + render ----
    async function buildCurrentSystem() {
        var spec = '';
        var kind = null;

        if (active_source === 'file' && file_state.url) {
            // extension is encoded in #fragment so phases can detect .cif/.csv
            spec = file_state.url + '#' + file_state.name;
            kind = file_state.kind;

            if (kind === 'cif') {
                var tiles = parseTilesSize();
                if (tiles) spec = applyCifSizeToSpec(spec, tiles);
            }
        } else {
            active_source = 'smiles';
            spec = (tbSmiles ? trimStr(tbSmiles.value) : '') || 'O';
            kind = 'smiles';
        }

        var built = await build_system_from_input(spec);
        atoms = built[0] || [];
        bonds = (built.length > 1 ? built[1] : null);
        setTitle(built.length > 2 ? built[2] : (kind || ''));
    }

    async function renderCurrent() {
        var w = clampInt(tbW ? tbW.value : 400, 128, 4096, 400);
        var h = clampInt(tbH ? tbH.value : 400, 128, 4096, 400);
        if (tbW) tbW.value = String(w);
        if (tbH) tbH.value = String(h);

        if (cvs.width !== w) cvs.width = w;
        if (cvs.height !== h) cvs.height = h;

        if (!atoms || atoms.length === 0) {
            clearCanvas(127);
            return;
        }

        // focal_z from atoms
        var focus_norm = rngFocus ? parseFloat(rngFocus.value) : 0.5;
        if (!Number.isFinite(focus_norm)) focus_norm = 0.5;

        var zmin = (atoms[0].z != null ? atoms[0].z : 0);
        var zmax = zmin;
        for (var i = 1; i < atoms.length; i++) {
            var z = (atoms[i].z != null ? atoms[i].z : 0);
            if (z < zmin) zmin = z;
            if (z > zmax) zmax = z;
        }
        var focal_z = zmin + focus_norm * (zmax - zmin);

        var opts = {
            bonds: bonds,
            img_size: [h, w],
            angstroms_per_pixel: rngZoom ? parseFloat(rngZoom.value) : 0.1,
            blur_sigma: rngBlur ? parseFloat(rngBlur.value) : 1.0,
            background_gray: 127,
            invert: cbInvert ? !!cbInvert.checked : false,
            noise_stddev: (cbNoise && cbNoise.checked && rngNoise) ? parseFloat(rngNoise.value) : 0.0,
            contrast: rngContrast ? parseFloat(rngContrast.value) : 1.0,
            compose_mode: 'sum',
            draw_bonds_flag: cbBonds ? !!cbBonds.checked : true,
            bond_wave_width_px: rngBwidth ? parseFloat(rngBwidth.value) : 6,
            bond_wave_amplitude: rngBamp ? parseFloat(rngBamp.value) : 0.4,
            low_clip: (el('tb-clip-lo') && el('tb-clip-lo').value !== '') ? parseFloat(el('tb-clip-lo').value) : null,
            high_clip: (el('tb-clip-hi') && el('tb-clip-hi').value !== '') ? parseFloat(el('tb-clip-hi').value) : null,
            focal_z: focal_z,
            dof_strength: rngDof ? parseFloat(rngDof.value) : 0.0,
            hide_front: cbHide ? !!cbHide.checked : false,
            show_scale_bar: cbScale ? !!cbScale.checked : false,
            scale_bar_corner: 'bl',
            scale_bar_margin_px: 12,
            canvasCtx: ctx
        };

        // renderer.js draws into canvasCtx
        render_image(atoms, opts);
    }

    var rafPending = false;
    function scheduleRender() {
        if (rafPending) return;
        rafPending = true;
        requestAnimationFrame(function () {
            rafPending = false;
            renderCurrent();
        });
    }

    var rebuildPending = false;
    var rebuildQueued = false;
    async function rebuildAndRender() {
        // Coalesce rebuild requests instead of dropping them.
        // This prevents the "typed SMILES -> temporary invalid -> fallback O -> never updates" issue
        // when RDKit build is still running.
        if (rebuildPending) {
            rebuildQueued = true;
            return;
        }

        rebuildPending = true;
        try {
            await buildCurrentSystem();
        } catch (e) {
            // Keep UI responsive even on parse errors.
            // (Minimal error handling; no debug spam.)
            atoms = [];
            bonds = null;
            setTitle('ERROR');
        } finally {
            rebuildPending = false;
        }

        await renderCurrent();

        if (rebuildQueued) {
            rebuildQueued = false;
            setTimeout(function () { rebuildAndRender(); }, 0);
        }
    }

    // ---- file load ----
    async function handleLocalFile(file) {
        if (!file) return;
        var kind = fileKindFromName(file.name);
        if (kind !== 'cif' && kind !== 'csv') return;

        revokeFileUrl();

        file_state.kind = kind;
        file_state.name = file.name;
        file_state.url = URL.createObjectURL(file);
        file_state.isBlob = true;

        active_source = 'file';
        await rebuildAndRender();
    }

    // ---- events ----

    // Render-only controls
    [rngZoom, rngContrast, rngBlur, rngNoise, rngFocus, rngDof, rngBwidth, rngBamp].forEach(function (x) {
        if (x) x.oninput = scheduleRender;
    });
    [cbInvert, cbNoise, cbBonds, cbHide, cbScale].forEach(function (x) {
        if (x) x.onchange = scheduleRender;
    });

    if (tbW) tbW.onchange = scheduleRender;
    if (tbH) tbH.onchange = scheduleRender;

    if (tbTiles) {
        tbTiles.onchange = function () {
            if (active_source === 'file' && file_state.kind === 'cif') rebuildAndRender();
        };
    }

    if (tbSmiles) {
        tbSmiles.value = trimStr(smiles_text);

        // Switching from file -> SMILES must be immediate.
        // We rebuild after a short debounce to avoid running RDKit on every keystroke.
        var smilesTimer = null;
        function scheduleSmilesRebuild() {
            active_source = 'smiles';
            if (smilesTimer) clearTimeout(smilesTimer);
            smilesTimer = setTimeout(function () {
                smilesTimer = null;
                rebuildAndRender();
            }, 250);
        }

        tbSmiles.oninput = scheduleSmilesRebuild;
        tbSmiles.onchange = function () {
            active_source = 'smiles';
            rebuildAndRender();
        };

        // Make plain Enter apply.
        // - For <textarea>: Enter applies; Shift+Enter inserts newline.
        // - For <input>: Enter applies.
        tbSmiles.addEventListener('keydown', function (ev) {
            if (ev.key !== 'Enter') return;

            var isTextarea = String(tbSmiles.tagName || '').toUpperCase() === 'TEXTAREA';
            if (isTextarea && ev.shiftKey) return; // allow newline

            ev.preventDefault();
            active_source = 'smiles';
            rebuildAndRender();
        });
    }

    if (fileInput) {    
        fileInput.onchange = function () {
            var files = fileInput.files;
            if (files && files.length) handleLocalFile(files[0]);
            fileInput.value = '';
        };
    }

    // Global drag&drop (drop anywhere)
    document.addEventListener('dragover', function (ev) { ev.preventDefault(); });
    document.addEventListener('drop', function (ev) {
        ev.preventDefault();
        var dt = ev.dataTransfer;
        if (!dt || !dt.files || !dt.files.length) return;
        handleLocalFile(dt.files[0]);
    });

    if (btnExport) {
        btnExport.onclick = function () {
            var a = document.createElement('a');
            a.download = 'tem.png';
            a.href = cvs.toDataURL('image/png');
            a.click();
        };
    }

    // ---- init ----
    active_source = 'smiles';
    await rebuildAndRender();
}
