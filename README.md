# SciGenEM / SciGenTEM Visualizer

<p align="center">
  <strong>Browser-based visualizer for atomic structures, molecules and microscopy-style educational illustrations.</strong>
</p>

<p align="center">
  <a href="https://yariktel54.github.io/SciGenEM/">Live demo</a> ·
  <a href="#quick-start">Quick start</a> ·
  <a href="#supported-inputs">Supported inputs</a> ·
  <a href="#microscopy-modes">Microscopy modes</a>
</p>

<p align="center">
  <img alt="JavaScript" src="https://img.shields.io/badge/JavaScript-browser%20first-f7df1e?style=for-the-badge&logo=javascript&logoColor=111" />
  <img alt="HTML5" src="https://img.shields.io/badge/HTML5-canvas-e34f26?style=for-the-badge&logo=html5&logoColor=fff" />
  <img alt="RDKit" src="https://img.shields.io/badge/RDKit-WASM-3366cc?style=for-the-badge" />
  <img alt="No backend required" src="https://img.shields.io/badge/no%20backend-required-4caf50?style=for-the-badge" />
</p>

---

## What is this?

**SciGenEM** is a lightweight web tool for turning molecular and crystal structure data into clean, microscopy-style illustrations.  
It is designed for educational materials, presentations, laboratory explanations and quick visual experiments where you need an atomic scene without opening a heavy desktop package.

The project runs directly in the browser: the page loads the structure, builds an atomic scene and renders the result on an HTML canvas. For SMILES processing it uses a browser chemistry backend based on **RDKit.js / OpenChemLib**, while file-based structures can be loaded directly into the scene.

> The goal is not to replace real TEM, STM or AFM analysis software.  
> The goal is to create understandable, controllable and visually convincing scientific illustrations.

---

## Why it is useful

- **Fast** — paste a SMILES string or drop a structure file and render immediately.
- **Browser-only** — no installation, no Python environment, no server-side processing.
- **Educational** — the result is suitable for slides, manuals, reports and teaching materials.
- **Flexible** — one interface for molecules, crystals, proteins and custom atomic scenes.
- **Export-ready** — save images, export coordinates and record simple GIF animations.

---

## Quick start

Open the live version:

```text
https://yariktel54.github.io/SciGenEM/
```

Or run locally:

```bash
git clone https://github.com/yariktel54/SciGenEM.git
cd SciGenEM
python -m http.server 8000
```

Then open:

```text
http://localhost:8000/
```

Do not open `index.html` directly from the file system. Browser modules, WASM files and sample loading work more reliably through a local web server.

---

## Supported inputs

SciGenEM can work with several common structure formats:

- **SMILES** — molecules from text formulas, for example `CCO`, `c1ccccc1`, `C`, `O`.
- **CIF** — crystal structures with optional tiling.
- **POSCAR / CONTCAR** — VASP-style crystal structures.
- **PDB / ENT** — protein and biomolecular structures.
- **XYZ** — simple atom-coordinate files.
- **MOL** — molecular connection table files.
- **CSV / JSON** — custom atomic scene formats.

You can use the input field, built-in samples or drag-and-drop file loading.

---

## Microscopy modes

SciGenEM currently includes three visual modes:

### TEM

Projection-style image formation with contrast, blur, noise, focus and depth-of-field controls. Useful for showing how atomic columns or molecular structures may appear in transmission-style visualization.

### STM

Surface-oriented mode with bright atomic protrusions, dark background and scanline-style effects. Useful for visualizing surface-like molecular scenes.

### AFM

Topography-inspired mode focused on surface shape, molecular outline and bond-like features. Useful for stylized molecular surface illustrations.

---

## Main features

- Interactive canvas rendering.
- Scene with multiple objects.
- Move, rotate and reset individual objects.
- CIF/POSCAR tiling for periodic structures.
- Z-layer clipping for cutting through a structure.
- Adjustable zoom, contrast, blur, noise and background.
- Per-element visual overrides.
- PNG export.
- XYZ export.
- GIF recording for simple animations.
- Ukrainian and English interface support.

---

## Project structure

```text
SciGenEM/
├── index.html              # main application page
├── help.html               # help page
├── translations.json       # UI translations
├── js/                     # application, loaders, builders and renderers
├── rdkit/                  # RDKit WASM files
├── samples/                # example structures
└── qr.png                  # QR image used by the interface
```

Important internal areas:

```text
js/app/       UI state, scene logic, input pipeline
js/builders/  structure builders for different formats
js/io/        file and SMILES input/output helpers
js/chem/      RDKit/OpenChemLib adapter layer
js/render/    TEM/STM/AFM rendering code
js/system/    shared system and provider contracts
```

---

## Example SMILES

Try these in the input field:

```text
O
C
CCO
c1ccccc1
Cn1cnc2n(C)c(=O)n(C)c(=O)c12
```

Examples for files are available in the `samples/` folder.

---

## Intended use

SciGenEM is especially useful for:

- illustrations for textbooks and lecture notes;
- visual material for chemistry, physics and materials science classes;
- quick previews of molecular or crystal structures;
- presentation figures where a real microscopy-style image is needed but full simulation is unnecessary;
- experiments with simplified TEM/STM/AFM-inspired rendering models.

---

## Notes on scientific accuracy

The renderer is built for **illustrative and educational visualization**. It uses physically inspired ideas such as atomic number dependence, blur, noise, focus and surface/topography-like rendering, but it should not be treated as a replacement for quantitative microscopy simulation software.

For teaching and explanation, this is a strength: the parameters are visible, adjustable and easy to connect with real physical concepts.

---

## Authors

Created by the SciGenEM / SciGenTEM project team.

Interface credits shown in the application include:

- Nataliia Yilmaz
- Yaroslav Teliashenko
- Pavlo Kozub

---

## Contributing

Ideas, issues and improvements are welcome. Useful directions include:

- better SMILES geometry handling;
- more realistic TEM/STM/AFM renderers;
- new import formats;
- improved educational presets;
- documentation and examples.

---

<p align="center">
  <strong>SciGenEM makes atomic structures visible, editable and useful for teaching.</strong>
</p>
