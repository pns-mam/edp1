# Site generation

`python3 tools/build.py` (needs pandoc >= 3) builds the static site into `docs/`:

- `README.md` gives the main page `docs/index.html`;
- each folder gets a page `docs/<folder>/index.html` made from `<folder>/<folder>.md` (or from its only `.md` file); the other `.md` files of the folder become extra pages linked from the folder page, and PDF / GIF / PNG files are copied;
- links to `.md` files are rewritten to the HTML pages; math is rendered by MathJax in the browser.

The workflow `.github/workflows/pages.yml` runs this script on every push to `main` and commits `docs/`. Do not edit `docs/` by hand. GitHub Pages is set to serve the branch `main`, folder `/docs`.

To keep a file off the site, list its path (one per line) in `tools/exclude.txt`.

Figures made by Julia code: the reusable workflow `.github/workflows/julia-figures.yml` extracts the Julia code of a corrigé (`tools/julia_blocks.py`), runs it in the TP folder (so that the GIF / PNG it writes land there), commits the files and, on `main`, starts `pages.yml` (a push made by a workflow does not start other workflows). It is called by `figures-tp1.yml`, `figures-tp2.yml` (Plots) and `figures-tp2-3d.yml` (GLMakie under a virtual display), each started by a change of its corrigé, by a pull request touching it, or by hand (Actions tab, `workflow_dispatch`, once the workflow is on `main`). The site itself only copies the existing GIF / PNG files: the corrigé can be kept off the site (`tools/exclude.txt`) while its figures are still generated, and a statement can include a figure made by the code of its corrigé.
