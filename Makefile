# Core entropy of the real quadratic family -- build everything from source.
#
#   make            data, both pages, the figures, the note, the previews
#   make data       core_entropy_data.csv (the sweep, via Julia, about 2 s)
#   make html       core_entropy.html, the Julia front end
#   make tent       tent_map.html, the Python page (entropy, bifurcation, tent map)
#   make figures    figs/fig_*.pdf for the note
#   make note       core_entropy_notes.pdf
#   make bands      bands.pdf (the companion note on bands)
#   make previews   preview*.png
#   make app        the native matplotlib window (needs a display)
#   make check      validation suite of both implementations
#   make clean      remove all of the above
#
# The toolchain: Julia (reference sweep), Python 3 with numpy and matplotlib,
# latexmk with the usual TeX packages, and pdftoppm for one preview.

JULIA   ?= julia
PYTHON  ?= python3
LATEXMK ?= latexmk

DATA     := core_entropy_data.csv
HTML     := core_entropy.html
TENT     := tent_map.html
NOTE     := core_entropy_notes.pdf
BANDS    := bands.pdf
FIGS     := $(addprefix figs/, fig_output.pdf fig_kneading.pdf fig_laps.pdf fig_model.pdf fig_bands.pdf \
                  fig_laps_period.pdf fig_julia.pdf \
                  fig_determinant.pdf fig_lapgrowth.pdf fig_plateau.pdf)
PREVIEWS := preview.png preview_bifurcation.png preview_kneading.png preview_tent.png
JULIA_SRC := core_entropy.jl src/CoreEntropy.jl src/checks.jl src/HTMLPage.jl

.PHONY: all data html tent figures note bands app previews check clean

all: $(DATA) $(HTML) $(TENT) $(FIGS) $(NOTE) $(BANDS) $(PREVIEWS)

# one run of the Julia front end writes the CSV and the HTML
data: $(DATA)
html: $(HTML)
$(DATA) $(HTML) &: $(JULIA_SRC)
	$(JULIA) --startup-file=no core_entropy.jl

# the Python page computes the sweep itself if the CSV is missing
tent: $(TENT)
$(TENT): tent_map.py core_entropy.py $(DATA)
	$(PYTHON) tent_map.py

figures: $(FIGS)
$(FIGS) &: figs/make_figures.py $(DATA)
	$(PYTHON) figs/make_figures.py

note: $(NOTE)
$(NOTE): core_entropy_notes.tex $(FIGS)
	$(LATEXMK) -pdf -interaction=nonstopmode core_entropy_notes.tex

# the companion note on bands and renormalization shares the figures
bands: $(BANDS)
$(BANDS): bands.tex $(FIGS)
	$(LATEXMK) -pdf -interaction=nonstopmode bands.tex

previews: $(PREVIEWS)
$(PREVIEWS) &: figs/make_previews.py $(DATA) figs/fig_kneading.pdf
	$(PYTHON) figs/make_previews.py

# the native window; not part of 'all', since it opens a GUI
app: $(DATA)
	$(PYTHON) tent_map_app.py

check: $(DATA)
	$(JULIA) --startup-file=no core_entropy.jl --check
	$(PYTHON) -c "import numpy as np, core_entropy as c; print('python: h(3.7) = %.9f' % c.core_entropy_one(3.7))"

clean:
	rm -f $(DATA) $(HTML) $(TENT) $(NOTE) $(BANDS) $(FIGS) $(PREVIEWS)
	rm -rf __pycache__
	$(LATEXMK) -c core_entropy_notes.tex >/dev/null 2>&1 || true
	$(LATEXMK) -c bands.tex >/dev/null 2>&1 || true
