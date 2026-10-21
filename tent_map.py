#!/usr/bin/env python3
"""
tent_map.py -- core entropy, bifurcation diagram and the semi-conjugate tent
map of the logistic family, on one interactive page.

Writes tent_map.html.  Everything is computed in Python:

  * core_entropy.py sweeps h(λ) from the Milnor-Thurston kneading determinant
    (the same algorithm as src/CoreEntropy.jl, and it agrees with it to 5e-11);
  * the bifurcation diagram is an occupancy grid of the attracting set,
    computed in the browser for the visible λ-range;
  * for each λ the tent map has slope s = exp(h(λ)) = 1/t₀, and Theorem 7.4 of
    Milnor-Thurston gives a continuous nondecreasing A with f_s(A(x)) = A(g_λ(x)).

    python3 tent_map.py            # self-contained page, plotly inlined (5 MB)
    python3 tent_map.py --cdn      # same page, plotly from cdn.plot.ly (0.6 MB)
"""
import os
import numpy as np
import core_entropy as ce

HERE = os.path.dirname(os.path.abspath(__file__))
PLOTLY_VERSION = "2.35.2"
CACHE = os.path.join(os.path.expanduser("~"), ".cache", "core_entropy", "plotly.min.js")


def plotly_script(cdn=False):
    """The plotly tag: inline (self-contained) or from the CDN (0.6 MB pages)."""
    if cdn:
        return ('<script src="https://cdn.plot.ly/plotly-%s.min.js" charset="utf-8">'
                '</script>' % PLOTLY_VERSION)
    return "<script>" + plotly_js() + "</script>"


def plotly_js():
    for p in (os.path.join(HERE, "plotly.min.js"), CACHE):
        if os.path.isfile(p):
            return open(p, encoding="utf-8").read()
    raise SystemExit("plotly.min.js not found; run core_entropy.jl once, or download it")


def curve():
    """The sweep, reused from core_entropy_data.csv when it matches the grid."""
    path = os.path.join(HERE, "core_entropy_data.csv")
    lams = ce.lambda_grid()
    if os.path.isfile(path):
        d = np.loadtxt(path, delimiter=",", skiprows=1)
        if d.shape[0] == len(lams) and np.allclose(d[:, 0], lams, atol=1e-9):
            print(f"reusing the sweep in {os.path.basename(path)}")
            return d[:, 0], d[:, 2]
    print("computing the sweep in Python (about 40 s) ...")
    hs = ce.compute_curve(lams)
    ce.write_csv(path, lams, hs)
    return lams, hs


def json_vec(v, fmt="%.9g"):
    return "[" + ",".join(fmt % x for x in v) + "]"


PAGE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Core entropy, bifurcation diagram and the tent map</title>
<style>
  body { font-family:-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
         color:#1b1b1b; margin:0; padding:18px 24px 26px; background:#fff; }
  h1 { font-size:19px; font-weight:600; margin:0 0 4px; }
  p.sub { color:#555; font-size:13px; line-height:1.5; margin:0 0 10px; max-width:1100px; }
  code, .mono { font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; }
  code { background:#f2f3f5; padding:1px 4px; border-radius:3px; font-size:12.5px; }
  #entropy { width:100%; height:300px; }
  .row { display:flex; gap:16px; align-items:flex-start; margin-top:4px; }
  .row > div { flex:1 1 0; min-width:0; }
  #logi, #tentp, #semi { height:300px; }
  .ctrl { display:flex; gap:18px; align-items:center; margin:6px 0 2px; }
  .ctrl input[type=range] { flex:1 1 auto; accent-color:#1c5f9e; }
  .ticks { display:flex; justify-content:space-between; color:#888; font-size:11px; margin:0 0 8px; }
  .info { display:flex; gap:30px; font-size:13px; line-height:1.55; }
  .kv { display:grid; grid-template-columns:max-content 1fr; gap:2px 14px;
        font-variant-numeric:tabular-nums; }
  .kv .k { color:#666; }
  .kv .v { font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; }
  .seq { word-break:break-all; line-height:1.45; font-size:12.5px; }
  .chk { font-size:12.5px; color:#333; white-space:nowrap; }
  button.btn { font:inherit; font-size:12.5px; padding:1px 8px; border-radius:4px;
               border:1px solid #c8ccd2; background:#f6f7f9; color:#222; cursor:pointer; }
  button.btn:hover { background:#eceef1; }
  button.step { padding:1px 6px; min-width:22px; }
  /* the entropy plot's margins are l:56, r:46, so the slider lines up with its axes */
  #lamctrl, #lamticks { padding-left:56px; padding-right:46px; }
  #lamctrl { margin-top:2px; }
  /* a range input insets its track by half a thumb at each end, so the box has to be
     pulled out by that much for the track to start and end where the axis does
     (Chromium's thumb is about 16 px) */
  #lamctrl input[type=range] { margin-left:-8px; margin-right:-8px; }
  .help { background:#f6f8fa; border-left:3px solid #1c5f9e; border-radius:0 4px 4px 0;
          padding:8px 14px 9px; margin:0 0 12px; max-width:1100px;
          font-size:12.5px; line-height:1.55; color:#333; }
  .help .hd { font-weight:600; color:#111; margin-bottom:2px; }
  .help ul { margin:0; padding-left:18px; }
  .help li { margin:1px 0; }
</style>
</head>
<body>

<h1>Core entropy, bifurcation diagram, and the semi-conjugate tent map</h1>
<p class="sub">
  Top: the core entropy <code>h(λ)</code> of <code>g_λ(x) = λx(1−x)</code>, computed from the
  Milnor–Thurston kneading determinant, with the bifurcation diagram of the same parameters
  behind it as a greyscale density (read on the right-hand axis).  The dashed line is the
  Feigenbaum parameter, the red one the current <code>λ</code>; zoom and pan act on
  <code>λ</code> alone.  Below: the cobweb of the critical orbit for <code>g_λ</code> and for
  the tent map <code>f_s</code> of slope <code>s = e^{h(λ)} = 1/t₀</code> — or, with
  <i>plot the iterates</i>, the graphs of <code>g_λ^n</code> and <code>f_s^n</code> for
  <code>n</code> up to 20, titled with their lap counts — and the map <code>A</code> of
  Theorem 7.4 of Milnor–Thurston, which satisfies <code>f_s(A(x)) = A(g_λ(x))</code>: the
  itinerary of <code>x</code>, read through <code>f_s</code>'s inverse branches, lands on the
  itinerary of <code>g_λ(x)</code> read the same way.  The arrow does not reverse, and
  <code>A</code> need not be injective.
</p>

<div class="help">
  <div class="hd">Using it</div>
  <ul>
    <li>Drag the <code>λ</code> slider, just under the plot, to move through the family.  Its two
        ends are the ends of the visible range, so it always spans what you are looking at.</li>
    <li>Zoom by dragging a box on the <code>λ</code> axis.  Only <code>λ</code> zooms: the entropy
        axis keeps its scale, and so do the three map panels.</li>
    <li>To come back, double-click the plot, press the home button in its mode bar, or use
        <i>reset zoom</i> (below the maps), or press <b>0</b>.</li>
    <li><b>+</b> and <b>&minus;</b>, or the buttons beside the slider, or the arrow keys on the
        slider itself, change the number of iterates, from 1 to 20.</li>
    <li>The checkboxes below the maps switch the bifurcation diagram, the iterates view (in which
        case the cobweb-length slider appears in its place) and the <code>A</code> panel.</li>
  </ul>
</div>

<div id="entropy"></div>

<div class="ctrl" id="lamctrl">
  <input type="range" id="slider" min="2" max="4" step="0.001" value="3.7">
</div>
<div class="ticks" id="lamticks"><span id="tick0">2</span><span>λ</span><span id="tick1">4</span></div>

<div class="row">
  <div id="logi"></div>
  <div id="tentp"></div>
  <div id="semi"></div>
</div>

<div class="ctrl">
  <label class="chk"><input type="checkbox" id="bif" checked> bifurcation diagram</label>
  <label class="chk"><input type="checkbox" id="iter" checked> plot the iterates
    <span class="mono">g<sup>n</sup>, f<sub>s</sub><sup>n</sup></span></label>
  <span class="chk" id="nwrap">n =
    <button type="button" class="btn step" id="nminus" title="one fewer iterate (the - key)">&minus;</button>
    <input type="range" id="nslider" min="1" max="20" step="1" value="1" style="width:90px;vertical-align:middle">
    <button type="button" class="btn step" id="nplus" title="one more iterate (the + key)">+</button>
    <span class="mono" id="nlabel">1</span></span>
  <span class="chk" id="fwrap" style="display:none">cobweb steps =
    <input type="range" id="fillslider" min="2" max="1000" step="1" value="34" style="width:90px;vertical-align:middle">
    <span class="mono" id="filllabel">34</span></span>
  <label class="chk"><input type="checkbox" id="ashow" checked> A panel</label>
  <button type="button" class="btn" id="resetzoom" title="back to the whole family (the 0 key)">reset zoom</button>
</div>
<div class="info">
  <div class="kv">
    <span class="k">λ</span><span class="v" id="v-lam"></span>
    <span class="k">h(λ)</span><span class="v" id="v-h"></span>
    <span class="k">c = λ/2 − λ²/4</span><span class="v" id="v-c"></span>
    <span class="k">s = e^{h(λ)}</span><span class="v" id="v-s"></span>
    <span class="k">h(λ)/log 2</span><span class="v" id="v-hr"></span>
    <span class="k">max |f_s(A) − A(g_λ)|</span><span class="v" id="v-err"></span>
  </div>
  <div>
    <div><b>itineraries of the critical orbits</b></div>
    <div class="seq mono" id="v-seqL"></div>
    <div class="seq mono" id="v-seqT"></div>
    <div style="color:#666;font-size:12px;margin-top:4px">
      L and R are the two sides of the critical point, * a return to it.  They need not
      agree: where <code>A</code> collapses an interval, the two itineraries can differ.
    </div>
  </div>
</div>

@PLOTLY@
<script>
@APP@
</script>
</body>
</html>
"""

APP = r"""
const LAM = @LAM@;
const H   = @H@;
const NPT = LAM.length;
const LOG2 = Math.log(2);

const g    = (x, lam) => lam * x * (1 - x);
const tent = (y, s)   => y < 0.5 ? s * y : s * (1 - y);

function idxOf(lam){
  if (lam <= LAM[0]) return 0;
  if (lam >= LAM[NPT-1]) return NPT-1;
  let lo = 0, hi = NPT-1;
  while (hi - lo > 1){ const m = (lo + hi) >> 1; if (LAM[m] <= lam) lo = m; else hi = m; }
  return (lam - LAM[lo] <= LAM[hi] - lam) ? lo : hi;
}
const hAt = lam => H[idxOf(lam)];
const sAt = lam => Math.exp(hAt(lam));

function orbit(f, n, x0){
  const xs = [x0];
  for (let k = 0; k < n; k++) xs.push(f(xs[xs.length-1]));
  return xs;
}
function cobweb(f, n, x0){
  const xs = [], ys = [];
  let x = x0;
  for (let k = 0; k < n; k++){
    const y = f(x);
    xs.push(x, x, y); ys.push(x, y, y);
    x = y;
  }
  return {x: xs, y: ys};
}
function itinerary(f, n, tol){
  let x = 0.5; const out = [];
  for (let k = 0; k < n; k++){
    x = f(x);
    out.push(Math.abs(x - 0.5) < tol ? '*' : (x < 0.5 ? 'L' : 'R'));
  }
  return out.join('');
}
function A_of(x, lam, s, n){
  let y = x; const w = [];
  for (let k = 0; k < n; k++){ w.push(y < 0.5 ? 1 : -1); y = lam*y*(1-y); }
  let u = 0.5;
  for (let k = w.length - 1; k >= 0; k--) u = w[k] > 0 ? u/s : 1 - u/s;
  return u;
}

// ------------------------------------------------ the bifurcation diagram
const BIF_COLS = 1000, BIF_ROWS = 400, BIF_TRANSIENT = 1500, BIF_KEEP = 1200;
const BIF_TOL = 1e-14;   // below the smallest gap between distinct preimages of 1/2
const LAP_CAP = 120000;  // stop counting turning points above this: 2^n at lambda = 4

function dedupe(v, tol){
  v.sort((a, b) => a - b);
  const out = [];
  for (const z of v) if (!out.length || z - out[out.length-1] > tol) out.push(z);
  return out;
}
function bifurcationGrid(x0, x1){
  const z = new Array(BIF_ROWS);
  for (let i = 0; i < BIF_ROWS; i++) z[i] = new Float64Array(BIF_COLS);
  for (let j = 0; j < BIF_COLS; j++){
    const lam = BIF_COLS < 2 ? x0 : x0 + (x1 - x0)*j/(BIF_COLS - 1);
    let x = 0.5;
    for (let k = 0; k < BIF_TRANSIENT; k++) x = g(x, lam);
    for (let k = 0; k < BIF_KEEP; k++){
      x = g(x, lam);
      let r = Math.floor(x*BIF_ROWS);
      if (r < 0) r = 0; else if (r >= BIF_ROWS) r = BIF_ROWS - 1;
      z[r][j] += 1;
    }
  }
  return z;
}
function normaliseGrid(z){
  let mx = 0;
  for (let i = 0; i < BIF_ROWS; i++){
    const row = z[i];
    for (let j = 0; j < BIF_COLS; j++) if (row[j] > mx) mx = row[j];
  }
  const sc = mx > 0 ? 1/Math.log1p(mx) : 0;
  const out = new Array(BIF_ROWS);
  for (let i = 0; i < BIF_ROWS; i++){
    const row = z[i], o = new Array(BIF_COLS);
    for (let j = 0; j < BIF_COLS; j++) o[j] = row[j] > 0 ? 0.6 + 0.4*Math.log1p(row[j])*sc : 0;
    out[i] = o;
  }
  return out;
}
function colCenters(x0, x1){
  const xs = new Array(BIF_COLS);
  for (let j = 0; j < BIF_COLS; j++) xs[j] = BIF_COLS < 2 ? x0 : x0 + (x1 - x0)*j/(BIF_COLS - 1);
  return xs;
}
function rowCenters(){
  const ys = new Array(BIF_ROWS);
  for (let i = 0; i < BIF_ROWS; i++) ys[i] = (i + 0.5)/BIF_ROWS;
  return ys;
}

// number of laps of the n-th iterate: 1 + the number of turning points, i.e.
// the distinct interior preimages of 1/2 under g^k, k = 0..n-1
function lapsG(lam, n){
  let S = [0.5], all = [0.5];
  for (let k = 0; k < n-1; k++){
    const nx = [];
    for (const y of S){
      const disc = 1 - 4*y/lam;
      if (disc > 0){ const r = Math.sqrt(disc); nx.push(0.5*(1+r), 0.5*(1-r)); }
    }
    S = dedupe(nx, BIF_TOL);
    all = dedupe(all.concat(S), BIF_TOL);
    if (all.length > LAP_CAP) return null;
  }
  return all.filter(p => p > 0 && p < 1).length + 1;
}
function lapsTent(s, n){
  let S = [0.5], all = [0.5];
  for (let k = 0; k < n-1; k++){
    const nx = [];
    for (const y of S){
      if (y <= s/2 + 1e-12){
        const a = y/s, b = 1 - y/s;
        if (a > 0 && a < 1) nx.push(a);
        if (b > 0 && b < 1) nx.push(b);
      }
    }
    S = dedupe(nx, BIF_TOL);
    all = dedupe(all.concat(S), BIF_TOL);
    if (all.length > LAP_CAP) return null;
  }
  return all.filter(p => p > 0 && p < 1).length + 1;
}

// ------------------------------------------------ plots
const font = {family: '-apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif',
              size: 11, color: '#222'};
const cfg = {responsive: true, displaylogo: false, doubleClick: 'reset',
             modeBarButtonsToRemove: ['lasso2d', 'select2d']};
const grey = '#b9c0c7';

const gdE = document.getElementById('entropy');
const gdL = document.getElementById('logi');
const gdT = document.getElementById('tentp');
const gdH = document.getElementById('semi');

const layoutE = {
  margin: {l: 56, r: 46, t: 24, b: 34},
  xaxis: {title: {text: 'λ', font: {size: 11}}, range: [2, 4], gridcolor: '#eef1f3',
          tickfont: {size: 10}},
  yaxis: {title: {text: 'core entropy  h(λ)', font: {size: 11}}, range: [-0.02, 0.70],
          gridcolor: '#eef1f3', tickfont: {size: 10}, fixedrange: true},
  yaxis2: {overlaying: 'y', side: 'right', range: [0, 1], showgrid: false, fixedrange: true,
           title: {text: 'x   (limit set)', font: {size: 10}}, tickfont: {size: 9},
           automargin: true},
  showlegend: false, font: font, paper_bgcolor: '#fff', plot_bgcolor: '#fff',
  hovermode: 'closest', dragmode: 'zoom',
  shapes: [
    {type: 'line', x0: 3.5699456718, x1: 3.5699456718, yref: 'paper', y0: 0, y1: 1,
     line: {color: '#777', width: 0.8, dash: 'dash'}},
    {type: 'line', x0: 3.7, x1: 3.7, yref: 'paper', y0: 0, y1: 1,
     line: {color: '#c1272d', width: 1.0}}
  ]
};
const trE = [
  {type: 'heatmap', x: [2], y: [0.5], z: [[0]], yaxis: 'y2',
   colorscale: [[0, 'rgba(255,255,255,0)'], [1, 'rgba(90,101,112,0.95)']],
   zmin: 0, zmax: 1, showscale: false, hoverinfo: 'skip'},
  {x: LAM, y: H, type: 'scatter', mode: 'lines', line: {color: '#1c5f9e', width: 1.2},
   hovertemplate: 'λ = %{x:.6f}<br>h = %{y:.6f}<extra></extra>'},
  {x: [3.7], y: [0], type: 'scatter', mode: 'markers',
   marker: {size: 10, color: '#c1272d', line: {color: '#fff', width: 1.2}}, hoverinfo: 'skip'}
];

function panel(title, xt, yt, colour){
  return {margin: {l: 44, r: 12, t: 28, b: 32},
          xaxis: {title: {text: xt, font: {size: 11}}, range: [-0.02, 1.04],
                  gridcolor: '#eef1f3', tickfont: {size: 10}},
          yaxis: {title: {text: yt, font: {size: 11}}, range: [-0.02, 1.04],
                  gridcolor: '#eef1f3', tickfont: {size: 10}},
          title: {text: title, font: {size: 11}}, showlegend: false, font: font,
          paper_bgcolor: '#fff', plot_bgcolor: '#fff'};
}
const trL = [
  {x: [0,1], y: [0,1], mode: 'lines', line: {color: grey, width: 1, dash: 'dot'}, hoverinfo: 'skip'},
  {x: [], y: [], mode: 'lines', line: {color: '#1c5f9e', width: 1.6}, hoverinfo: 'skip'},
  {x: [], y: [], mode: 'lines', line: {color: '#c0392b', width: 0.35}, hoverinfo: 'skip'},
  {x: [], y: [], mode: 'markers', marker: {size: 4, color: '#c0392b'}, hoverinfo: 'skip'}
];
const trT = [
  {x: [0,1], y: [0,1], mode: 'lines', line: {color: grey, width: 1, dash: 'dot'}, hoverinfo: 'skip'},
  {x: [], y: [], mode: 'lines', line: {color: '#2f6b4f', width: 1.6}, hoverinfo: 'skip'},
  {x: [], y: [], mode: 'lines', line: {color: '#c0392b', width: 0.35}, hoverinfo: 'skip'},
  {x: [], y: [], mode: 'markers', marker: {size: 4, color: '#c0392b'}, hoverinfo: 'skip'}
];
const trH = [
  {x: [0,1], y: [0,1], mode: 'lines', line: {color: grey, width: 1, dash: 'dot'}, hoverinfo: 'skip'},
  {x: [], y: [], mode: 'lines', line: {color: '#333', width: 1.2}, hoverinfo: 'skip'},
  {x: [], y: [], mode: 'markers', marker: {size: 6, color: '#c0392b', line: {color: '#fff', width: 1}},
   hoverinfo: 'skip'}
];

Plotly.newPlot(gdE, trE, layoutE, cfg);
Plotly.newPlot(gdL, trL, panel('g_λ on its interval', 'x', 'g_λ(x)'), cfg);
Plotly.newPlot(gdT, trT, panel('the tent map f_s', 'y', 'f_s(y)'), cfg);
Plotly.newPlot(gdH, trH, panel('Milnor–Thurston A:  f_s ∘ A = A ∘ g_λ',
                               'x  (interval of g_λ)', 'A(x)  (tent coordinate)'), cfg);

// ------------------------------------------------ update everything
const NSYM = 48, NG = 240;   // symbols shown for each critical orbit
const fillSlider = document.getElementById('fillslider');
const fillLabel = document.getElementById('filllabel');
const fillSteps = () => parseInt(fillSlider.value, 10);
const cbBif = document.getElementById('bif');
let bifTimer = null;

const lamSlider = document.getElementById('slider');
const tick0 = document.getElementById('tick0');
const tick1 = document.getElementById('tick1');
let syncing = false;

function syncSliderToView(){
  if (syncing) return;
  const r = visibleRange();
  const lo = Math.max(2, r[0]), hi = Math.min(4, r[1]);
  if (!(hi - lo > 1e-6)) return;
  syncing = true;
  lamSlider.min = lo.toFixed(5);
  lamSlider.max = hi.toFixed(5);
  lamSlider.step = ((hi - lo)/2000).toFixed(8);
  tick0.textContent = lo.toFixed(4);
  tick1.textContent = hi.toFixed(4);
  let v = parseFloat(lamSlider.value);
  if (v < lo || v > hi){
    v = Math.min(hi, Math.max(lo, v));
    lamSlider.value = v;
    syncing = false;
    setLambda(v);            // keep lambda inside the visible range
    return;
  }
  syncing = false;
}

function visibleRange(){
  const r = gdE.layout.xaxis.range;
  return (r && isFinite(r[0]) && isFinite(r[1]) && r[1] > r[0]) ? r : [2, 4];
}
function drawBifurcation(){
  const r = visibleRange();
  const z = normaliseGrid(bifurcationGrid(r[0], r[1]));
  Plotly.restyle(gdE, {x: [colCenters(r[0], r[1])], y: [rowCenters()], z: [z]}, [0]);
}

const iterBox = document.getElementById('iter');
const nSlider = document.getElementById('nslider');
const nLabel  = document.getElementById('nlabel');
const nWrap   = document.getElementById('nwrap');

// either the cobweb of the critical orbit, or the graph of the n-th iterate
// itself, which is what shows the lap count exploding
function drawMaps(lam, s){
  const orbL = orbit(x => g(x, lam), fillSteps(), 0.5);
  if (iterBox.checked){
    const n = Math.max(1, Math.min(20, parseInt(nSlider.value, 10)));
    const lG = lapsG(lam, n), lT = lapsTent(s, n);
    // sample by the number of laps: a fixed grid aliases a curve with thousands of them
    const M = Math.min(40001, Math.max(2401, 30*(lG || LAP_CAP))), xs = [], yl = [], yt = [];
    for (let i = 0; i <= M; i++){
      const x = i/M;
      let a = x, b = x;
      for (let k = 0; k < n; k++){ a = g(a, lam); b = tent(b, s); }
      xs.push(x); yl.push(a); yt.push(b);
    }
    const lw = Math.max(0.5, 1.6 - 0.1*n);
    Plotly.restyle(gdL, {line: {width: lw}}, [1]);
    Plotly.restyle(gdT, {line: {width: lw}}, [1]);
    Plotly.restyle(gdL, {x: [xs, [], []], y: [yl, [], []]}, [1, 2, 3]);
    Plotly.restyle(gdT, {x: [xs, [], []], y: [yt, [], []]}, [1, 2, 3]);
    const lgTxt = lG === null ? '&gt;' + LAP_CAP : lG;
    const ltTxt = lT === null ? '&gt;' + LAP_CAP : lT;
    Plotly.relayout(gdL, {title: {text: 'g_λ^' + n + '   with ℓ_' + n + ' = ' + lgTxt},
                          'xaxis.title.text': 'x', 'yaxis.title.text': 'g_λ^' + n + '(x)'});
    Plotly.relayout(gdT, {title: {text: 'f_s^' + n + '   with ℓ_' + n + ' = ' + ltTxt},
                          'xaxis.title.text': 'y', 'yaxis.title.text': 'f_s^' + n + '(y)'});
    nLabel.textContent = n;
  } else {
    Plotly.restyle(gdL, {line: {width: 1.6}}, [1]);
    Plotly.restyle(gdT, {line: {width: 1.6}}, [1]);
    const cwL = cobweb(x => g(x, lam), fillSteps(), 0.5);
    const cwT = cobweb(y => tent(y, s), fillSteps(), 0.5);
    const orbT = orbit(y => tent(y, s), fillSteps(), 0.5);
    const xs = [], ys = [];
    for (let i = 0; i <= 220; i++){ const x = i/220; xs.push(x); ys.push(g(x, lam)); }
    Plotly.restyle(gdL, {x: [xs, cwL.x, orbL.slice(0, 7)], y: [ys, cwL.y, orbL.slice(0, 7)]}, [1, 2, 3]);
    Plotly.restyle(gdT, {x: [[0, 0.5, 1], cwT.x, orbT.slice(0, 7)],
                         y: [[0, s/2, 0], cwT.y, orbT.slice(0, 7)]}, [1, 2, 3]);
    Plotly.relayout(gdL, {title: {text: 'g_λ on its interval'},
                          'xaxis.title.text': 'x', 'yaxis.title.text': 'g_λ(x)'});
    Plotly.relayout(gdT, {title: {text: 'the tent map f_s'},
                          'xaxis.title.text': 'y', 'yaxis.title.text': 'f_s(y)'});
  }
  return orbL;
}

let curLam = 3.7;

function setLambda(lam){
  lam = Math.max(2, Math.min(4, lam));   // the whole family; below lambda_inf, h = 0
                                          // and the model is the tent map of slope 1
  const h = hAt(lam), s = Math.exp(h);
  curLam = lam;
  Plotly.restyle(gdE, {x: [[lam]], y: [[h]]}, [2]);

  const orbL = drawMaps(lam, s);

  const gx = [], gy = [];
  let err = 0;
  for (let i = 0; i <= NG; i++){
    const x = i/NG, a = A_of(x, lam, s, 64);
    gx.push(x); gy.push(a);
    const e = Math.abs(A_of(g(x, lam), lam, s, 64) - tent(a, s));
    if (e > err) err = e;
  }
  const oy = orbL.slice(0, 7).map(x => A_of(x, lam, s, 64));
  Plotly.restyle(gdH, {x: [gx, orbL.slice(0, 7)], y: [gy, oy]}, [1, 2]);
  // below the Feigenbaum parameter h = 0: the model has slope 1 and A collapses to a
  // point (A = 1/2 is a semi-conjugacy, with residual 0, but a degenerate one)
  Plotly.relayout(gdH, {title: {text: h > 1e-12
    ? 'Milnor–Thurston A:  f_s ∘ A = A ∘ g_λ'
    : 'Milnor–Thurston A:  h = 0 here, so f_s has slope 1 and A is constant'}});

  Plotly.relayout(gdE, {shapes: [
    {type: 'line', x0: 3.5699456718, x1: 3.5699456718, yref: 'paper', y0: 0, y1: 1,
     line: {color: '#777', width: 0.8, dash: 'dash'}},
    {type: 'line', x0: lam, x1: lam, yref: 'paper', y0: 0, y1: 1,
     line: {color: '#c1272d', width: 1.0}}
  ]});
  document.getElementById('v-lam').textContent = lam.toFixed(6);
  document.getElementById('v-c').textContent = (lam/2 - lam*lam/4).toFixed(6);
  document.getElementById('v-hr').textContent = (h/LOG2).toFixed(4);
  document.getElementById('v-h').textContent = h.toFixed(6);
  document.getElementById('v-s').textContent = s.toFixed(6);
  document.getElementById('v-err').textContent = err.toExponential(1);
  document.getElementById('v-seqL').textContent = 'g_λ : ' + itinerary(x => g(x, lam), NSYM, 1e-6);
  document.getElementById('v-seqT').textContent = 'f_s : ' + itinerary(y => tent(y, s), NSYM, 1e-6);
}

document.getElementById('slider').addEventListener('input', e => setLambda(parseFloat(e.target.value)));
iterBox.addEventListener('change', e => {
  nWrap.style.display = e.target.checked ? '' : 'none';
  document.getElementById('fwrap').style.display = e.target.checked ? 'none' : '';
  drawMaps(curLam, Math.exp(hAt(curLam)));
});
function stepN(d){
  nSlider.value = Math.max(1, Math.min(20, parseInt(nSlider.value, 10) + d));
  nLabel.textContent = nSlider.value;
  if (iterBox.checked) drawMaps(curLam, Math.exp(hAt(curLam)));
}
nSlider.addEventListener('input', () => {
  nLabel.textContent = nSlider.value;
  if (iterBox.checked) drawMaps(curLam, Math.exp(hAt(curLam)));
});
document.getElementById('nminus').addEventListener('click', () => stepN(-1));
document.getElementById('nplus').addEventListener('click', () => stepN(+1));

// the + and - keys do what the buttons do.  Sliders and checkboxes keep the keys:
// after dragging the slider it still has the focus, and swallowing them there would
// make the shortcuts look broken.
document.addEventListener('keydown', e => {
  const el = e.target || {};
  const tag = el.tagName || '';
  const type = (el.type || '').toLowerCase();
  if (tag === 'TEXTAREA') return;
  if (tag === 'INPUT' && type !== 'range' && type !== 'checkbox') return;
  if (e.key === '+' || e.key === '=') { stepN(+1); e.preventDefault(); }
  else if (e.key === '-' || e.key === '_') { stepN(-1); e.preventDefault(); }
  else if (e.key === '0') { resetZoom(); e.preventDefault(); }
});
fillSlider.addEventListener('input', () => {
  fillLabel.textContent = fillSlider.value;
  if (!iterBox.checked) drawMaps(curLam, Math.exp(hAt(curLam)));
});
document.getElementById('ashow').addEventListener('change', e => {
  document.getElementById('semi').style.display = e.target.checked ? '' : 'none';
  Plotly.Plots.resize(gdL); Plotly.Plots.resize(gdT);
});
function resetZoom(){
  Plotly.relayout(gdE, {'xaxis.autorange': true});
  setTimeout(syncSliderToView, 0);
  [gdL, gdT, gdH].forEach(gd =>
    Plotly.relayout(gd, {'xaxis.autorange': true, 'yaxis.autorange': true}));
}
document.getElementById('resetzoom').addEventListener('click', () => {
  resetZoom();
});
cbBif.addEventListener('change', e => {
  if (e.target.checked) drawBifurcation();
  Plotly.restyle(gdE, {visible: [e.target.checked]}, [0]);
});
gdE.on('plotly_relayout', () => {
  syncSliderToView();
  if (!cbBif.checked) return;
  clearTimeout(bifTimer);
  bifTimer = setTimeout(drawBifurcation, 140);
});

drawBifurcation();
setLambda(parseFloat(document.getElementById('slider').value));
syncSliderToView();
drawMaps(curLam, Math.exp(hAt(curLam)));
"""


def main():
    import sys
    lams, hs = curve()
    html = (PAGE.replace("@PLOTLY@", plotly_script("--cdn" in sys.argv))
                .replace("@APP@", APP.replace("@LAM@", json_vec(lams)).replace("@H@", json_vec(hs))))
    out = os.path.join(HERE, "tent_map.html")
    with open(out, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"wrote {out} ({os.path.getsize(out)/2**20:.1f} MB)")
    print(f"h(3.7) = {ce.h_at(lams, hs, 3.7):.6f}, s = {np.exp(ce.h_at(lams, hs, 3.7)):.6f}")


if __name__ == "__main__":
    main()
