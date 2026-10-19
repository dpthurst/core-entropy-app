"""
    HTMLPage

The interactive Plotly page: a self-contained HTML file with the entropy curve,
a λ slider and a panel showing g_λ and its critical orbit.  Built on the
constants and the solved curve of `CoreEntropy`.
"""
module HTMLPage

using Printf
using Downloads
using ..CoreEntropy: LAMBDA_INF, LAMBDA_P3

const PLOTLY_VERSION = "2.35.2"
const PLOTLY_URL = "https://cdn.plot.ly/plotly-$(PLOTLY_VERSION).min.js"

export write_html

# ================================================================ html

function json_vec(v::AbstractVector{Float64}, fmt::String)
    io = IOBuffer()
    print(io, '[')
    f = Printf.Format(fmt)
    for i in eachindex(v)
        i > 1 && print(io, ',')
        print(io, Printf.format(f, v[i]))
    end
    print(io, ']')
    return String(take!(io))
end

"""plotly.js, from the script's directory, a cache, or the CDN (in that order)."""
function plotly_source()
    for p in (joinpath(@__DIR__, "plotly.min.js"),
              joinpath(homedir(), ".cache", "core_entropy", "plotly.min.js"))
        isfile(p) && return read(p, String)
    end
    cache = joinpath(homedir(), ".cache", "core_entropy", "plotly.min.js")
    mkpath(dirname(cache))
    try
        Downloads.download(PLOTLY_URL, cache)
        @printf("downloaded plotly.js %s to %s\n", PLOTLY_VERSION, cache)
        return read(cache, String)
    catch err
        @warn "could not obtain plotly.js; the HTML will load it from the CDN" err
        return nothing
    end
end

const PAGE = raw"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Core entropy of the real quadratic family</title>
<style>
  html, body { background:#fff; }
  body { font-family:-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
         color:#1b1b1b; margin:0; padding:22px 30px 34px; }
  h1 { font-size:20px; font-weight:600; margin:0 0 6px; letter-spacing:-0.01em; }
  p.sub { color:#555; font-size:13.5px; line-height:1.5; margin:0 0 14px; max-width:980px; }
  code, .mono { font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; }
  code { background:#f2f3f5; padding:1px 4px; border-radius:3px; font-size:12.5px; }
  #curve { width:100%; height:470px; }
  .row { display:flex; gap:26px; align-items:flex-start; margin-top:6px; }
  #map { flex:0 0 430px; height:430px; }
  .panel { flex:1 1 auto; min-width:280px; font-size:13.5px; line-height:1.55; }
  .kv { display:grid; grid-template-columns:max-content 1fr; gap:3px 16px;
        font-variant-numeric:tabular-nums; margin:8px 0 14px; }
  .kv .k { color:#666; }
  .kv .v { font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; }
  .slider-wrap { margin:10px 0 4px; }
  input[type=range] { width:100%; accent-color:#1c5f9e; }
  .btns { margin:6px 0 12px; display:flex; gap:8px; flex-wrap:wrap; }
  .chks { display:flex; flex-direction:column; gap:5px; margin:0 0 12px; }
  label.chk { display:flex; align-items:center; gap:6px; font-size:12.5px;
              color:#333; margin:0; cursor:pointer; user-select:none; }
  label.chk input { accent-color:#1c5f9e; }
  button { font:inherit; font-size:12.5px; padding:4px 10px; border:1px solid #cfd4d9;
           background:#fafbfc; border-radius:5px; cursor:pointer; color:#222; }
  button:hover { background:#f0f3f6; }
  .marks { font-size:13px; line-height:1.75; }
  .marks a { color:#1c5f9e; text-decoration:none; border-bottom:1px dotted #9db8d2; }
  .marks a:hover { border-bottom-style:solid; }
  .note { color:#666; font-size:12.5px; margin-top:12px; }
  .seq { word-break:break-all; line-height:1.4; }
</style>
</head>
<body>

<h1>Core entropy of the real quadratic family</h1>
<p class="sub">
  <b>h(λ)</b> is the topological entropy of <code>g_λ(x) = λx(1−x)</code> restricted to its core
  (the Hubbard tree), for <code>λ ∈ [2,4]</code> — equivalently of <code>z ↦ z²+c</code> with
  <code>c = λ/2 − λ²/4 ∈ [−2,0]</code>. It is continuous and nondecreasing, flat on every
  hyperbolic component, and runs from <code>h(2)=0</code> to <code>h(4)=log 2</code>
  (Tiozzo 2015). Values are computed from the Milnor–Thurston kneading determinant of the
  critical orbit. Hover or zoom the curve; click it (or drag the slider) to inspect the map.
</p>

<div id="curve"></div>

<div class="row">
  <div id="map"></div>
  <div class="panel">
    <div class="slider-wrap">
      <input type="range" id="slider" min="2" max="4" step="0.0002" value="3.8318741">
    </div>
    <div class="btns">
      <button id="zoom-full">full range [2,4]</button>
      <button id="zoom-stair">zoom to the staircase [3.55,4]</button>
    </div>
    <div class="chks">
      <label class="chk"><input type="checkbox" id="curveOn" checked> entropy curve</label>
      <label class="chk"><input type="checkbox" id="bif"> bifurcation diagram
        (x on the right axis)</label>
    </div>
    <div class="kv">
      <span class="k">λ</span><span class="v" id="v-lam"></span>
      <span class="k">c = λ/2 − λ²/4</span><span class="v" id="v-c"></span>
      <span class="k">h(λ)</span><span class="v" id="v-h"></span>
      <span class="k">h / log 2</span><span class="v" id="v-hn"></span>
      <span class="k">critical orbit</span><span class="v seq" id="v-seq"></span>
      <span class="k">eventually periodic</span><span class="v" id="v-per"></span>
    </div>
    <div class="marks">
      <b>landmarks</b><br>
      <a href="#" data-lam="2">λ = 2</a>, c = 0 — critical point fixed, h = 0<br>
      <a href="#" data-lam="3">λ = 3</a>, c = −3/4 — parabolic fixed point, h = 0<br>
      <a href="#" data-lam="3.2360679775">λ = 1+√5</a>, c = −1 — basilica centre, h = 0<br>
      <a href="#" data-lam="3.5699456718">λ∞ ≈ 3.5699457</a>, c ≈ −1.4011552 — Feigenbaum
      accumulation of period doublings, h = 0<br>
      <a href="#" data-lam="3.8318740553">λ ≈ 3.8318741</a>, c ≈ −1.7548777 — period-3 centre,
      h = log φ ≈ 0.4812118<br>
      <a href="#" data-lam="4">λ = 4</a>, c = −2 — Chebyshev map, h = log 2 ≈ 0.6931472
    </div>
    <div class="note">
      The dotted lines mark the Feigenbaum parameter and the period-3 centre. The shaded band in
      the lower panel is the hull of the critical orbit (the core).
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

const APP = raw"""
const LAM = @LAM@;
const H   = @H@;
const NPT = LAM.length;
const LOG2 = Math.log(2);
const cOf  = l => l/2 - l*l/4;
const lamOf = c => 1 + Math.sqrt(Math.max(0, 1 - 4*c));

// nearest grid parameter
function idxOf(lam){
  if (lam <= LAM[0]) return 0;
  if (lam >= LAM[NPT-1]) return NPT-1;
  let lo = 0, hi = NPT-1;
  while (hi - lo > 1){ const m = (lo + hi) >> 1; if (LAM[m] <= lam) lo = m; else hi = m; }
  return (lam - LAM[lo] <= LAM[hi] - lam) ? lo : hi;
}
const hAt = lam => H[idxOf(lam)];

const CD = new Array(NPT);
for (let i = 0; i < NPT; i++) CD[i] = [cOf(LAM[i]), H[i]/LOG2];

// ---------------------------------------------------------------- main curve
const gdCurve = document.getElementById('curve');
const cTicks = [-2, -1.75, -1.5, -1.25, -1, -0.75, -0.5, -0.25, 0];
const LAM_INF = @LAM_INF@, LAM_P3 = @LAM_P3@;

const layoutCurve = {
  margin: {l: 66, r: 26, t: 48, b: 46},
  xaxis: {title: {text: 'λ'}, range: [2, 4], zeroline: false,
          gridcolor: '#eceff2', rangeslider: {visible: true, thickness: 0.07, bgcolor: '#f5f7f9'}},
  yaxis: {title: {text: 'core entropy   h(λ)'}, range: [-0.02, 0.73], zeroline: false,
          gridcolor: '#eceff2'},
  yaxis2: {overlaying: 'y', side: 'right', range: [0, 1], visible: false,
           title: {text: 'x   (limit set)'}, showgrid: false, zeroline: false,
           automargin: true},
  xaxis2: {overlaying: 'x', side: 'top', range: [2, 4], showgrid: false, zeroline: false,
           title: {text: 'c   (z ↦ z² + c)'},
           tickvals: cTicks.map(lamOf), ticktext: cTicks.map(String)},
  hovermode: 'closest', showlegend: false,
  font: {family: '-apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif',
         size: 12, color: '#222'},
  paper_bgcolor: '#fff', plot_bgcolor: '#fff',
  shapes: [
    {type: 'line', x0: LAM_INF, x1: LAM_INF, y0: 0, y1: 1, yref: 'paper',
     line: {color: '#c9ced4', width: 1, dash: 'dot'}},
    {type: 'line', x0: LAM_P3, x1: LAM_P3, y0: 0, y1: 1, yref: 'paper',
     line: {color: '#c9ced4', width: 1, dash: 'dot'}}
  ],
  annotations: [
    {x: LAM_INF, y: 1, yref: 'paper', text: 'λ∞', showarrow: false,
     yshift: 6, xshift: -10, font: {size: 11, color: '#8a9099'}},
    {x: LAM_P3, y: 1, yref: 'paper', text: 'period 3', showarrow: false,
     yshift: 6, xshift: 16, font: {size: 11, color: '#8a9099'}}
  ]
};

// Bifurcation diagram: the attracting set of g_λ, drawn behind the curve and
// read on the right-hand axis (x ∈ [0,1]).  Instead of single dots it is an
// occupancy grid --- many iterates per parameter deposited into bins --- shown
// as a greyscale density, which is what makes the branches of the periodic
// orbits (the fold lines) appear as continuous curves.  Recomputed in the
// browser for the visible λ-range, so zooming sharpens it.
const trBif = {type: 'heatmap', x: [2], y: [0.5], z: [[0]], yaxis: 'y2',
               colorscale: [[0, 'rgba(255,255,255,0)'], [0.10, 'rgba(60,60,60,0.28)'],
                            [0.40, 'rgba(25,25,25,0.80)'], [1, 'rgba(0,0,0,0.95)']],
               zmin: 0, zmax: 1, showscale: false, hoverinfo: 'skip', visible: false};

const trCurve = {
  x: LAM, y: H, customdata: CD, type: 'scatter', mode: 'lines',
  line: {color: '#1c5f9e', width: 1.3},
  hovertemplate: 'λ = %{x:.6f}<br>c = %{customdata[0]:.6f}<br>' +
                 'h = %{y:.6f}<br>h / log 2 = %{customdata[1]:.4f}<extra></extra>'
};
const trMark = {x: [3.8318741], y: [hAt(3.8318741)], type: 'scatter', mode: 'markers',
                marker: {size: 9, color: '#c1272d', line: {color: '#fff', width: 1.2}},
                hoverinfo: 'skip'};
// Plotly only creates a secondary axis if some trace lives on it, so keep an
// invisible one on x2: it is what makes the top c-axis exist.
const trX2 = {x: [2, 4], y: [0, 0], xaxis: 'x2', yaxis: 'y', type: 'scatter', mode: 'markers',
              marker: {size: 2, opacity: 0, color: '#ffffff'},
              showlegend: false, hoverinfo: 'skip'};
Plotly.newPlot(gdCurve, [trBif, trCurve, trMark, trX2], layoutCurve,
               {responsive: true, displaylogo: false,
                modeBarButtonsToRemove: ['lasso2d', 'select2d']});

// ---------------------------------------------------------------- map panel
const gdMap = document.getElementById('map');
const layoutMap = {
  margin: {l: 48, r: 14, t: 30, b: 40},
  xaxis: {title: {text: 'x'}, range: [0, 1], zeroline: false, gridcolor: '#eef1f3'},
  yaxis: {title: {text: 'g_λ(x)'}, range: [0, 1.02], zeroline: false, gridcolor: '#eef1f3',
          scaleanchor: 'x', scaleratio: 1},
  showlegend: false,
  font: {family: '-apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif',
         size: 12, color: '#222'},
  paper_bgcolor: '#fff', plot_bgcolor: '#fff'
};

function mapTraces(lam){
  const M = 240, xs = [], ys = [];
  for (let i = 0; i <= M; i++){ const x = i/M; xs.push(x); ys.push(lam*x*(1-x)); }
  let x = 0.5; const ox = [], oy = [];
  for (let k = 0; k < 160; k++){ const y = lam*x*(1-x); ox.push(x); oy.push(y); x = y; }
  return [
    {x: [0, 1], y: [0, 1], type: 'scatter', mode: 'lines',
     line: {color: '#b9c0c7', width: 1, dash: 'dot'}, hoverinfo: 'skip'},
    {x: xs, y: ys, type: 'scatter', mode: 'lines', line: {color: '#1c5f9e', width: 2},
     hovertemplate: 'x = %{x:.4f}<br>g_λ(x) = %{y:.4f}<extra></extra>'},
    {x: ox, y: oy, type: 'scatter', mode: 'markers',
     marker: {size: 5, color: '#e08a1e', opacity: 0.75}, hoverinfo: 'skip'}
  ];
}

function coreHull(lam){
  let z = 0.5, mn = 1, mx = 0;
  for (let k = 0; k < 6000; k++){
    z = lam*z*(1-z);
    if (k > 2500){ if (z < mn) mn = z; if (z > mx) mx = z; }
  }
  return [mn, mx];
}

function drawMap(lam){
  const hull = coreHull(lam);
  const mn = hull[0], mx = hull[1];
  const lay = Object.assign({}, layoutMap, {
    title: {text: 'g_λ on [0,1]', font: {size: 13}, x: 0.5, xanchor: 'center'},
    shapes: [
      {type: 'rect', x0: mn, x1: mx, y0: 0, y1: 1, yref: 'paper',
       fillcolor: 'rgba(224,138,30,0.09)', line: {width: 0}, layer: 'below'},
      {type: 'line', x0: 0.5, x1: 0.5, y0: 0, y1: 1, yref: 'paper',
       line: {color: '#c9ced4', width: 1, dash: 'dot'}}
    ]
  });
  Plotly.react(gdMap, mapTraces(lam), lay, {responsive: true, displaylogo: false,
      modeBarButtonsToRemove: ['lasso2d', 'select2d']});
}

// ---------------------------------------------------------------- readout
function itinerary(lam, maxn){
  let x = 0.5; const xs = [x], sym = [];
  for (let k = 0; k < maxn; k++){
    const y = lam*x*(1-x);
    sym.push(y < 0.5 ? 'L' : (y > 0.5 ? 'R' : '*'));
    xs.push(y); x = y;
  }
  let per = null, pre = null;
  for (let n = 1; n < xs.length && per === null; n++){
    for (let m = 0; m < n; m++){
      if (Math.abs(xs[n] - xs[m]) < 1e-12){ per = n - m; pre = m; break; }
    }
  }
  return {sym: sym, per: per, pre: pre};
}

let curLam = 3.8318741;
const fmt = (x, d) => x.toFixed(d);

function updateReadout(lam, h){
  document.getElementById('v-lam').textContent = fmt(lam, 6);
  document.getElementById('v-c').textContent = fmt(cOf(lam), 6);
  document.getElementById('v-h').textContent = fmt(h, 6);
  document.getElementById('v-hn').textContent = fmt(h/LOG2, 4);
  const it = itinerary(lam, 120);
  document.getElementById('v-seq').textContent =
      it.sym.slice(0, 42).join('') + (it.sym.length > 42 ? '…' : '');
  document.getElementById('v-per').textContent =
      (it.per === null) ? 'no (within 120 iterates)'
                        : 'period ' + it.per + (it.pre > 0 ? ' (after ' + it.pre + ')' : '');
}

function setLambda(lam, fromSlider){
  lam = Math.max(2, Math.min(4, lam));
  curLam = lam;
  const h = hAt(lam);
  Plotly.restyle(gdCurve, {x: [[lam]], y: [[h]]}, [2]);
  drawMap(lam);
  updateReadout(lam, h);
  if (!fromSlider) document.getElementById('slider').value = lam;
}

// ---------------------------------------------------------------- wiring
document.getElementById('slider').addEventListener('input', e => {
  setLambda(parseFloat(e.target.value), true);
});
gdCurve.on('plotly_click', e => {
  if (e.points && e.points.length) setLambda(e.points[0].x, false);
});
document.querySelectorAll('a[data-lam]').forEach(a => {
  a.addEventListener('click', ev => {
    ev.preventDefault();
    setLambda(parseFloat(a.getAttribute('data-lam')), false);
  });
});
document.getElementById('zoom-full').addEventListener('click', () => {
  Plotly.relayout(gdCurve, {'xaxis.range': [2, 4], 'xaxis2.range': [2, 4]});
});
document.getElementById('zoom-stair').addEventListener('click', () => {
  Plotly.relayout(gdCurve, {'xaxis.range': [3.55, 4], 'xaxis2.range': [3.55, 4]});
});

// ---------------------------------------------------------------- bifurcation overlay
const BIF_COLS = 1000, BIF_ROWS = 400, BIF_TRANSIENT = 1500, BIF_KEEP = 1200;

// occupancy grid: z[row][col], row = bin of x, col = bin of λ
function bifurcationGrid(x0, x1){
  const z = new Array(BIF_ROWS);
  for (let i = 0; i < BIF_ROWS; i++) z[i] = new Float64Array(BIF_COLS);
  for (let j = 0; j < BIF_COLS; j++){
    const lam = BIF_COLS < 2 ? x0 : x0 + (x1 - x0)*j/(BIF_COLS - 1);
    let x = 0.5;
    for (let k = 0; k < BIF_TRANSIENT; k++) x = lam*x*(1 - x);
    for (let k = 0; k < BIF_KEEP; k++){
      x = lam*x*(1 - x);
      let r = Math.floor(x*BIF_ROWS);
      if (r < 0) r = 0; else if (r >= BIF_ROWS) r = BIF_ROWS - 1;
      z[r][j] += 1;
    }
  }
  return z;
}

// log-scaled so that a densely visited point (a periodic orbit) is dark and a
// thinly visited one (a chaotic band) is pale
function normaliseGrid(z){
  let mx = 0;
  for (let i = 0; i < BIF_ROWS; i++){
    const row = z[i];
    for (let j = 0; j < BIF_COLS; j++) if (row[j] > mx) mx = row[j];
  }
  const s = mx > 0 ? 1/Math.log1p(mx) : 0;
  const out = new Array(BIF_ROWS);
  for (let i = 0; i < BIF_ROWS; i++){
    const row = z[i], o = new Array(BIF_COLS);
    for (let j = 0; j < BIF_COLS; j++) o[j] = row[j] > 0 ? Math.log1p(row[j])*s : 0;
    out[i] = o;
  }
  return out;
}

function columnCenters(x0, x1){
  const xs = new Array(BIF_COLS);
  for (let j = 0; j < BIF_COLS; j++)
    xs[j] = BIF_COLS < 2 ? x0 : x0 + (x1 - x0)*j/(BIF_COLS - 1);
  return xs;
}

function rowCenters(){
  const ys = new Array(BIF_ROWS);
  for (let i = 0; i < BIF_ROWS; i++) ys[i] = (i + 0.5)/BIF_ROWS;
  return ys;
}

const cbBif = document.getElementById('bif');
let bifTimer = null;

function visibleRange(){
  const r = gdCurve.layout.xaxis.range;
  return (r && isFinite(r[0]) && isFinite(r[1]) && r[1] > r[0]) ? r : [2, 4];
}

function drawBifurcation(){
  const r = visibleRange();
  const z = normaliseGrid(bifurcationGrid(r[0], r[1]));
  Plotly.restyle(gdCurve, {x: [columnCenters(r[0], r[1])], y: [rowCenters()], z: [z]}, [0]);
}

cbBif.addEventListener('change', e => {
  const on = e.target.checked;
  if (on) drawBifurcation();
  Plotly.restyle(gdCurve, {visible: [on]}, [0]);
  Plotly.relayout(gdCurve, {'yaxis2.visible': on});
});

gdCurve.on('plotly_relayout', () => {
  if (!cbBif.checked) return;
  clearTimeout(bifTimer);
  bifTimer = setTimeout(drawBifurcation, 120);
});

// hiding the entropy curve leaves the bifurcation diagram on its own axis
const cbCurve = document.getElementById('curveOn');
cbCurve.addEventListener('change', e => {
  const on = e.target.checked;
  Plotly.restyle(gdCurve, {visible: [on, on]}, [1, 2]);   // curve and its marker
  Plotly.relayout(gdCurve, {'yaxis.visible': on});
});

setLambda(curLam, false);
"""

function write_html(path::String, lams::Vector{Float64}, hs::Vector{Float64};
                    cdn::Bool = false)
    lamjson = json_vec(lams, "%.9g")
    hjson = json_vec(hs, "%.8g")
    app = replace(APP, "@LAM@" => lamjson, "@H@" => hjson,
                  "@LAM_INF@" => string(LAMBDA_INF),
                  "@LAM_P3@" => string(LAMBDA_P3))
    js = cdn ? nothing : plotly_source()
    script = if js === nothing
        "<script src=\"$(PLOTLY_URL)\"></script>"
    else
        "<script>$(js)</script>"
    end
    html = replace(PAGE, "@PLOTLY@" => script, "@APP@" => app)
    open(path, "w") do io
        write(io, html)
    end
    return path
end


end # module HTMLPage
