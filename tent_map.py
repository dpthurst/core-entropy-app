#!/usr/bin/env python3
"""
tent_map.py -- g_lambda and its semi-conjugate tent map, side by side.

Writes tent_map.html, a self-contained interactive page.  For a logistic
parameter lambda the tent map is taken with slope

    s = exp(h(lambda)) = 1/t_0 ,

the exponential of the core entropy, i.e. the reciprocal of the smallest zero
of the kneading determinant.  With that slope the tent map has the same
kneading sequence as g_lambda, so the two are semi-conjugate: there is a
continuous nondecreasing h: [0,1] -> [0,1] with h o f_s = g_lambda o h.

The page shows both maps with the cobweb of the critical orbit, the common
itinerary, and the semi-conjugacy h itself, which is drawn by matching the
turning points of f_s^n with those of g_lambda^n.

    python3 tent_map.py        # writes tent_map.html
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PLOTLY_VERSION = "2.35.2"
PLOTLY_URL = f"https://cdn.plot.ly/plotly-{PLOTLY_VERSION}.min.js"


def plotly_js():
    for p in (os.path.join(HERE, "plotly.min.js"),
              os.path.expanduser("~/.cache/core_entropy/plotly.min.js")):
        if os.path.isfile(p):
            return open(p, encoding="utf-8").read()
    raise SystemExit("plotly.min.js not found; run core_entropy.jl once, or download it")


def json_vec(v, fmt="%.9g"):
    return "[" + ",".join(fmt % x for x in v) + "]"


PAGE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>g_λ and its semi-conjugate tent map</title>
<style>
  body { font-family:-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
         color:#1b1b1b; margin:0; padding:20px 26px 30px; background:#fff; }
  h1 { font-size:19px; font-weight:600; margin:0 0 4px; }
  p.sub { color:#555; font-size:13px; line-height:1.5; margin:0 0 12px; max-width:1000px; }
  code, .mono { font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; }
  code { background:#f2f3f5; padding:1px 4px; border-radius:3px; font-size:12.5px; }
  .row { display:flex; gap:18px; align-items:flex-start; }
  .row > div { flex:1 1 0; min-width:0; }
  #logi, #tentp { height:340px; }
  #semi { height:330px; flex:2 1 0; }
  .info { flex:1 1 0; font-size:13px; line-height:1.6; padding-top:6px; }
  .kv { display:grid; grid-template-columns:max-content 1fr; gap:2px 14px;
        font-variant-numeric:tabular-nums; margin-bottom:10px; }
  .kv .k { color:#666; }
  .kv .v { font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; }
  .seq { word-break:break-all; line-height:1.45; font-size:12.5px; }
  input[type=range] { width:100%; accent-color:#1c5f9e; margin:14px 0 2px; }
  .ticks { display:flex; justify-content:space-between; color:#888; font-size:11px; }
</style>
</head>
<body>

<h1>g<sub>λ</sub> and its semi-conjugate tent map</h1>
<p class="sub">
  For <code>g_λ(x) = λx(1−x)</code> take the tent map <code>f_s</code> of slope
  <code>s = e^{h(λ)} = 1/t₀</code>, the exponential of the core entropy: it always has the
  same entropy as <code>g_λ</code>.  Where the two also share a kneading sequence they are
  semi-conjugate, and a continuous nondecreasing <code>h</code> with
  <code>h ∘ f_s = g_λ ∘ h</code> carries one orbit to the other.  That is the case at 86% of
  the parameters here; elsewhere <code>g_λ</code> sits inside a hyperbolic component, whose
  attracting cycle no tent map can have, and then no tent map shares its sequence.  Below:
  the two maps with the cobweb of the critical orbit, the itineraries, and the map
  <code>h</code> itself, drawn by matching the turning points of the iterates.
</p>

<div class="row">
  <div id="logi"></div>
  <div id="tentp"></div>
</div>
<div class="row" style="margin-top:6px">
  <div id="semi"></div>
  <div class="info">
    <div class="kv">
      <span class="k">λ</span><span class="v" id="v-lam"></span>
      <span class="k">h(λ)</span><span class="v" id="v-h"></span>
      <span class="k">s = e^{h(λ)} = 1/t₀</span><span class="v" id="v-s"></span>
      <span class="k">turning points</span><span class="v" id="v-n"></span>
    </div>
    <div><b>itinerary</b> <span class="v" id="v-cmp"></span></div>
    <div class="seq mono" id="v-seqL"></div>
    <div class="seq mono" id="v-seqT"></div>
    <div style="color:#666;font-size:12px;margin-top:8px">
      first line: g<sub>λ</sub>; second: f<sub>s</sub>.  L and R are the two sides of the
      critical point, * a return to it.
    </div>
  </div>
</div>

<input type="range" id="slider" min="3.6" max="4" step="0.0002" value="3.7">
<div class="ticks"><span>3.6</span><span>λ</span><span>4</span></div>

@PLOTLY@
<script>
@APP@
</script>
</body>
</html>
"""

APP = r"""
const LAM = @LAM@;
const S   = @S@;
const NPT = LAM.length;
const TOL = 1e-9;

const g    = (x, lam) => lam * x * (1 - x);
const tent = (y, s)   => y < 0.5 ? s * y : s * (1 - y);

function sOf(lam){
  if (lam <= LAM[0]) return S[0];
  if (lam >= LAM[NPT-1]) return S[NPT-1];
  let lo = 0, hi = NPT-1;
  while (hi - lo > 1){ const m = (lo + hi) >> 1; if (LAM[m] <= lam) lo = m; else hi = m; }
  const t = (lam - LAM[lo]) / (LAM[hi] - LAM[lo]);
  return S[lo] + t * (S[hi] - S[lo]);
}

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

function dedupe(v, tol){
  v.sort((a, b) => a - b);
  const out = [];
  for (const z of v){
    if (!out.length || z - out[out.length-1] > tol) out.push(z);
  }
  return out;
}

// turning points of the two maps: union of preimages of 1/2 up to level n.
// Both trees are pruned by the same kneading sequence, so the two sorted sets
// have the same length and correspond point by point under the semi-conjugacy.
function turningPoints(lam, s, maxLevel, target){
  let SL = [0.5], ST = [0.5];
  let AL = [0.5], AT = [0.5];
  let level = 0, ok = true;
  for (; level < maxLevel; level++){
    const nl = [], nt = [];
    for (const y of SL){
      const disc = 1 - 4 * y / lam;
      if (disc > 0){ const r = Math.sqrt(disc); nl.push(0.5*(1+r), 0.5*(1-r)); }
    }
    for (const y of ST){
      if (y <= s/2 + 1e-12){
        const a = y / s, b = 1 - y / s;
        if (a > 0 && a < 1) nt.push(a);
        if (b > 0 && b < 1) nt.push(b);
      }
    }
    const dl = dedupe(nl, TOL), dt = dedupe(nt, TOL);
    if (dl.length !== dt.length){ ok = false; break; }   // numerical mismatch
    SL = dl; ST = dt;
    const pl = AL, pt = AT;
    AL = dedupe(AL.concat(SL), TOL);
    AT = dedupe(AT.concat(ST), TOL);
    if (AL.length !== AT.length){ AL = pl; AT = pt; ok = false; break; }
    if (AL.length >= target || SL.length === 0) break;
  }
  return {logistic: AL, tent: AT, level: level, ok: ok};
}

// ---------------------------------------------------------------- plots
const colL = '#1c5f9e', colT = '#2f6b4f', colC = '#c0392b', grey = '#b9c0c7';
const font = {family: '-apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif',
              size: 11, color: '#222'};
const cfg = {responsive: true, displaylogo: false,
             modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d']};

function ax(title){
  return {title: {text: title, font: {size: 11}}, range: [-0.02, 1.04], zeroline: false,
          gridcolor: '#eef1f3', tickfont: {size: 10}};
}
function baseLayout(xt, yt){
  return {margin: {l: 44, r: 12, t: 30, b: 34}, xaxis: ax(xt), yaxis: ax(yt),
          showlegend: false, font: font, paper_bgcolor: '#fff', plot_bgcolor: '#fff'};
}

const gdL = document.getElementById('logi');
const gdT = document.getElementById('tentp');
const gdH = document.getElementById('semi');

const trL = [
  {x: [0,1], y: [0,1], mode: 'lines', line: {color: grey, width: 1, dash: 'dot'}, hoverinfo: 'skip'},
  {x: [], y: [], mode: 'lines', line: {color: colL, width: 1.8}, hoverinfo: 'skip'},
  {x: [], y: [], mode: 'lines', line: {color: colC, width: 0.7, simplify: false}, hoverinfo: 'skip'},
  {x: [], y: [], mode: 'markers', marker: {size: 4, color: colC}, hoverinfo: 'skip'}
];
const trT = [
  {x: [0,1], y: [0,1], mode: 'lines', line: {color: grey, width: 1, dash: 'dot'}, hoverinfo: 'skip'},
  {x: [], y: [], mode: 'lines', line: {color: colT, width: 1.8}, hoverinfo: 'skip'},
  {x: [], y: [], mode: 'lines', line: {color: colC, width: 0.7, simplify: false}, hoverinfo: 'skip'},
  {x: [], y: [], mode: 'markers', marker: {size: 4, color: colC}, hoverinfo: 'skip'}
];
const trH = [
  {x: [0,1], y: [0,1], mode: 'lines', line: {color: grey, width: 1, dash: 'dot'}, hoverinfo: 'skip'},
  {x: [], y: [], mode: 'lines', line: {color: '#444', width: 1.2}, hoverinfo: 'skip'},
  {x: [], y: [], mode: 'markers', marker: {size: 6, color: colC, line: {color: '#fff', width: 1}},
   hovertemplate: 'f_s^k(½) = %{x:.4f}<br>g_λ^k(½) = %{y:.4f}<extra></extra>'}
];

Plotly.newPlot(gdL, trL, baseLayout('x', 'g_λ(x)'), cfg);
Plotly.newPlot(gdT, trT, baseLayout('y', 'f_s(y)'), cfg);
Plotly.newPlot(gdH, trH, baseLayout('y   (tent coordinate)', 'x = h(y)'), cfg);

// ---------------------------------------------------------------- update
const NFILL = 34;          // cobweb steps
const NSYM  = 15;          // itinerary symbols compared
function setLambda(lam){
  const s = sOf(lam);
  const h = Math.log(s);
  const orbL = orbit(x => g(x, lam), NFILL, 0.5);
  const orbT = orbit(y => tent(y, s), NFILL, 0.5);
  const cwL = cobweb(x => g(x, lam), NFILL, 0.5);
  const cwT = cobweb(y => tent(y, s), NFILL, 0.5);
  const xs = [], ys = [];
  for (let i = 0; i <= 220; i++){ const x = i/220; xs.push(x); ys.push(g(x, lam)); }

  Plotly.restyle(gdL, {x: [xs, cwL.x, orbL.slice(0, 7)], y: [ys, cwL.y, orbL.slice(0, 7)]}, [1, 2, 3]);
  Plotly.restyle(gdT, {x: [[0, 0.5, 1], cwT.x, orbT.slice(0, 7)],
                       y: [[0, s/2, 0], cwT.y, orbT.slice(0, 7)]}, [1, 2, 3]);

  const tp = turningPoints(lam, s, 40, 3000);
  const ok = tp.ok && tp.logistic.length >= 30;
  const hx = [0].concat(tp.tent, [1]);
  const hy = [0].concat(tp.logistic, [1]);
  Plotly.restyle(gdH, {x: [hx, orbT.slice(0, 7)], y: [hy, orbL.slice(0, 7)],
                       visible: [ok, ok]}, [1, 2]);
  Plotly.relayout(gdH, {annotations: ok ? [] : [{
      x: 0.5, y: 0.55, xref: 'paper', yref: 'paper', showarrow: false,
      text: 'no tent map shares this kneading sequence',
      font: {size: 11, color: '#999'}}]});

  const sL = itinerary(x => g(x, lam), NSYM, 1e-6);
  const sT = itinerary(y => tent(y, s), NSYM, 1e-6);
  document.getElementById('v-lam').textContent = lam.toFixed(6);
  document.getElementById('v-h').textContent = h.toFixed(6);
  document.getElementById('v-s').textContent = s.toFixed(6);
  document.getElementById('v-n').textContent = tp.logistic.length + '  (level ' + tp.level + ')';
  document.getElementById('v-seqL').textContent = 'g_λ : ' + sL;
  document.getElementById('v-seqT').textContent = 'f_s : ' + sT;
  const same = sL === sT;
  document.getElementById('v-cmp').textContent = same
      ? 'identical on the first ' + NSYM + ' symbols'
      : 'differ: hyperbolic component, no tent partner';
  document.getElementById('v-cmp').style.color = same ? '#2f6b4f' : '#c0392b';
}

document.getElementById('slider').addEventListener('input', e => setLambda(parseFloat(e.target.value)));
setLambda(parseFloat(document.getElementById('slider').value));
"""


def main():
    d = np.loadtxt(os.path.join(HERE, "core_entropy_data.csv"), delimiter=",", skiprows=1)
    lam, h = d[:, 0], d[:, 2]
    s = np.exp(h)
    html = (PAGE.replace("@PLOTLY@", "<script>" + plotly_js() + "</script>")
                .replace("@APP@", APP.replace("@LAM@", json_vec(lam)).replace("@S@", json_vec(s))))
    out = os.path.join(HERE, "tent_map.html")
    with open(out, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"wrote {out} ({os.path.getsize(out)/2**20:.1f} MB)")


if __name__ == "__main__":
    main()
