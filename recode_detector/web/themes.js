/* Ten looks. A look = CSS variables + a layout class + a background animation.
   Pick one in the header ("Look") or open the page with ?look=<id>. The choice is remembered. */
const TH = {};                                   // colours the canvas plots read; refreshed by applyLook()
const LOOKS = [
  { id: 'ultraviolet', name: 'Ultraviolet', bg: 'lcms', layout: 'l-float',
    vars: {},
    pal: ['#ff5ca8', '#8b7bff', '#ffb86b', '#5ec8ff', '#c6f26b', '#f87171', '#e879f9', '#94a3b8'],
    heat: [[27, 19, 64], [139, 123, 255], [255, 92, 168]], div: [[27, 19, 64], [255, 122, 89], [94, 168, 255]], ns: 'rgba(154,147,196,.4)' },

  { id: 'solar', name: 'Solar Flare', bg: 'spectrum', layout: 'l-flat l-uc',
    vars: {
      '--bg1': '#130905', '--bg2': '#2b1307', '--ink': '#fff1e2', '--mute': '#c9a58a', '--lilac': '#ffd0a8',
      '--a': '#ff9a3c', '--b': '#ff4d2e', '--c': '#ffd166', '--d': '#ffb4a2', '--e': '#e9ff70',
      '--a-rgb': '255,154,60', '--b-rgb': '255,77,46', '--c-rgb': '255,209,102', '--d-rgb': '255,180,162', '--line-rgb': '255,170,110',
      '--up': '#ff5a36', '--down': '#4aa3ff', '--card1': 'rgba(60,28,10,.78)', '--card2': 'rgba(24,11,4,.82)',
      '--header-bg': 'rgba(28,13,5,.8)', '--input-bg': 'rgba(18,8,3,.7)', '--thead-bg': '#2a1307', '--code-bg': 'rgba(255,154,60,.16)',
      '--code-ink': '#ffcf9e', '--tip-bg': 'rgba(24,11,4,.96)', '--detail-bg': 'rgba(14,6,2,.75)',
      '--veil': 'radial-gradient(ellipse at 30% 0%,rgba(19,9,5,.15),rgba(19,9,5,.6) 68%)', '--btn-ink': '#1c0a02',
      '--h2-grad': 'linear-gradient(95deg,#fff3e2 10%,#ffd08a 50%,#ff7a45)', '--num-grad': 'linear-gradient(95deg,#fff3e2,#ffd08a)',
      '--body-ink': '#f0d6bf', '--toggle-off': '#4a2410',
      '--font-sans': '"Futura","Century Gothic","Trebuchet MS",sans-serif', '--font-head': '"Futura","Century Gothic","Trebuchet MS",sans-serif',
      '--head-weight': '700', '--head-case': 'uppercase', '--head-track': '.01em', '--radius': '8px', '--radius-sm': '6px', '--pill': '8px',
      '--grad': 'linear-gradient(135deg,#ff4d2e,#ff9a3c)' },
    pal: ['#ff9a3c', '#ff4d2e', '#ffd166', '#ffb4a2', '#e9ff70', '#f87171', '#fb923c', '#a8a29e'],
    heat: [[43, 17, 8], [255, 154, 60], [255, 235, 160]], div: [[43, 17, 8], [255, 90, 54], [74, 163, 255]], ns: 'rgba(201,165,138,.4)' },

  { id: 'paper', name: 'Paper Lab', bg: 'ridges', layout: 'l-flat', light: true,
    vars: {
      '--bg1': '#f7f1e5', '--bg2': '#ece0c8', '--ink': '#1d2a3a', '--mute': '#6d6657', '--lilac': '#1f3a93',
      '--a': '#1f3a93', '--b': '#d9452b', '--c': '#c98a1b', '--d': '#2f7d6d', '--e': '#7a9a2b',
      '--a-rgb': '31,58,147', '--b-rgb': '217,69,43', '--c-rgb': '201,138,27', '--d-rgb': '47,125,109', '--line-rgb': '60,50,30',
      '--up': '#d9452b', '--down': '#1f5fbf', '--bad': '#b42318', '--card1': 'rgba(255,253,247,.88)', '--card2': 'rgba(252,246,232,.88)',
      '--header-bg': 'rgba(247,241,229,.92)', '--input-bg': 'rgba(255,255,255,.85)', '--thead-bg': '#efe5cf', '--code-bg': 'rgba(217,69,43,.1)',
      '--code-ink': '#a63219', '--tip-bg': 'rgba(255,252,244,.98)', '--detail-bg': 'rgba(245,236,216,.92)',
      '--veil': 'linear-gradient(180deg,rgba(247,241,229,.1),rgba(247,241,229,.25))', '--btn-ink': '#fff', '--btn-alt-ink': '#3a2400',
      '--btn-alt-bg': 'linear-gradient(135deg,#f0c36a,#d9962b)', '--h2-grad': 'linear-gradient(#1d2a3a,#1d2a3a)',
      '--num-grad': 'linear-gradient(95deg,#1d2a3a,#1f3a93)', '--body-ink': '#3a3f4a', '--toggle-off': '#cdbf9f',
      '--shadow': '0 14px 30px -20px rgba(60,40,10,.45),0 0 0 1px rgba(60,50,30,.04)',
      '--font-head': '"Iowan Old Style","Palatino Linotype",Georgia,serif', '--head-weight': '600', '--head-track': '-.02em',
      '--radius': '6px', '--radius-sm': '4px', '--pill': '6px', '--grad': 'linear-gradient(135deg,#d9452b,#c23a22)' },
    pal: ['#1f3a93', '#d9452b', '#c98a1b', '#2f7d6d', '#7a9a2b', '#b42318', '#7c3aed', '#6d6657'],
    heat: [[250, 244, 230], [31, 58, 147], [217, 69, 43]], div: [[250, 244, 230], [217, 69, 43], [31, 95, 191]], ns: 'rgba(110,100,80,.35)' },

  { id: 'synthwave', name: 'Neon Grid', bg: 'grid', layout: 'l-float l-uc',
    vars: {
      '--bg1': '#0a0614', '--bg2': '#240a3d', '--ink': '#fdf4ff', '--mute': '#b69ad6', '--lilac': '#e8b4ff',
      '--a': '#00e5ff', '--b': '#ff2bd6', '--c': '#ffe14d', '--d': '#7b61ff', '--e': '#39ff88',
      '--a-rgb': '0,229,255', '--b-rgb': '255,43,214', '--c-rgb': '255,225,77', '--d-rgb': '123,97,255', '--line-rgb': '0,229,255',
      '--up': '#ff6a3d', '--down': '#22d3ee', '--card1': 'rgba(24,8,48,.8)', '--card2': 'rgba(10,4,24,.84)', '--header-bg': 'rgba(14,5,32,.84)',
      '--input-bg': 'rgba(8,3,20,.8)', '--thead-bg': '#1b0a38', '--code-bg': 'rgba(255,43,214,.16)', '--code-ink': '#ff9df0',
      '--tip-bg': 'rgba(10,4,24,.97)', '--detail-bg': 'rgba(6,2,16,.8)',
      '--veil': 'radial-gradient(ellipse at 50% 0%,rgba(10,6,20,.05),rgba(10,6,20,.45) 75%)', '--btn-ink': '#0a0614',
      '--h2-grad': 'linear-gradient(95deg,#fff 5%,#7dfcff 40%,#ff7be9 90%)', '--num-grad': 'linear-gradient(95deg,#fff,#7dfcff)',
      '--body-ink': '#e4d3f5', '--toggle-off': '#2b1055',
      '--shadow': '0 0 0 1px rgba(0,229,255,.12),0 0 30px -8px rgba(0,229,255,.35)',
      '--font-head': '"Avenir Next Condensed","Futura","Trebuchet MS",sans-serif', '--head-weight': '800', '--head-case': 'uppercase', '--head-track': '.03em',
      '--radius': '4px', '--radius-sm': '4px', '--pill': '4px', '--grad': 'linear-gradient(90deg,#00e5ff,#ff2bd6)' },
    pal: ['#00e5ff', '#ff2bd6', '#ffe14d', '#7b61ff', '#39ff88', '#ff6a3d', '#c084fc', '#94a3b8'],
    heat: [[14, 8, 30], [123, 97, 255], [255, 43, 214]], div: [[14, 8, 30], [255, 106, 61], [34, 211, 238]], ns: 'rgba(182,154,214,.4)' },

  { id: 'forest', name: 'Mycelium', bg: 'network', layout: 'l-float',
    vars: {
      '--bg1': '#06130d', '--bg2': '#0e2a1c', '--ink': '#eaf7ee', '--mute': '#8fb59f', '--lilac': '#b8e6c8',
      '--a': '#3ddc84', '--b': '#f2b705', '--c': '#ff8a5b', '--d': '#7bd0ff', '--e': '#c6f26b',
      '--a-rgb': '61,220,132', '--b-rgb': '242,183,5', '--c-rgb': '255,138,91', '--d-rgb': '123,208,255', '--line-rgb': '120,220,160',
      '--up': '#ff8a5b', '--down': '#7bd0ff', '--card1': 'rgba(14,48,32,.78)', '--card2': 'rgba(6,22,14,.82)', '--header-bg': 'rgba(7,24,16,.8)',
      '--input-bg': 'rgba(4,16,10,.7)', '--thead-bg': '#0f2c1d', '--code-bg': 'rgba(61,220,132,.14)', '--code-ink': '#a5f0c4',
      '--tip-bg': 'rgba(5,18,11,.96)', '--detail-bg': 'rgba(3,12,8,.75)',
      '--veil': 'radial-gradient(ellipse at 28% 0%,rgba(6,19,13,.2),rgba(6,19,13,.6) 66%)', '--btn-ink': '#04150d',
      '--h2-grad': 'linear-gradient(95deg,#fff 10%,#c9f5d8 50%,#f6d56a)', '--num-grad': 'linear-gradient(95deg,#fff,#c9f5d8)',
      '--body-ink': '#cfe8d8', '--toggle-off': '#1d4a33',
      '--font-sans': '"Gill Sans","Optima","Trebuchet MS",sans-serif', '--font-head': '"Gill Sans","Optima","Trebuchet MS",sans-serif',
      '--head-weight': '700', '--head-track': '-.01em', '--radius': '14px', '--radius-sm': '10px', '--grad': 'linear-gradient(135deg,#f2b705,#3ddc84)' },
    pal: ['#3ddc84', '#f2b705', '#ff8a5b', '#7bd0ff', '#c6f26b', '#f87171', '#e879f9', '#94a3b8'],
    heat: [[6, 30, 18], [61, 220, 132], [242, 183, 5]], div: [[8, 32, 20], [255, 138, 91], [123, 208, 255]], ns: 'rgba(143,181,159,.4)' },

  { id: 'mono', name: 'Black Box', bg: 'aarain', layout: 'l-flat l-brutal l-uc',
    vars: {
      '--bg1': '#000', '--bg2': '#0c0c0c', '--ink': '#fff', '--mute': '#9b9b9b', '--lilac': '#d8d8d8',
      '--a': '#c6ff3d', '--b': '#ff3d3d', '--c': '#ffe14d', '--d': '#6bd6ff', '--e': '#c6ff3d',
      '--a-rgb': '198,255,61', '--b-rgb': '255,61,61', '--c-rgb': '255,225,77', '--d-rgb': '107,214,255', '--line-rgb': '255,255,255',
      '--up': '#ff6b3d', '--down': '#6bd6ff', '--card1': 'rgba(16,16,16,.92)', '--card2': 'rgba(6,6,6,.92)', '--header-bg': 'rgba(0,0,0,.92)',
      '--input-bg': '#000', '--thead-bg': '#111', '--code-bg': 'rgba(198,255,61,.14)', '--code-ink': '#c6ff3d', '--tip-bg': '#000',
      '--detail-bg': 'rgba(0,0,0,.92)', '--veil': 'linear-gradient(180deg,rgba(0,0,0,.15),rgba(0,0,0,.4))', '--btn-ink': '#000',
      '--btn-alt-bg': '#fff', '--btn-alt-ink': '#000', '--h2-grad': 'linear-gradient(#fff,#fff)', '--num-grad': 'linear-gradient(#c6ff3d,#c6ff3d)',
      '--body-ink': '#d9d9d9', '--toggle-off': '#333', '--shadow': 'none',
      '--font-sans': 'ui-monospace,"SF Mono",Menlo,monospace', '--font-head': 'ui-monospace,"SF Mono",Menlo,monospace',
      '--head-weight': '700', '--head-case': 'uppercase', '--head-track': '0', '--radius': '0px', '--radius-sm': '0px', '--pill': '0px',
      '--grad': 'linear-gradient(#c6ff3d,#c6ff3d)' },
    pal: ['#c6ff3d', '#ff3d3d', '#ffe14d', '#6bd6ff', '#ffffff', '#ff9d3d', '#d38bff', '#9b9b9b'],
    heat: [[8, 8, 8], [110, 140, 40], [198, 255, 61]], div: [[14, 14, 14], [255, 107, 61], [107, 214, 255]], ns: 'rgba(155,155,155,.4)' },

  { id: 'rose', name: 'Rose Quartz', bg: 'aurora', layout: 'l-float', light: true,
    vars: {
      '--bg1': '#fff3f8', '--bg2': '#ece6ff', '--ink': '#2b1b3d', '--mute': '#7b6a93', '--lilac': '#6b4fd0',
      '--a': '#7c5cff', '--b': '#ff5fa2', '--c': '#ffa94d', '--d': '#4cc3ff', '--e': '#5fd39a',
      '--a-rgb': '124,92,255', '--b-rgb': '255,95,162', '--c-rgb': '255,169,77', '--d-rgb': '76,195,255', '--line-rgb': '124,92,255',
      '--up': '#ff6b6b', '--down': '#5b8bff', '--card1': 'rgba(255,255,255,.74)', '--card2': 'rgba(255,244,250,.74)', '--header-bg': 'rgba(255,255,255,.72)',
      '--input-bg': 'rgba(255,255,255,.88)', '--thead-bg': '#f3ecff', '--code-bg': 'rgba(255,95,162,.12)', '--code-ink': '#c0306f',
      '--tip-bg': 'rgba(255,255,255,.97)', '--detail-bg': 'rgba(250,240,255,.92)',
      '--veil': 'linear-gradient(180deg,transparent,rgba(255,243,248,.2))', '--btn-ink': '#fff', '--btn-alt-ink': '#3a1a00',
      '--h2-grad': 'linear-gradient(95deg,#2b1b3d 10%,#7c5cff 60%,#ff5fa2)', '--num-grad': 'linear-gradient(95deg,#2b1b3d,#7c5cff)',
      '--body-ink': '#4a3c60', '--toggle-off': '#d9cff5',
      '--shadow': '0 24px 50px -28px rgba(124,92,255,.45),inset 0 1px 0 rgba(255,255,255,.8)', '--radius': '26px', '--radius-sm': '14px' },
    pal: ['#7c5cff', '#ff5fa2', '#ffa94d', '#4cc3ff', '#5fd39a', '#f87171', '#c084fc', '#94a3b8'],
    heat: [[255, 245, 250], [124, 92, 255], [255, 95, 162]], div: [[255, 245, 250], [255, 107, 107], [91, 139, 255]], ns: 'rgba(123,106,147,.35)' },

  { id: 'copper', name: 'Copper Helix', bg: 'helix', layout: 'l-center',
    vars: {
      '--bg1': '#0e0e11', '--bg2': '#1d1712', '--ink': '#f4ede3', '--mute': '#ab9c8a', '--lilac': '#e3c4a6',
      '--a': '#d9895b', '--b': '#e8b923', '--c': '#8fb8a8', '--d': '#7fa6d9', '--e': '#c7d96a',
      '--a-rgb': '217,137,91', '--b-rgb': '232,185,35', '--c-rgb': '143,184,168', '--d-rgb': '127,166,217', '--line-rgb': '217,160,110',
      '--up': '#e0734a', '--down': '#6aa0d9', '--card1': 'rgba(36,29,24,.82)', '--card2': 'rgba(18,15,13,.86)', '--header-bg': 'rgba(20,17,14,.84)',
      '--input-bg': 'rgba(12,10,8,.75)', '--thead-bg': '#241d17', '--code-bg': 'rgba(217,137,91,.16)', '--code-ink': '#f0b48c',
      '--tip-bg': 'rgba(14,12,10,.96)', '--detail-bg': 'rgba(10,8,6,.8)',
      '--veil': 'radial-gradient(ellipse at 50% 0%,rgba(14,14,17,.15),rgba(14,14,17,.62) 72%)', '--btn-ink': '#1a0f08',
      '--btn-alt-bg': 'linear-gradient(135deg,#bfd9cd,#8fb8a8)', '--btn-alt-ink': '#0f1d17',
      '--h2-grad': 'linear-gradient(95deg,#fff4e4 10%,#e9c9a6 55%,#d9895b)', '--num-grad': 'linear-gradient(95deg,#fff4e4,#e9c9a6)',
      '--body-ink': '#dccfbf', '--toggle-off': '#3a2f25',
      '--font-head': '"Didot","Bodoni 72","Iowan Old Style",Georgia,serif', '--head-weight': '500', '--head-track': '-.01em',
      '--radius': '10px', '--radius-sm': '8px', '--grad': 'linear-gradient(135deg,#e8b923,#d9895b)' },
    pal: ['#d9895b', '#e8b923', '#8fb8a8', '#7fa6d9', '#c7d96a', '#e0734a', '#d48ad9', '#a39a8f'],
    heat: [[22, 18, 14], [217, 137, 91], [255, 226, 150]], div: [[24, 19, 15], [224, 115, 74], [106, 160, 217]], ns: 'rgba(171,156,138,.4)' },

  { id: 'skylab', name: 'Sky Lab', bg: 'waves', layout: 'l-flat', light: true,
    vars: {
      '--bg1': '#f4f9ff', '--bg2': '#dcebff', '--ink': '#0f2342', '--mute': '#5b7090', '--lilac': '#1d4ed8',
      '--a': '#2563eb', '--b': '#ff6b6b', '--c': '#f5a524', '--d': '#0ea5e9', '--e': '#22c55e',
      '--a-rgb': '37,99,235', '--b-rgb': '255,107,107', '--c-rgb': '245,165,36', '--d-rgb': '14,165,233', '--line-rgb': '37,99,235',
      '--up': '#ef4444', '--down': '#2563eb', '--card1': 'rgba(255,255,255,.9)', '--card2': 'rgba(240,247,255,.9)', '--header-bg': 'rgba(255,255,255,.88)',
      '--input-bg': '#fff', '--thead-bg': '#e8f1ff', '--code-bg': 'rgba(37,99,235,.1)', '--code-ink': '#1d4ed8', '--tip-bg': 'rgba(255,255,255,.98)',
      '--detail-bg': 'rgba(235,244,255,.95)', '--veil': 'linear-gradient(180deg,transparent,rgba(244,249,255,.15))', '--btn-ink': '#fff',
      '--btn-alt-ink': '#3a2400', '--h2-grad': 'linear-gradient(#0f2342,#0f2342)', '--num-grad': 'linear-gradient(95deg,#0f2342,#2563eb)',
      '--body-ink': '#2a3f5f', '--toggle-off': '#c7d9f5',
      '--shadow': '0 18px 36px -24px rgba(15,35,66,.4),0 0 0 1px rgba(37,99,235,.06)',
      '--font-sans': '"Helvetica Neue","Avenir Next",Arial,sans-serif', '--font-head': '"Helvetica Neue","Avenir Next",Arial,sans-serif',
      '--head-weight': '800', '--head-track': '-.03em', '--radius': '16px', '--radius-sm': '10px', '--grad': 'linear-gradient(135deg,#2563eb,#0ea5e9)' },
    pal: ['#2563eb', '#ff6b6b', '#f5a524', '#0ea5e9', '#22c55e', '#ef4444', '#8b5cf6', '#64748b'],
    heat: [[236, 244, 255], [14, 165, 233], [255, 107, 107]], div: [[244, 249, 255], [239, 68, 68], [37, 99, 235]], ns: 'rgba(91,112,144,.35)' },

  { id: 'sunset', name: 'Spectrum Sunset', bg: 'mountains', layout: 'l-float',
    vars: {
      '--bg1': '#1b0b33', '--bg2': '#6a1b5d', '--ink': '#fff3ea', '--mute': '#e6b8c8', '--lilac': '#ffd1b0',
      '--a': '#ffb347', '--b': '#ff4f8b', '--c': '#ffe08a', '--d': '#8ec5ff', '--e': '#b6f27a',
      '--a-rgb': '255,179,71', '--b-rgb': '255,79,139', '--c-rgb': '255,224,138', '--d-rgb': '142,197,255', '--line-rgb': '255,190,150',
      '--up': '#ffa14a', '--down': '#7fb2ff', '--card1': 'rgba(60,20,70,.72)', '--card2': 'rgba(30,10,45,.76)', '--header-bg': 'rgba(36,12,56,.72)',
      '--input-bg': 'rgba(22,8,36,.7)', '--thead-bg': '#3a1450', '--code-bg': 'rgba(255,79,139,.18)', '--code-ink': '#ffb3d0',
      '--tip-bg': 'rgba(26,8,40,.96)', '--detail-bg': 'rgba(20,6,32,.78)',
      '--veil': 'linear-gradient(180deg,rgba(27,11,51,.12),rgba(27,11,51,.12) 60%,rgba(27,11,51,.3))', '--btn-ink': '#2a0d2a',
      '--h2-grad': 'linear-gradient(95deg,#fff3ea 10%,#ffd1a0 50%,#ff7ab0)', '--num-grad': 'linear-gradient(95deg,#fff3ea,#ffd1a0)',
      '--body-ink': '#f0d6df', '--toggle-off': '#5a2a6a', '--radius': '22px', '--grad': 'linear-gradient(135deg,#ff4f8b,#ffb347)' },
    pal: ['#ff4f8b', '#ffb347', '#ffe08a', '#8ec5ff', '#b6f27a', '#f87171', '#e879f9', '#c9a9b8'],
    heat: [[40, 12, 60], [255, 79, 139], [255, 224, 138]], div: [[40, 14, 60], [255, 161, 74], [127, 178, 255]], ns: 'rgba(230,184,200,.4)' },
];

const rgbOf = s => { const m = (s || '').match(/\d+/g); return m ? m.slice(0, 3).map(Number) : [255, 255, 255]; };
const hexRgb = h => { h = h.replace('#', ''); if (h.length === 3) h = h.split('').map(x => x + x).join(''); return [0, 2, 4].map(i => parseInt(h.substr(i, 2), 16)); };

function applyLook(id, quiet) {
  const L = LOOKS.find(l => l.id === id) || LOOKS[0], root = document.documentElement;
  (applyLook.set || []).forEach(k => root.style.removeProperty(k));
  applyLook.set = Object.keys(L.vars);
  applyLook.set.forEach(k => root.style.setProperty(k, L.vars[k]));
  document.body.className = L.layout;
  const css = getComputedStyle(root), v = k => css.getPropertyValue(k).trim();
  Object.assign(TH, {
    id: L.id, name: L.name, light: !!L.light, bgMode: L.bg, pal: L.pal, heat: L.heat, div: L.div, ns: L.ns,
    ink: v('--ink'), mute: v('--mute'), a: v('--a'), b: v('--b'), c: v('--c'), d: v('--d'), e: v('--e'),
    up: v('--up'), down: v('--down'), bg1: v('--bg1'), bg2: v('--bg2'),
    aRGB: rgbOf(v('--a-rgb')), bRGB: rgbOf(v('--b-rgb')), cRGB: rgbOf(v('--c-rgb')), dRGB: rgbOf(v('--d-rgb')), lineRGB: rgbOf(v('--line-rgb')),
    bg1RGB: hexRgb(v('--bg1')), bg2RGB: hexRgb(v('--bg2')),
    grid: `rgba(${v('--line-rgb')},${L.light ? .18 : .14})`, border: `rgba(${v('--line-rgb')},.4)`, cAlpha: `rgba(${v('--c-rgb')},.6)`,
  });
  const sel = document.getElementById('look'); if (sel) sel.value = L.id;
  try { localStorage.setItem('rd-look', L.id); } catch (e) {}
  if (!quiet && window.onLookChange) window.onLookChange();
}
(function () {
  const sel = document.getElementById('look');
  if (sel) { sel.innerHTML = LOOKS.map(l => `<option value="${l.id}">${l.name}</option>`).join(''); sel.onchange = () => applyLook(sel.value); }
  let id = 'paper';
  try { id = localStorage.getItem('rd-look') || id; } catch (e) {}
  const q = new URLSearchParams(location.search).get('look'); if (q) id = q;
  applyLook(id, true);
})();
