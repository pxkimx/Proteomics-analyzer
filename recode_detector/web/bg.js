/* Animated backgrounds, one per look. All are drawn from the active look's colours (TH, see themes.js).
   Each mode has build(W,H) (set up its state) and draw(s,W,H) (draw one frame at time s seconds). */
(function () {
  const cv = document.getElementById('bg'), c = cv.getContext('2d');
  const reduce = matchMedia('(prefers-reduced-motion: reduce)').matches;
  let W = 0, H = 0, running = true, last = 0, mode = null, S = {};
  const rnd = (a, b) => a + Math.random() * (b - a);
  const rgba = (rgb, a) => `rgba(${rgb[0]},${rgb[1]},${rgb[2]},${a})`;
  const mix = (p, q, t) => p.map((v, i) => Math.round(v + (q[i] - v) * t));
  const base = () => { const g = c.createLinearGradient(0, 0, W * .6, H); g.addColorStop(0, TH.bg1); g.addColorStop(1, TH.bg2); c.globalCompositeOperation = 'source-over'; c.fillStyle = g; c.fillRect(0, 0, W, H); };
  const glow = () => TH.light ? 'source-over' : 'lighter';

  const MODES = {
    /* ---- Ultraviolet: LC-MS feature map with an acquisition sweep ---- */
    lcms: {
      build() {
        S.f = []; const n = Math.max(90, Math.floor(W * H / 7000));
        for (let i = 0; i < n; i++) S.f.push({ x: rnd(.02, .98), y: rnd(.08, .86), I: Math.pow(Math.random(), 2.2), env: Math.random() < .6 ? 3 + Math.floor(Math.random() * 2) : 1, w: rnd(16, 64), ph: rnd(0, 6.28), sp: rnd(.2, .7), z: 1 + Math.floor(Math.random() * 3) });
        S.tr = new Float32Array(Math.ceil(W / 3) + 1);
        for (const f of S.f) { const cx = f.x * S.tr.length, sd = Math.max(2, f.w / 6); for (let k = Math.max(0, Math.floor(cx - 4 * sd)); k < Math.min(S.tr.length, cx + 4 * sd); k++) S.tr[k] = Math.max(S.tr[k], f.I * Math.exp(-((k - cx) ** 2) / (2 * sd * sd))); }
      },
      draw(s) {
        base(); const RUN = 34, sweep = reduce ? .62 : ((s % RUN) / RUN) * 1.16 - .04, L = 54, B = H - 56, T = 70, mh = B - T, ww = W - L - 20;
        c.strokeStyle = rgba(TH.lineRGB, .12); c.lineWidth = 1; c.beginPath();
        for (let i = 0; i <= 12; i++) { const x = L + ww * i / 12; c.moveTo(x, T); c.lineTo(x, B); }
        for (let j = 0; j <= 7; j++) { const y = T + mh * j / 7; c.moveTo(L, y); c.lineTo(W - 20, y); } c.stroke();
        c.fillStyle = rgba(TH.lineRGB, .45); c.font = '10px ui-monospace,Menlo,monospace'; c.textAlign = 'center';
        for (let i = 0; i <= 12; i += 2) c.fillText(String(i * 10), L + ww * i / 12, B + 16);
        c.fillText('retention time (min)', L + ww / 2, B + 30); c.textAlign = 'right';
        for (let j = 0; j <= 7; j++) c.fillText(String(1500 - j * 170), L - 8, T + mh * j / 7 + 3);
        c.save(); c.translate(14, T + mh / 2); c.rotate(-Math.PI / 2); c.textAlign = 'center'; c.fillText('m/z', 0, 0); c.restore();
        c.globalCompositeOperation = glow();
        const stops = [TH.aRGB, TH.bRGB, TH.cRGB];
        const heat = (t, a) => { const k = Math.min(.999, Math.max(0, t)) * 2, i = Math.floor(k); return rgba(mix(stops[i], stops[i + 1], k - i), a); };
        for (const f of S.f) {
          const px = L + f.x * ww, py = T + (1 - f.y) * mh, lit = f.x <= sweep ? Math.min(1, (sweep - f.x) * 14 + .3) : .2, near = Math.max(0, 1 - Math.abs(f.x - sweep) / .03);
          const pulse = reduce ? 1 : .85 + .15 * Math.sin(s * f.sp + f.ph), a = (.16 + f.I * .75) * lit * pulse + near * .5 * (f.I + .25);
          for (let e = 0; e < f.env; e++) {
            const w = f.w * (1 - e * .08), y = py - e * (7 / Math.max(1, f.z)) * 1.6, ai = a * [1, .8, .46, .2][e], g = c.createRadialGradient(px, y, 0, px, y, w);
            g.addColorStop(0, heat(Math.min(1, f.I * 1.3 + near * .3), Math.min(1, ai * 1.7))); g.addColorStop(1, heat(.1, 0));
            c.fillStyle = g; c.save(); c.translate(px, y); c.scale(1, .2); c.translate(-px, -y); c.beginPath(); c.arc(px, y, w, 0, 6.3); c.fill(); c.restore();
          }
        }
        c.globalCompositeOperation = 'source-over';
        const last = Math.min(S.tr.length - 1, Math.floor(sweep * S.tr.length)); c.beginPath(); c.moveTo(L, H - 30);
        for (let k = 0; k <= last; k++) c.lineTo(L + k / S.tr.length * ww, H - 30 - S.tr[k] * 34);
        c.strokeStyle = rgba(TH.bRGB, .55); c.lineWidth = 1.4; c.stroke();
        const sx = L + sweep * ww; if (sx > L && sx < W - 20) { const sg = c.createLinearGradient(sx - 60, 0, sx, 0); sg.addColorStop(0, rgba(TH.aRGB, 0)); sg.addColorStop(1, rgba(TH.aRGB, .14)); c.fillStyle = sg; c.fillRect(sx - 60, T, 60, mh); c.fillStyle = rgba(TH.bRGB, .6); c.fillRect(sx, T, 1, mh + 28); }
      } },

    /* ---- Solar Flare: a stick spectrum that burns ---- */
    spectrum: {
      build() { S.p = []; const n = Math.floor(W / 7); for (let i = 0; i < n; i++) { const x = (i + .5) / n, env = Math.random() < .35 ? 4 : 1, b = Math.pow(Math.random(), 2.3) * .85 + .05; for (let k = 0; k < env; k++) S.p.push({ x: x + k * .0035, b: b * [1, .8, .45, .18][k], ph: rnd(0, 6.28), sp: rnd(.5, 1.6), a: rnd(.1, .3) }); } },
      draw(s) {
        base(); const bl = H * .9, ph = H * .7, sweep = reduce ? .5 : (s * .05) % 1.2 - .1;
        c.globalCompositeOperation = glow(); c.lineCap = 'round';
        for (const p of S.p) {
          const h = (p.b + (reduce ? 0 : Math.sin(s * p.sp + p.ph) * p.a * p.b)) * ph, x = p.x * W, near = Math.max(0, 1 - Math.abs(p.x - sweep) / .06);
          const g = c.createLinearGradient(0, bl, 0, bl - h); g.addColorStop(0, rgba(TH.bRGB, .05 + near * .4)); g.addColorStop(.5, rgba(TH.aRGB, .35 + near * .4)); g.addColorStop(1, rgba(TH.cRGB, .7 + near * .3));
          c.strokeStyle = g; c.lineWidth = 2 + p.b * 3; c.beginPath(); c.moveTo(x, bl); c.lineTo(x, bl - h); c.stroke();
        }
        c.globalCompositeOperation = 'source-over'; const g2 = c.createLinearGradient(0, bl - 120, 0, H); g2.addColorStop(0, rgba(TH.bRGB, 0)); g2.addColorStop(1, rgba(TH.bRGB, .18)); c.fillStyle = g2; c.fillRect(0, bl - 120, W, H);
        c.fillStyle = rgba(TH.cRGB, .5); c.font = '10px ui-monospace,Menlo,monospace'; c.textAlign = 'left'; for (let i = 0; i <= 12; i += 2) c.fillText(String(200 + i * 150), i / 12 * (W - 40) + 10, H - 12);
      } },

    /* ---- Paper Lab: ridge-line plot (a stack of drifting spectra) ---- */
    ridges: {
      build() { S.r = []; const n = 30; for (let i = 0; i < n; i++) { const pk = []; for (let k = 0; k < 4; k++) pk.push({ x: rnd(.1, .9), a: rnd(.3, 1), w: rnd(.015, .05), d: rnd(-.004, .004), ph: rnd(0, 6.28) }); S.r.push(pk); } },
      draw(s) {
        base(); const top = H * .3, bot = H * .94, n = S.r.length, step = 5;
        for (let i = 0; i < n; i++) {
          const y0 = top + (bot - top) * i / (n - 1), pk = S.r[i];
          c.beginPath(); c.moveTo(0, H);
          const pts = [];
          for (let x = 0; x <= W; x += step) {
            let v = 0; for (const p of pk) { const px = p.x + (reduce ? 0 : Math.sin(s * .15 + p.ph) * .05 + p.d * s), d = (x / W - px) / p.w; v += p.a * Math.exp(-d * d / 2); }
            const env = Math.exp(-(((x / W) - .5) ** 2) / .16); pts.push([x, y0 - v * env * H * .16 - 2 * Math.sin(x * .05 + i + (reduce ? 0 : s))]);
          }
          c.lineTo(0, pts[0][1]); for (const [x, y] of pts) c.lineTo(x, y); c.lineTo(W, H); c.closePath();
          c.fillStyle = rgba(TH.bg1RGB, .96); c.fill();
          c.beginPath(); pts.forEach(([x, y], k) => k ? c.lineTo(x, y) : c.moveTo(x, y));
          c.strokeStyle = i === Math.floor(n * .55) ? rgba(TH.bRGB, .85) : rgba(TH.lineRGB, .45); c.lineWidth = i === Math.floor(n * .55) ? 1.6 : 1; c.stroke();
        }
      } },

    /* ---- Neon Grid: perspective grid under a striped sun ---- */
    grid: {
      build() {},
      draw(s) {
        base(); const hz = H * .58, off = reduce ? .3 : (s * .22) % 1;
        const sky = c.createLinearGradient(0, 0, 0, hz); sky.addColorStop(0, rgba(TH.bg1RGB, 0)); sky.addColorStop(1, rgba(TH.bRGB, .28)); c.fillStyle = sky; c.fillRect(0, 0, W, hz);
        const sr = Math.min(W, H) * .2, sy = hz - sr * .55, sg = c.createLinearGradient(0, sy - sr, 0, sy + sr); sg.addColorStop(0, TH.c); sg.addColorStop(1, TH.b);
        c.save(); c.beginPath(); c.arc(W / 2, sy, sr, 0, 6.3); c.clip(); c.fillStyle = sg; c.fillRect(W / 2 - sr, sy - sr, sr * 2, sr * 2);
        c.fillStyle = TH.bg1; for (let k = 0; k < 7; k++) { const y = sy + sr * (.05 + k * .14) + ((s * 6) % 6) * .0, th = 1.5 + k * 1.6; c.fillRect(W / 2 - sr, y, sr * 2, th); } c.restore();
        c.fillStyle = rgba(TH.bg1RGB, .92); c.fillRect(0, hz, W, H - hz);
        c.strokeStyle = rgba(TH.aRGB, .8); c.shadowColor = TH.a; c.shadowBlur = 8; c.lineWidth = 1.2;
        for (let k = 0; k < 14; k++) { const z = (k + off) / 14, y = hz + (H - hz) * z * z; c.globalAlpha = .25 + z * .75; c.beginPath(); c.moveTo(0, y); c.lineTo(W, y); c.stroke(); }
        c.globalAlpha = 1; for (let i = -16; i <= 16; i++) { c.beginPath(); c.moveTo(W / 2 + i * 6, hz); c.lineTo(W / 2 + i * W * .09, H); c.stroke(); }
        c.shadowBlur = 0; c.strokeStyle = TH.b; c.lineWidth = 1.5; c.beginPath(); c.moveTo(0, hz); c.lineTo(W, hz); c.stroke();
      } },

    /* ---- Mycelium: a drifting protein-interaction network with signals ---- */
    network: {
      build() { const n = Math.min(110, Math.floor(W * H / 14000)); S.n = Array.from({ length: n }, () => ({ x: rnd(0, W), y: rnd(0, H), vx: rnd(-.12, .12), vy: rnd(-.1, .1), r: Math.random() < .12 ? rnd(4, 6.5) : rnd(1.4, 3), ph: rnd(0, 6.28) })); S.pulse = []; S.t = 0; },
      draw(s) {
        base(); const D = 170; c.lineWidth = 1;
        for (const p of S.n) { if (!reduce) { p.x += p.vx; p.y += p.vy; if (p.x < -20) p.x = W + 20; if (p.x > W + 20) p.x = -20; if (p.y < -20) p.y = H + 20; if (p.y > H + 20) p.y = -20; } }
        const edges = [];
        for (let i = 0; i < S.n.length; i++) for (let j = i + 1; j < S.n.length; j++) { const a = S.n[i], b = S.n[j], d = Math.hypot(a.x - b.x, a.y - b.y); if (d < D) { edges.push([a, b]); c.strokeStyle = rgba(TH.lineRGB, (1 - d / D) * .55); c.beginPath(); c.moveTo(a.x, a.y); c.lineTo(b.x, b.y); c.stroke(); } }
        if (!reduce && edges.length && s - S.t > .35) { S.t = s; S.pulse.push({ e: edges[Math.floor(Math.random() * edges.length)], k: 0 }); if (S.pulse.length > 14) S.pulse.shift(); }
        for (const p of S.n) { const r = p.r * (reduce ? 1 : .9 + .1 * Math.sin(s * 1.3 + p.ph)); if (p.r > 4) { const g = c.createRadialGradient(p.x, p.y, 0, p.x, p.y, r * 5); g.addColorStop(0, rgba(TH.aRGB, .35)); g.addColorStop(1, rgba(TH.aRGB, 0)); c.fillStyle = g; c.beginPath(); c.arc(p.x, p.y, r * 5, 0, 6.3); c.fill(); } c.fillStyle = rgba(p.r > 4 ? TH.bRGB : TH.aRGB, .85); c.beginPath(); c.arc(p.x, p.y, r, 0, 6.3); c.fill(); }
        S.pulse = S.pulse.filter(q => (q.k += .02) < 1); for (const q of S.pulse) { const [a, b] = q.e, x = a.x + (b.x - a.x) * q.k, y = a.y + (b.y - a.y) * q.k; c.fillStyle = rgba(TH.bRGB, 1 - Math.abs(q.k - .5)); c.shadowColor = TH.b; c.shadowBlur = 10; c.beginPath(); c.arc(x, y, 2.6, 0, 6.3); c.fill(); c.shadowBlur = 0; }
      } },

    /* ---- Black Box: falling amino-acid letters; R flips to W ---- */
    aarain: {
      persist: true,
      build() { S.cols = []; S.cw = 22; const n = Math.ceil(W / S.cw); for (let i = 0; i < n; i++) S.cols.push({ y: rnd(-H, 0), v: rnd(30, 90), row: -1 }); c.globalCompositeOperation = 'source-over'; c.fillStyle = TH.bg1; c.fillRect(0, 0, W, H); S.last = 0; },
      draw(s) {
        const dt = Math.min(.1, s - S.last || .04); S.last = s; c.globalCompositeOperation = 'source-over';
        c.fillStyle = rgba(TH.bg1RGB, .09); c.fillRect(0, 0, W, H); c.font = '700 14px ui-monospace,Menlo,monospace'; c.textAlign = 'center';
        const AA = 'ACDEFGHIKLMNPQRSTVYR';
        S.cols.forEach((col, i) => {
          col.y += col.v * dt * (reduce ? 0 : 1); const row = Math.floor(col.y / 20);
          if (row !== col.row && col.y > 0) {
            col.row = row; let ch = AA[Math.floor(Math.random() * AA.length)], hot = false;
            if (ch === 'R' && Math.random() < .55) { ch = 'W'; hot = true; }
            c.fillStyle = hot ? TH.b : rgba(TH.aRGB, .9); if (hot) { c.shadowColor = TH.b; c.shadowBlur = 12; }
            c.fillText(ch, i * S.cw + S.cw / 2, row * 20); c.shadowBlur = 0;
          }
          if (col.y > H + 40) { col.y = rnd(-300, 0); col.row = -1; col.v = rnd(30, 90); }
        });
      } },

    /* ---- Rose Quartz: slow aurora blobs ---- */
    aurora: {
      build() { S.b = [0, 1, 2, 3, 4].map(i => ({ ph: rnd(0, 6.28), fx: rnd(.05, .12), fy: rnd(.04, .1), cx: rnd(.15, .85), cy: rnd(.2, .8), r: rnd(.28, .45), col: [TH.aRGB, TH.bRGB, TH.cRGB, TH.dRGB, TH.aRGB][i] })); },
      draw(s) {
        base(); c.globalCompositeOperation = glow();
        for (const b of S.b) { const x = (b.cx + (reduce ? 0 : Math.sin(s * b.fx + b.ph) * .18)) * W, y = (b.cy + (reduce ? 0 : Math.cos(s * b.fy + b.ph) * .15)) * H, r = b.r * Math.max(W, H) * .8, g = c.createRadialGradient(x, y, 0, x, y, r); g.addColorStop(0, rgba(b.col, TH.light ? .5 : .42)); g.addColorStop(1, rgba(b.col, 0)); c.fillStyle = g; c.fillRect(0, 0, W, H); }
        c.globalCompositeOperation = 'source-over'; c.fillStyle = rgba(TH.lineRGB, .08); for (let i = 0; i < 70; i++) { const x = (i * 97.3) % W, y = (i * 53.7) % H; c.beginPath(); c.arc(x, y, 1.2, 0, 6.3); c.fill(); }
      } },

    /* ---- Copper Helix: DNA strands carrying codons ---- */
    helix: {
      build() { S.h = [{ cy: .3, A: 62, k: .012, sp: .7, sc: 1 }, { cy: .62, A: 46, k: .017, sp: -.5, sc: .8 }, { cy: .86, A: 34, k: .022, sp: .9, sc: .6 }]; },
      draw(s) {
        base(); c.globalCompositeOperation = glow(); const L = 'ACGT', cod = ['CGG', 'TGG'];
        S.h.forEach((h, hi) => {
          const cy = h.cy * H, ph = reduce ? 0 : s * h.sp;
          for (let x = 0; x < W; x += 5) {
            const t1 = x * h.k + ph, y1 = cy + Math.sin(t1) * h.A, y2 = cy - Math.sin(t1) * h.A, d = Math.cos(t1);
            c.fillStyle = rgba(TH.aRGB, .28 + .5 * (d + 1) / 2 * h.sc); c.beginPath(); c.arc(x, y1, 1.6 + d * .8 * h.sc, 0, 6.3); c.fill();
            c.fillStyle = rgba(TH.cRGB, .28 + .5 * (1 - d) / 2 * h.sc); c.beginPath(); c.arc(x, y2, 1.6 - d * .8 * h.sc, 0, 6.3); c.fill();
            if (x % 25 === 0) { c.strokeStyle = rgba(TH.lineRGB, .18 * h.sc); c.lineWidth = 1; c.beginPath(); c.moveTo(x, y1); c.lineTo(x, y2); c.stroke(); if (Math.abs(d) < .35) { c.fillStyle = rgba(TH.bRGB, .55); c.font = '10px ui-monospace,Menlo,monospace'; c.textAlign = 'center'; c.fillText(L[(x / 25 + hi) % 4 | 0], x, cy + 3); } }
          }
        });
        c.globalCompositeOperation = 'source-over'; c.fillStyle = rgba(TH.bRGB, .5); c.font = '700 11px ui-monospace,Menlo,monospace'; c.textAlign = 'left';
        const t = reduce ? 0 : s; c.fillText(Math.floor(t / 3) % 2 ? 'TGG  Trp' : 'CGG  Arg', W * .72, H * .2);
      } },

    /* ---- Sky Lab: layered signal waves ---- */
    waves: {
      build() {},
      draw(s) {
        base(); const cols = [TH.aRGB, TH.dRGB, TH.bRGB, TH.aRGB, TH.dRGB];
        for (let i = 0; i < 5; i++) {
          const by = H * (.6 + i * .075), A = 14 + i * 7, k = .004 + i * .0016, sp = (reduce ? 0 : s) * (.5 + i * .12) * (i % 2 ? -1 : 1);
          c.beginPath(); c.moveTo(0, H); for (let x = 0; x <= W; x += 6) c.lineTo(x, by + Math.sin(x * k + sp) * A + Math.sin(x * k * 2.3 + sp * 1.4) * A * .35);
          c.lineTo(W, H); c.closePath(); c.fillStyle = rgba(cols[i], .1 + i * .02); c.fill();
          c.beginPath(); for (let x = 0; x <= W; x += 6) { const y = by + Math.sin(x * k + sp) * A + Math.sin(x * k * 2.3 + sp * 1.4) * A * .35; x ? c.lineTo(x, y) : c.moveTo(x, y); } c.strokeStyle = rgba(cols[i], .35); c.lineWidth = 1.3; c.stroke();
        }
        for (let j = 0; j < 3; j++) { c.beginPath(); for (let x = 0; x <= W; x += 5) { const y = H * (.18 + j * .08) + Math.sin(x * .01 + (reduce ? 0 : s) * (.8 + j * .3)) * 16 + Math.sin(x * .027 - (reduce ? 0 : s) * 1.2) * 7; x ? c.lineTo(x, y) : c.moveTo(x, y); } c.strokeStyle = rgba(TH.aRGB, .16); c.lineWidth = 1; c.stroke(); }
      } },

    /* ---- Spectrum Sunset: mass-spec peaks as a mountain range ---- */
    mountains: {
      build() {
        S.layers = [0, 1, 2, 3].map(i => { const peaks = []; const n = 60 + i * 15; for (let k = 0; k < n; k++) peaks.push({ x: k / n, h: Math.pow(Math.random(), 1.6) * (.16 + .05 * i) + .03, w: rnd(.004, .014) }); return { peaks, sp: 4 + i * 9, y: .6 + i * .09 }; });
        S.stars = Array.from({ length: 70 }, () => ({ x: Math.random(), y: Math.random() * .45, p: rnd(0, 6.28) }));
      },
      draw(s) {
        c.globalCompositeOperation = 'source-over'; const hz = H * .72, sky = c.createLinearGradient(0, 0, 0, H);
        sky.addColorStop(0, TH.bg1); sky.addColorStop(.45, TH.bg2); sky.addColorStop(.74, rgba(TH.cRGB, 1)); sky.addColorStop(.78, rgba(TH.bRGB, 1)); sky.addColorStop(1, TH.bg1); c.fillStyle = sky; c.fillRect(0, 0, W, H);
        c.fillStyle = '#fff'; for (const st of S.stars) { c.globalAlpha = .3 + .5 * (reduce ? .5 : .5 + .5 * Math.sin(s * 1.4 + st.p)); c.fillRect(st.x * W, st.y * H, 1.3, 1.3); } c.globalAlpha = 1;
        const sx = W * .68, sy = hz - H * .02, sr = Math.min(W, H) * .13, sg = c.createRadialGradient(sx, sy, 0, sx, sy, sr * 2.6); sg.addColorStop(0, rgba(TH.cRGB, .95)); sg.addColorStop(.35, rgba(TH.bRGB, .5)); sg.addColorStop(1, rgba(TH.bRGB, 0)); c.fillStyle = sg; c.beginPath(); c.arc(sx, sy, sr * 2.6, 0, 6.3); c.fill();
        c.fillStyle = rgba(TH.cRGB, 1); c.beginPath(); c.arc(sx, sy, sr, 0, 6.3); c.fill();
        S.layers.forEach((L, i) => {
          const off = reduce ? 0 : (s * L.sp) % W, shade = mix([70, 22, 96], [14, 5, 30], i / 3), by = H * L.y + H * .1;
          c.fillStyle = `rgb(${shade[0]},${shade[1]},${shade[2]})`;
          for (let rep = -1; rep <= 1; rep++) {
            c.beginPath(); c.moveTo(rep * W - off, H);
            for (const p of L.peaks) { const x = rep * W + p.x * W - off; c.lineTo(x - p.w * W, by); c.lineTo(x, by - p.h * H); c.lineTo(x + p.w * W, by); }
            c.lineTo(rep * W + W - off, H); c.closePath(); c.fill();
          }
        });
      } },
  };

  function size() {
    const d = Math.min(devicePixelRatio || 1, 1.5);
    W = innerWidth; H = innerHeight; cv.width = W * d; cv.height = H * d; c.setTransform(d, 0, 0, d, 0, 0);
    start();
  }
  function start() {
    mode = MODES[TH.bgMode] || MODES.lcms; S = {}; c.shadowBlur = 0; c.globalAlpha = 1;
    mode.build(); frame(performance.now());
  }
  function frame(t) { c.shadowBlur = 0; c.globalAlpha = 1; mode.draw(t / 1000, W, H); }
  function loop(t) {
    if (!running) return;
    if (t - last > (mode.persist ? 45 : 40)) { last = t; frame(t); }
    requestAnimationFrame(loop);
  }
  window.setBackgroundAnimation = function (on) { running = on && !reduce; if (running) requestAnimationFrame(loop); else frame(performance.now()); };
  window.restartBackground = function () { start(); window.setBackgroundAnimation(document.getElementById('anim').checked); };
  addEventListener('resize', size);
  document.addEventListener('visibilitychange', () => { if (!document.hidden) window.setBackgroundAnimation(document.getElementById('anim').checked); });
  size();
  window.setBackgroundAnimation(true);
})();
