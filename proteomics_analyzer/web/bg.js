/* Animated background: a live mass spectrum. Stick peaks breathe, a detector sweep lights the
   peaks it passes, isotope envelopes drift, and chromatographic traces slide underneath. */
(function () {
  const cv = document.getElementById('bg'), c = cv.getContext('2d');
  const reduce = matchMedia('(prefers-reduced-motion: reduce)').matches;
  let W = 0, H = 0, peaks = [], traces = [], running = true, last = 0;
  const rnd = (a, b) => a + Math.random() * (b - a);

  function build() {
    peaks = [];
    const n = Math.max(40, Math.floor(W / 11));
    for (let i = 0; i < n; i++) {
      const x = rnd(0.01, 0.99), base = Math.pow(Math.random(), 2.4) * 0.8 + 0.04;
      const env = Math.random() < 0.45 ? 4 : 1;               // isotope cluster
      for (let k = 0; k < env; k++)
        peaks.push({ x: x + k * 0.0045, base: base * (k === 0 ? 1 : [0, .78, .38, .14][k]),
                     ph: rnd(0, 6.28), sp: rnd(.25, .9), amp: rnd(.05, .2) });
    }
    traces = Array.from({ length: 4 }, (_, i) => ({ x: rnd(0, 1), w: rnd(.03, .07), h: rnd(.08, .2), v: rnd(.006, .016) * (i % 2 ? 1 : .7) }));
  }
  function size() {
    const d = Math.min(devicePixelRatio || 1, 1.5);
    W = innerWidth; H = innerHeight; cv.width = W * d; cv.height = H * d; c.setTransform(d, 0, 0, d, 0, 0);
    build(); draw(performance.now());
  }
  function draw(t) {
    c.clearRect(0, 0, W, H);
    const g = c.createLinearGradient(0, 0, 0, H);
    g.addColorStop(0, '#04161b'); g.addColorStop(1, '#073038'); c.fillStyle = g; c.fillRect(0, 0, W, H);
    const base = H * 0.86, ph = H * 0.62, s = t / 1000;
    // m/z axis
    c.strokeStyle = 'rgba(94,234,212,.22)'; c.lineWidth = 1; c.beginPath(); c.moveTo(0, base); c.lineTo(W, base); c.stroke();
    c.fillStyle = 'rgba(94,234,212,.28)'; c.font = '10px ui-monospace,Menlo,monospace';
    for (let i = 0; i <= 12; i++) {
      const x = i / 12 * W; c.fillRect(x, base, 1, 6); if (i % 2 === 0) c.fillText(String(200 + i * 150), x + 4, base + 17);
    }
    // chromatographic traces
    traces.forEach((tr, i) => {
      if (!reduce) tr.x = (tr.x + tr.v * 0.016 + 1) % 1;
      c.beginPath();
      for (let px = 0; px <= W; px += 6) {
        const u = px / W, d = ((u - tr.x + 1.5) % 1) - 0.5;
        const y = base + 40 + 60 - Math.exp(-d * d / (2 * tr.w * tr.w)) * tr.h * 220 + Math.sin(u * 40 + i) * 1.2;
        px ? c.lineTo(px, y) : c.moveTo(px, y);
      }
      c.strokeStyle = `rgba(45,212,191,${.1 + i * .03})`; c.lineWidth = 1.4; c.stroke();
    });
    // detector sweep
    const sweep = reduce ? .42 : ((s * 0.045) % 1.15) - 0.075;
    // peaks
    for (const p of peaks) {
      const px = p.x * W, h = (p.base + (reduce ? 0 : Math.sin(s * p.sp + p.ph) * p.amp * p.base)) * ph;
      const near = Math.max(0, 1 - Math.abs(p.x - sweep) / 0.05);
      const a = .2 + near * .6;
      c.strokeStyle = near > .15 ? `rgba(251,191,36,${a})` : `rgba(94,234,212,${.16 + p.base * .35})`;
      c.lineWidth = 1.6; c.beginPath(); c.moveTo(px, base); c.lineTo(px, base - h); c.stroke();
      if (near > .5 && h > ph * 0.22) { c.fillStyle = `rgba(251,191,36,${near * .8})`; c.beginPath(); c.arc(px, base - h, 2.2, 0, 6.3); c.fill(); }
    }
    const sx = sweep * W, sg = c.createLinearGradient(sx - 40, 0, sx + 2, 0);
    sg.addColorStop(0, 'rgba(45,212,191,0)'); sg.addColorStop(1, 'rgba(45,212,191,.12)');
    c.fillStyle = sg; c.fillRect(sx - 40, 0, 42, base);
    c.fillStyle = 'rgba(94,234,212,.5)'; c.fillRect(sx, 0, 1, base);
  }
  function loop(t) {
    if (!running) return;
    if (t - last > 33) { last = t; draw(t); }          // ~30 fps is plenty for a backdrop
    requestAnimationFrame(loop);
  }
  window.setBackgroundAnimation = function (on) {
    running = on && !reduce;
    if (running) requestAnimationFrame(loop); else draw(performance.now());
  };
  addEventListener('resize', size);
  document.addEventListener('visibilitychange', () => { if (!document.hidden) window.setBackgroundAnimation(document.getElementById('anim').checked); });
  size();
  window.setBackgroundAnimation(true);
})();
