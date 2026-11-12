/**
 * charts.js
 * Publication-Grade Interactive Graphical Suite for FragPipe -> saintPy Studio
 * Includes:
 *  1. ProHits-style Bait x Prey DotPlot Matrix (Interactive crosshairs, selection ring + full Vector SVG export)
 *  2. Interactive Volcano & Significance Scatter Plot (Bait color legend, pill callouts + full Vector SVG export)
 *  3. Force-Directed Bait-Prey Interactome Network (Interactive physics, hub halos, drag + Vector SVG/PNG export)
 *  4. Replicate-Level Breakdown & SAINT Posterior Fit (eta_i, d_i, Z_ijr, P(Z=1|y)) Inspector (+ SVG/PNG export)
 *  5. Global SAINT Score Bimodal Distribution & ICM Log-Likelihood Convergence Plot (+ Vector SVG export)
 */

(function (root, factory) {
  root.SaintPyCharts = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  // --- Publication Color Palettes & Scales ---
  const BAIT_PALETTE = [
    "#2563eb", // Royal Blue
    "#059669", // Emerald
    "#d97706", // Amber
    "#7c3aed", // Violet
    "#db2777", // Rose Pink
    "#0891b2", // Cyan
    "#dc2626", // Crimson
    "#4f46e5", // Indigo
    "#65a30d", // Lime
    "#ea580c", // Deep Orange
    "#0d9488", // Teal
    "#9333ea", // Purple
  ];

  // Continuous color maps for SaintScore (0 -> 1)
  const COLORMAPS = {
    viridis: [
      [68, 1, 84],
      [59, 82, 139],
      [33, 145, 140],
      [94, 201, 98],
      [253, 231, 37],
    ],
    plasma: [
      [13, 8, 135],
      [126, 3, 168],
      [204, 71, 120],
      [248, 149, 64],
      [240, 249, 33],
    ],
    ylgnbu: [
      [237, 248, 177],
      [127, 205, 187],
      [44, 127, 184],
      [37, 52, 148],
      [8, 29, 88],
    ],
    crimson: [
      [254, 229, 217],
      [252, 174, 145],
      [251, 106, 74],
      [222, 45, 38],
      [165, 15, 21],
    ],
  };

  function interpolateColor(t, mapName = "plasma") {
    const stops = COLORMAPS[mapName] || COLORMAPS.plasma;
    const clamped = Math.max(0, Math.min(1, Number.isNaN(t) ? 0 : t));
    const scaled = clamped * (stops.length - 1);
    const idx = Math.min(stops.length - 2, Math.floor(scaled));
    const frac = scaled - idx;
    const c0 = stops[idx];
    const c1 = stops[idx + 1];
    const r = Math.round(c0[0] + frac * (c1[0] - c0[0]));
    const g = Math.round(c0[1] + frac * (c1[1] - c0[1]));
    const b = Math.round(c0[2] + frac * (c1[2] - c0[2]));
    return `rgb(${r}, ${g}, ${b})`;
  }

  function setupHiDPICanvas(canvas, width, height) {
    const dpr = Math.max(2, window.devicePixelRatio || 1);
    canvas.width = Math.round(width * dpr);
    canvas.height = Math.round(height * dpr);
    canvas.style.width = `${width}px`;
    canvas.style.height = `${height}px`;
    const ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return ctx;
  }

  function drawRoundedRect(ctx, x, y, w, h, r) {
    const radius = Math.min(r, w / 2, h / 2);
    ctx.beginPath();
    ctx.moveTo(x + radius, y);
    ctx.arcTo(x + w, y, x + w, y + h, radius);
    ctx.arcTo(x + w, y + h, x, y + h, radius);
    ctx.arcTo(x, y + h, x, y, radius);
    ctx.arcTo(x, y, x + w, y, radius);
    ctx.closePath();
  }

  function escapeXml(str) {
    return String(str ?? "")
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }

  function downloadBlob(content, filename, mimeType) {
    const blob = new Blob([content], { type: mimeType });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  function exportCanvasPNG(canvas, filename = "saintpy-figure.png") {
    const dataUrl = canvas.toDataURL("image/png");
    const a = document.createElement("a");
    a.href = dataUrl;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  }

  // ============================================================================
  // 1. PROHITS-STYLE BAIT x PREY DOTPLOT MATRIX
  // ============================================================================

  function prepareDotPlotMatrix(rows, options = {}) {
    const {
      minScore = 0.5,
      maxBfdr = 0.2,
      topNPerBait = 12,
      maxTotalPreys = 42,
      selectedBaits = null,
      preySearch = "",
    } = options;

    const allBaitsSet = new Set();
    for (const r of rows) allBaitsSet.add(r.Bait);
    const baits = selectedBaits && selectedBaits.length > 0 ? selectedBaits : Array.from(allBaitsSet);
    if (baits.length === 0) return { baits: [], preys: [], cellMap: new Map(), maxAbundance: 1 };

    const searchLower = preySearch.trim().toLowerCase();

    const chosenPreySet = new Set();
    const byBait = new Map(baits.map((b) => [b, []]));

    for (const r of rows) {
      if (!byBait.has(r.Bait)) continue;
      if (searchLower) {
        if (
          !r.PreyGene.toLowerCase().includes(searchLower) &&
          !r.Prey.toLowerCase().includes(searchLower)
        ) {
          continue;
        }
      }
      byBait.get(r.Bait).push(r);
    }

    for (const b of baits) {
      const list = byBait.get(b) || [];
      list.sort((a, z) => z.SaintScore - a.SaintScore || z.AvgSpec - a.AvgSpec);
      let count = 0;
      for (const item of list) {
        if (item.SaintScore >= minScore && (Number.isNaN(item.computedBFDR) || item.computedBFDR <= maxBfdr)) {
          chosenPreySet.add(item.Prey);
          count++;
          if (count >= topNPerBait) break;
        }
      }
    }

    // Fallback: if strict threshold yielded 0 preys, take top preys overall so chart is never blank
    if (chosenPreySet.size === 0) {
      const sortedAll = [...rows].sort((a, z) => z.SaintScore - a.SaintScore || z.AvgSpec - a.AvgSpec);
      for (const item of sortedAll) {
        if (searchLower) {
          if (
            !item.PreyGene.toLowerCase().includes(searchLower) &&
            !item.Prey.toLowerCase().includes(searchLower)
          ) {
            continue;
          }
        }
        chosenPreySet.add(item.Prey);
        if (chosenPreySet.size >= Math.min(20, maxTotalPreys)) break;
      }
    }

    let preysList = Array.from(chosenPreySet);
    const cellMap = new Map();
    const preyGeneMap = new Map();
    let maxAbundance = 0;

    for (const r of rows) {
      if (chosenPreySet.has(r.Prey) && byBait.has(r.Bait)) {
        cellMap.set(`${r.Bait}|||${r.Prey}`, r);
        preyGeneMap.set(r.Prey, r.PreyGene || r.Prey);
        if (r.AvgSpec > maxAbundance) maxAbundance = r.AvgSpec;
      }
    }

    // Sort preys by primary bait cluster, then shared bait count, then max SaintScore
    preysList.sort((pA, pB) => {
      let bestBaitA = 999,
        bestScoreA = -1,
        sigCountA = 0;
      let bestBaitB = 999,
        bestScoreB = -1,
        sigCountB = 0;
      baits.forEach((b, bIdx) => {
        const cA = cellMap.get(`${b}|||${pA}`);
        if (cA) {
          if (cA.SaintScore >= minScore) sigCountA++;
          if (cA.SaintScore > bestScoreA) {
            bestScoreA = cA.SaintScore;
            bestBaitA = bIdx;
          }
        }
        const cB = cellMap.get(`${b}|||${pB}`);
        if (cB) {
          if (cB.SaintScore >= minScore) sigCountB++;
          if (cB.SaintScore > bestScoreB) {
            bestScoreB = cB.SaintScore;
            bestBaitB = bIdx;
          }
        }
      });
      if (bestBaitA !== bestBaitB) return bestBaitA - bestBaitB;
      if (sigCountB !== sigCountA) return sigCountB - sigCountA;
      return bestScoreB - bestScoreA;
    });

    if (preysList.length > maxTotalPreys) {
      preysList = preysList.slice(0, maxTotalPreys);
    }

    return {
      baits,
      preys: preysList.map((pid) => {
        let sigCount = 0;
        for (const b of baits) {
          const c = cellMap.get(`${b}|||${pid}`);
          if (c && c.SaintScore >= minScore) sigCount++;
        }
        return {
          prey_id: pid,
          gene_name: preyGeneMap.get(pid) || pid,
          isShared: sigCount >= 2,
        };
      }),
      cellMap,
      maxAbundance: maxAbundance > 0 ? maxAbundance : 1,
    };
  }

  function renderDotPlot(canvas, rows, options = {}, callbacks = {}) {
    const matrix = prepareDotPlotMatrix(rows, options);
    const { baits, preys, cellMap, maxAbundance } = matrix;
    const colormap = options.colormap || "plasma";
    const sizeScale = options.sizeScale || "log";
    const selectedRow = options.selectedRow || null;

    const cellW = Math.max(46, Math.min(74, Math.floor(460 / Math.max(1, baits.length))));
    const cellH = 30;
    const leftPad = 148;
    const topPad = 102;
    const rightLegendW = 215;
    const bottomPad = 38;

    const plotW = Math.max(200, baits.length * cellW);
    const plotH = Math.max(140, preys.length * cellH);
    const width = leftPad + plotW + rightLegendW;
    const height = Math.max(350, topPad + plotH + bottomPad);

    const ctx = setupHiDPICanvas(canvas, width, height);

    if (preys.length === 0 || baits.length === 0) {
      ctx.fillStyle = "#ffffff";
      ctx.fillRect(0, 0, width, height);
      ctx.fillStyle = "#64748b";
      ctx.font = "500 14px Inter, system-ui, sans-serif";
      ctx.textAlign = "center";
      ctx.fillText("No bait-prey pairs match the current filter thresholds.", width / 2, height / 2);
      return { matrix, hitRegions: [], toSVG: () => "" };
    }

    const maxRadius = Math.min(cellW, cellH) * 0.42;
    const minRadius = 3.0;

    function calcRadius(val) {
      if (!val || val <= 0) return minRadius;
      let norm;
      if (sizeScale === "log") {
        norm = Math.log2(val + 1) / Math.log2(maxAbundance + 1);
      } else {
        norm = val / maxAbundance;
      }
      norm = Math.max(0.1, Math.min(1, norm));
      return minRadius + norm * (maxRadius - minRadius);
    }

    const hitRegions = [];
    for (let r = 0; r < preys.length; r++) {
      const cy = topPad + r * cellH + cellH / 2;
      for (let c = 0; c < baits.length; c++) {
        const cx = leftPad + c * cellW + cellW / 2;
        const cell = cellMap.get(`${baits[c]}|||${preys[r].prey_id}`);
        if (!cell || cell.AvgSpec <= 0) continue;
        const rad = calcRadius(cell.AvgSpec);
        hitRegions.push({
          cx,
          cy,
          r: Math.max(rad + 4, 10),
          rad,
          rowIdx: r,
          colIdx: c,
          row: cell,
        });
      }
    }

    let hoveredHit = null;

    function drawScene() {
      ctx.fillStyle = "#ffffff";
      ctx.fillRect(0, 0, width, height);

      // Title banner at top-left
      ctx.fillStyle = "#0f172a";
      ctx.font = "700 13px Inter, system-ui, sans-serif";
      ctx.textAlign = "left";
      ctx.textBaseline = "top";
      ctx.fillText("ProHits Bait × Prey Interaction Matrix", 16, 14);
      ctx.fillStyle = "#64748b";
      ctx.font = "500 11px Inter, system-ui, sans-serif";
      ctx.fillText(
        `${baits.length} Baits × ${preys.length} Top Preys · Color = SaintScore · Size = Abundance`,
        16,
        32
      );

      const gridW = baits.length * cellW;
      const gridH = preys.length * cellH;

      // Alternating row zebra stripes & hover highlight
      for (let r = 0; r < preys.length; r++) {
        const y = topPad + r * cellH;
        const isHoveredRow = hoveredHit && hoveredHit.rowIdx === r;
        if (isHoveredRow) {
          ctx.fillStyle = "rgba(37, 99, 235, 0.07)";
          ctx.fillRect(leftPad, y, gridW, cellH);
        } else if (r % 2 === 0) {
          ctx.fillStyle = "#f8fafc";
          ctx.fillRect(leftPad, y, gridW, cellH);
        }
        ctx.strokeStyle = isHoveredRow ? "#93c5fd" : "#e2e8f0";
        ctx.lineWidth = isHoveredRow ? 1.2 : 1;
        ctx.beginPath();
        ctx.moveTo(leftPad, y + cellH / 2);
        ctx.lineTo(leftPad + gridW, y + cellH / 2);
        ctx.stroke();
      }

      // Vertical bait grid lines & hover column highlight
      for (let c = 0; c < baits.length; c++) {
        const colLeft = leftPad + c * cellW;
        const x = colLeft + cellW / 2;
        const isHoveredCol = hoveredHit && hoveredHit.colIdx === c;
        if (isHoveredCol) {
          ctx.fillStyle = "rgba(37, 99, 235, 0.05)";
          ctx.fillRect(colLeft, topPad, cellW, gridH);
        }
        ctx.strokeStyle = isHoveredCol ? "#93c5fd" : "#e2e8f0";
        ctx.lineWidth = isHoveredCol ? 1.2 : 1;
        ctx.beginPath();
        ctx.moveTo(x, topPad);
        ctx.lineTo(x, topPad + gridH);
        ctx.stroke();
      }

      // Outer matrix border
      ctx.strokeStyle = "#cbd5e1";
      ctx.lineWidth = 1.4;
      ctx.strokeRect(leftPad, topPad, gridW, gridH);

      // Column labels (Baits, rotated -38 deg)
      ctx.save();
      for (let c = 0; c < baits.length; c++) {
        const x = leftPad + c * cellW + cellW / 2;
        const y = topPad - 12;
        const isHoveredCol = hoveredHit && hoveredHit.colIdx === c;
        ctx.save();
        ctx.translate(x, y);
        ctx.rotate(-Math.PI / 4.5);
        ctx.font = isHoveredCol
          ? "700 12.5px Inter, system-ui, sans-serif"
          : "600 12px Inter, system-ui, sans-serif";
        ctx.fillStyle = isHoveredCol ? "#1d4ed8" : "#0f172a";
        ctx.textAlign = "left";
        ctx.textBaseline = "middle";
        ctx.fillText(baits[c], 0, 0);
        ctx.restore();
      }
      ctx.restore();

      // Row labels (Preys + shared complex dot indicator)
      ctx.textAlign = "right";
      ctx.textBaseline = "middle";
      for (let r = 0; r < preys.length; r++) {
        const y = topPad + r * cellH + cellH / 2;
        const isHoveredRow = hoveredHit && hoveredHit.rowIdx === r;
        const p = preys[r];
        ctx.font =
          isHoveredRow || p.isShared
            ? "600 12px 'JetBrains Mono', Inter, monospace"
            : "500 12px 'JetBrains Mono', Inter, monospace";
        ctx.fillStyle = isHoveredRow ? "#1d4ed8" : p.isShared ? "#0f172a" : "#334155";
        ctx.fillText(p.gene_name.slice(0, 16), leftPad - (p.isShared ? 18 : 10), y);

        if (p.isShared) {
          ctx.beginPath();
          ctx.arc(leftPad - 9, y, 3.2, 0, Math.PI * 2);
          ctx.fillStyle = "#2563eb";
          ctx.fill();
        }
      }

      // Draw interaction dots
      for (const h of hitRegions) {
        const cell = h.row;
        const rad = h.rad;
        const fillCol = interpolateColor(cell.SaintScore, colormap);
        const isSelected =
          selectedRow && selectedRow.Bait === cell.Bait && selectedRow.Prey === cell.Prey;
        const isHovered = hoveredHit === h;

        if (isSelected || isHovered) {
          ctx.beginPath();
          ctx.arc(h.cx, h.cy, rad + 4.5, 0, Math.PI * 2);
          ctx.fillStyle = isSelected ? "rgba(37, 99, 235, 0.18)" : "rgba(15, 23, 42, 0.12)";
          ctx.fill();
        }

        ctx.beginPath();
        ctx.arc(h.cx, h.cy, rad, 0, Math.PI * 2);
        ctx.fillStyle = fillCol;
        ctx.fill();

        // BFDR Significance Border Ring
        const bfdr = cell.computedBFDR;
        ctx.beginPath();
        ctx.arc(h.cx, h.cy, rad + 0.8, 0, Math.PI * 2);
        if (!Number.isNaN(bfdr) && bfdr <= 0.01) {
          ctx.strokeStyle = "#0f172a";
          ctx.lineWidth = 2.1;
          ctx.setLineDash([]);
          ctx.stroke();
        } else if (!Number.isNaN(bfdr) && bfdr <= 0.05) {
          ctx.strokeStyle = "#334155";
          ctx.lineWidth = 1.6;
          ctx.setLineDash([2.5, 2]);
          ctx.stroke();
          ctx.setLineDash([]);
        } else {
          ctx.strokeStyle = "rgba(100,116,139,0.35)";
          ctx.lineWidth = 0.8;
          ctx.setLineDash([]);
          ctx.stroke();
        }

        if (isSelected) {
          ctx.beginPath();
          ctx.arc(h.cx, h.cy, rad + 3.2, 0, Math.PI * 2);
          ctx.strokeStyle = "#2563eb";
          ctx.lineWidth = 1.8;
          ctx.stroke();
        }
      }

      // --- Right-hand ProHits Legend Card ---
      const legCardX = leftPad + gridW + 20;
      const legCardY = topPad;
      const legCardW = rightLegendW - 32;
      const legCardH = 262;

      drawRoundedRect(ctx, legCardX, legCardY, legCardW, legCardH, 8);
      ctx.fillStyle = "#f8fafc";
      ctx.fill();
      ctx.strokeStyle = "#e2e8f0";
      ctx.lineWidth = 1;
      ctx.stroke();

      const legX = legCardX + 14;
      let legY = legCardY + 14;

      // 1. SaintScore Colorbar
      ctx.fillStyle = "#0f172a";
      ctx.font = "700 11px Inter, system-ui, sans-serif";
      ctx.textAlign = "left";
      ctx.textBaseline = "top";
      ctx.fillText("SaintScore (AvgP)", legX, legY);
      legY += 18;

      const barW = legCardW - 28;
      const barH = 12;
      for (let i = 0; i <= barW; i++) {
        ctx.fillStyle = interpolateColor(i / barW, colormap);
        ctx.fillRect(legX + i, legY, 1.5, barH);
      }
      ctx.strokeStyle = "#cbd5e1";
      ctx.lineWidth = 1;
      ctx.strokeRect(legX, legY, barW, barH);

      ctx.fillStyle = "#475569";
      ctx.font = "600 10px 'JetBrains Mono', Inter, monospace";
      ctx.fillText("0.0", legX, legY + barH + 4);
      ctx.textAlign = "center";
      ctx.fillText("0.5", legX + barW / 2, legY + barH + 4);
      ctx.textAlign = "right";
      ctx.fillText("1.0", legX + barW, legY + barH + 4);

      // 2. Relative Abundance Size Legend
      legY += barH + 26;
      ctx.textAlign = "left";
      ctx.fillStyle = "#0f172a";
      ctx.font = "700 11px Inter, system-ui, sans-serif";
      ctx.fillText("Quantitative Abundance", legX, legY);
      legY += 18;

      const sizeTiers = [0.1, 0.5, 1.0];
      let curX = legX + 16;
      for (const frac of sizeTiers) {
        const val = maxAbundance * frac;
        const rad = calcRadius(val);
        ctx.beginPath();
        ctx.arc(curX, legY + maxRadius, rad, 0, Math.PI * 2);
        ctx.fillStyle = "#cbd5e1";
        ctx.fill();
        ctx.strokeStyle = "#334155";
        ctx.lineWidth = 1.2;
        ctx.stroke();

        ctx.fillStyle = "#475569";
        ctx.font = "600 10px 'JetBrains Mono', Inter, monospace";
        ctx.textAlign = "center";
        const valTxt = val >= 10 ? Math.round(val) : val.toFixed(1);
        ctx.fillText(String(valTxt), curX, legY + maxRadius * 2 + 5);
        curX += 50;
      }

      // 3. BFDR Border Ring Legend
      legY += maxRadius * 2 + 28;
      ctx.textAlign = "left";
      ctx.fillStyle = "#0f172a";
      ctx.font = "700 11px Inter, system-ui, sans-serif";
      ctx.fillText("Significance (BFDR)", legX, legY);
      legY += 18;

      const ringItems = [
        { label: "BFDR ≤ 1% (High)", dash: [], width: 2.1, stroke: "#0f172a" },
        { label: "BFDR ≤ 5% (Mod)", dash: [2.5, 2], width: 1.6, stroke: "#334155" },
        { label: "BFDR > 5% (Low)", dash: [], width: 0.8, stroke: "#94a3b8" },
      ];
      for (const item of ringItems) {
        ctx.beginPath();
        ctx.arc(legX + 8, legY + 6, 6, 0, Math.PI * 2);
        ctx.fillStyle = "#ffffff";
        ctx.fill();
        ctx.strokeStyle = item.stroke;
        ctx.lineWidth = item.width;
        ctx.setLineDash(item.dash);
        ctx.stroke();
        ctx.setLineDash([]);

        ctx.fillStyle = "#334155";
        ctx.font = "500 11px Inter, system-ui, sans-serif";
        ctx.fillText(item.label, legX + 22, legY + 1);
        legY += 20;
      }

      // Shared prey indicator note
      legY += 4;
      ctx.beginPath();
      ctx.arc(legX + 8, legY + 5, 3.2, 0, Math.PI * 2);
      ctx.fillStyle = "#2563eb";
      ctx.fill();
      ctx.fillStyle = "#475569";
      ctx.font = "500 10.5px Inter, system-ui, sans-serif";
      ctx.fillText("Multi-bait shared prey", legX + 18, legY);
    }

    drawScene();

    // Custom hover + click interactivity with crosshair redraw
    canvas.onmousemove = (ev) => {
      const rect = canvas.getBoundingClientRect();
      const mx = ev.clientX - rect.left;
      const my = ev.clientY - rect.top;
      let hit = null;
      for (let i = hitRegions.length - 1; i >= 0; i--) {
        const h = hitRegions[i];
        const dx = mx - h.cx;
        const dy = my - h.cy;
        if (dx * dx + dy * dy <= h.r * h.r) {
          hit = h;
          break;
        }
      }
      if (hit !== hoveredHit) {
        hoveredHit = hit;
        drawScene();
      }
      canvas.style.cursor = hit ? "pointer" : "default";
      if (hit && callbacks.onHoverRow) {
        callbacks.onHoverRow(hit.row, ev.clientX, ev.clientY);
      } else if (!hit && callbacks.onLeave) {
        callbacks.onLeave();
      }
    };

    canvas.onclick = (ev) => {
      const rect = canvas.getBoundingClientRect();
      const mx = ev.clientX - rect.left;
      const my = ev.clientY - rect.top;
      for (let i = hitRegions.length - 1; i >= 0; i--) {
        const h = hitRegions[i];
        const dx = mx - h.cx;
        const dy = my - h.cy;
        if (dx * dx + dy * dy <= h.r * h.r) {
          if (callbacks.onSelectRow) callbacks.onSelectRow(h.row);
          break;
        }
      }
    };

    canvas.onmouseleave = () => {
      if (hoveredHit) {
        hoveredHit = null;
        drawScene();
      }
      canvas.style.cursor = "default";
      if (callbacks.onLeave) callbacks.onLeave();
    };

    // Complete Publication Vector SVG Generator (including Right-Hand Legend!)
    function toSVG() {
      const gridW = baits.length * cellW;
      const gridH = preys.length * cellH;
      const legCardX = leftPad + gridW + 20;
      const legCardY = topPad;
      const legCardW = rightLegendW - 32;
      const legX = legCardX + 14;
      const barW = legCardW - 28;

      const svgParts = [
        `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" font-family="Inter, -apple-system, BlinkMacSystemFont, sans-serif">`,
        `<defs>`,
        `  <linearGradient id="dotplot-cbar" x1="0%" y1="0%" x2="100%" y2="0%">`,
      ];
      for (let s = 0; s <= 10; s++) {
        const pct = s * 10;
        svgParts.push(`    <stop offset="${pct}%" stop-color="${interpolateColor(s / 10, colormap)}"/>`);
      }
      svgParts.push(
        `  </linearGradient>`,
        `</defs>`,
        `<rect width="${width}" height="${height}" fill="#ffffff"/>`,
        `<text x="16" y="26" font-size="13" font-weight="700" fill="#0f172a">ProHits Bait × Prey Interaction Matrix</text>`,
        `<text x="16" y="43" font-size="11" font-weight="500" fill="#64748b">${escapeXml(
          `${baits.length} Baits × ${preys.length} Top Preys · Color = SaintScore · Size = Abundance`
        )}</text>`
      );

      for (let r = 0; r < preys.length; r++) {
        const y = topPad + r * cellH;
        if (r % 2 === 0) {
          svgParts.push(
            `<rect x="${leftPad}" y="${y}" width="${gridW}" height="${cellH}" fill="#f8fafc"/>`
          );
        }
        svgParts.push(
          `<line x1="${leftPad}" y1="${y + cellH / 2}" x2="${leftPad + gridW}" y2="${
            y + cellH / 2
          }" stroke="#e2e8f0" stroke-width="1"/>`
        );
      }
      for (let c = 0; c < baits.length; c++) {
        const x = leftPad + c * cellW + cellW / 2;
        svgParts.push(
          `<line x1="${x}" y1="${topPad}" x2="${x}" y2="${
            topPad + gridH
          }" stroke="#e2e8f0" stroke-width="1"/>`,
          `<text x="${x}" y="${topPad - 12}" transform="rotate(-40 ${x} ${
            topPad - 12
          })" font-size="12" font-weight="600" fill="#0f172a">${escapeXml(baits[c])}</text>`
        );
      }
      svgParts.push(
        `<rect x="${leftPad}" y="${topPad}" width="${gridW}" height="${gridH}" fill="none" stroke="#cbd5e1" stroke-width="1.4"/>`
      );
      for (let r = 0; r < preys.length; r++) {
        const y = topPad + r * cellH + cellH / 2 + 4;
        const p = preys[r];
        svgParts.push(
          `<text x="${leftPad - (p.isShared ? 18 : 10)}" y="${y}" text-anchor="end" font-size="12" font-weight="${
            p.isShared ? 600 : 500
          }" font-family="JetBrains Mono, monospace" fill="#1e293b">${escapeXml(p.gene_name)}</text>`
        );
        if (p.isShared) {
          svgParts.push(`<circle cx="${leftPad - 9}" cy="${y - 4}" r="3.2" fill="#2563eb"/>`);
        }
      }
      for (const h of hitRegions) {
        const cell = h.row;
        const rad = h.rad;
        const fillCol = interpolateColor(cell.SaintScore, colormap);
        const bfdr = cell.computedBFDR;
        let strokeAttr = `stroke="rgba(100,116,139,0.35)" stroke-width="0.8"`;
        if (!Number.isNaN(bfdr) && bfdr <= 0.01) {
          strokeAttr = `stroke="#0f172a" stroke-width="2.1"`;
        } else if (!Number.isNaN(bfdr) && bfdr <= 0.05) {
          strokeAttr = `stroke="#334155" stroke-width="1.6" stroke-dasharray="2.5,2"`;
        }
        svgParts.push(
          `<circle cx="${h.cx.toFixed(1)}" cy="${h.cy.toFixed(1)}" r="${rad.toFixed(
            2
          )}" fill="${fillCol}" ${strokeAttr}><title>${escapeXml(
            `${cell.Bait} → ${cell.PreyGene}: SaintScore=${cell.SaintScore.toFixed(3)}, BFDR=${
              Number.isNaN(bfdr) ? "N/A" : bfdr.toFixed(3)
            }, AvgSpec=${cell.AvgSpec.toFixed(2)}`
          )}</title></circle>`
        );
      }

      // Legend SVG elements
      let ly = legCardY + 14;
      svgParts.push(
        `<rect x="${legCardX}" y="${legCardY}" width="${legCardW}" height="262" rx="8" fill="#f8fafc" stroke="#e2e8f0" stroke-width="1"/>`,
        `<text x="${legX}" y="${ly + 10}" font-size="11" font-weight="700" fill="#0f172a">SaintScore (AvgP)</text>`
      );
      ly += 18;
      svgParts.push(
        `<rect x="${legX}" y="${ly}" width="${barW}" height="12" fill="url(#dotplot-cbar)" stroke="#cbd5e1" stroke-width="1"/>`,
        `<text x="${legX}" y="${ly + 24}" font-size="10" fill="#475569">0.0</text>`,
        `<text x="${legX + barW / 2}" y="${ly + 24}" text-anchor="middle" font-size="10" fill="#475569">0.5</text>`,
        `<text x="${legX + barW}" y="${ly + 24}" text-anchor="end" font-size="10" fill="#475569">1.0</text>`
      );
      ly += 38;
      svgParts.push(
        `<text x="${legX}" y="${ly + 10}" font-size="11" font-weight="700" fill="#0f172a">Quantitative Abundance</text>`
      );
      ly += 18;
      let cxLeg = legX + 16;
      for (const frac of [0.1, 0.5, 1.0]) {
        const val = maxAbundance * frac;
        const rad = calcRadius(val);
        const valTxt = val >= 10 ? Math.round(val) : val.toFixed(1);
        svgParts.push(
          `<circle cx="${cxLeg}" cy="${ly + maxRadius}" r="${rad.toFixed(
            1
          )}" fill="#cbd5e1" stroke="#334155" stroke-width="1.2"/>`,
          `<text x="${cxLeg}" y="${
            ly + maxRadius * 2 + 14
          }" text-anchor="middle" font-size="10" fill="#475569">${valTxt}</text>`
        );
        cxLeg += 50;
      }
      ly += maxRadius * 2 + 28;
      svgParts.push(
        `<text x="${legX}" y="${ly + 10}" font-size="11" font-weight="700" fill="#0f172a">Significance (BFDR)</text>`
      );
      ly += 18;
      const svgRings = [
        { label: "BFDR ≤ 1% (High)", attr: `stroke="#0f172a" stroke-width="2.1"` },
        { label: "BFDR ≤ 5% (Mod)", attr: `stroke="#334155" stroke-width="1.6" stroke-dasharray="2.5,2"` },
        { label: "BFDR > 5% (Low)", attr: `stroke="#94a3b8" stroke-width="0.8"` },
      ];
      for (const sr of svgRings) {
        svgParts.push(
          `<circle cx="${legX + 8}" cy="${ly + 6}" r="6" fill="#ffffff" ${sr.attr}/>`,
          `<text x="${legX + 22}" y="${ly + 10}" font-size="11" fill="#334155">${escapeXml(sr.label)}</text>`
        );
        ly += 20;
      }
      svgParts.push(`</svg>`);
      return svgParts.join("\n");
    }

    return { matrix, hitRegions, toSVG };
  }

  // ============================================================================
  // 2. INTERACTIVE VOLCANO & SIGNIFICANCE SCATTER PLOT
  // ============================================================================

  function renderVolcanoPlot(canvas, rows, options = {}, callbacks = {}) {
    const {
      xAxisMode = "log2FC",
      yAxisMode = "SaintScore",
      scoreThreshold = 0.8,
      bfdrThreshold = 0.05,
      fcThreshold = 2.0,
      baitFilter = "ALL",
      maxLabels = 20,
      selectedRow = null,
    } = options;

    const width = Math.max(600, canvas.parentElement ? canvas.parentElement.clientWidth - 24 : 700);
    const height = 450;
    const ctx = setupHiDPICanvas(canvas, width, height);

    const padL = 68;
    const padR = 32;
    const padT = 50;
    const padB = 56;
    const plotW = width - padL - padR;
    const plotH = height - padT - padB;

    const filtered = baitFilter === "ALL" ? rows : rows.filter((r) => r.Bait === baitFilter);
    if (filtered.length === 0) {
      ctx.fillStyle = "#ffffff";
      ctx.fillRect(0, 0, width, height);
      ctx.fillStyle = "#64748b";
      ctx.font = "500 14px Inter, system-ui, sans-serif";
      ctx.textAlign = "center";
      ctx.fillText("No interactions available for Volcano Plot.", width / 2, height / 2);
      return { hitRegions: [], toSVG: () => "" };
    }

    const uniqueBaits = Array.from(new Set(rows.map((r) => r.Bait)));
    const baitColorMap = new Map(
      uniqueBaits.map((b, idx) => [b, BAIT_PALETTE[idx % BAIT_PALETTE.length]])
    );

    function getX(r) {
      if (xAxisMode === "logOddsScore") {
        return Math.max(-20, Math.min(50, Number.isNaN(r.logOddsScore) ? 0 : r.logOddsScore));
      }
      if (xAxisMode === "log10Spec") {
        return Math.log10(Math.max(0, r.AvgSpec) + 1);
      }
      return Number.isNaN(r.log2FC) ? 0 : Math.max(-6, Math.min(12, r.log2FC));
    }

    function getY(r) {
      if (yAxisMode === "negLog10BFDR") {
        const b = Number.isNaN(r.computedBFDR) ? 1.0 : Math.max(1e-4, r.computedBFDR);
        return -Math.log10(b);
      }
      if (yAxisMode === "MaxP") {
        return Number.isNaN(r.MaxP) ? 0 : r.MaxP;
      }
      return Number.isNaN(r.SaintScore) ? 0 : r.SaintScore;
    }

    const pts = filtered.map((r) => {
      const x = getX(r);
      const y = getY(r);
      const passScore = r.SaintScore >= scoreThreshold;
      const passBfdr = Number.isNaN(r.computedBFDR) || r.computedBFDR <= bfdrThreshold;
      const passFc = xAxisMode === "log2FC" ? r.log2FC >= Math.log2(Math.max(1, fcThreshold)) : true;
      const isSig = passScore && passBfdr && passFc;
      return { row: r, x, y, isSig };
    });

    let minX = Infinity,
      maxX = -Infinity,
      minY = 0,
      maxY = -Infinity;
    for (const p of pts) {
      if (p.x < minX) minX = p.x;
      if (p.x > maxX) maxX = p.x;
      if (p.y > maxY) maxY = p.y;
    }
    if (xAxisMode === "log2FC") {
      minX = Math.min(-1.5, minX - 0.6);
      maxX = Math.max(4.0, maxX + 0.9);
    } else {
      const span = Math.max(1, maxX - minX);
      minX -= span * 0.08;
      maxX += span * 0.08;
    }
    if (yAxisMode === "SaintScore" || yAxisMode === "MaxP") {
      minY = -0.04;
      maxY = 1.07;
    } else {
      minY = -0.1;
      maxY = Math.max(2.5, maxY * 1.12);
    }

    const toPxX = (x) => padL + ((x - minX) / (maxX - minX || 1)) * plotW;
    const toPxY = (y) => padT + plotH - ((y - minY) / (maxY - minY || 1)) * plotH;

    const xCutVal = xAxisMode === "log2FC" ? Math.log2(Math.max(1, fcThreshold)) : 0;
    const yCutVal =
      yAxisMode === "negLog10BFDR"
        ? -Math.log10(Math.max(1e-4, bfdrThreshold))
        : scoreThreshold;

    const xCutPx = Math.max(padL, Math.min(padL + plotW, toPxX(xCutVal)));
    const yCutPx = Math.max(padT, Math.min(padT + plotH, toPxY(yCutVal)));

    const xLabels = {
      log2FC: "Log₂ Fold Change (Regularized Test vs. Control)",
      logOddsScore: "SAINT Log-Odds Score (min across replicates)",
      log10Spec: "Log₁₀(Avg Abundance + 1)",
    };
    const yLabels = {
      SaintScore: "SaintScore (AvgP Posterior Probability)",
      negLog10BFDR: "-Log₁₀(Bayesian FDR)",
      MaxP: "MaxP (Maximum Replicate Probability)",
    };

    const hitRegions = [];
    const sigPts = [];
    for (const p of pts) {
      p.cx = toPxX(p.x);
      p.cy = toPxY(p.y);
      if (p.isSig) sigPts.push(p);
    }

    // Pre-compute smart non-overlapping gene callout labels
    sigPts.sort((a, b) => b.row.SaintScore - a.row.SaintScore || b.x - a.x);
    ctx.font = "600 11px Inter, system-ui, sans-serif";
    const placedLabels = [];
    for (const p of sigPts) {
      if (placedLabels.length >= maxLabels) break;
      const txt =
        uniqueBaits.length > 1 && baitFilter === "ALL" && uniqueBaits.length <= 6
          ? `${p.row.PreyGene} (${p.row.Bait})`
          : p.row.PreyGene;
      const tw = ctx.measureText(txt).width + 10;
      const th = 16;
      const candidates = [
        { lx: p.cx + 9, ly: p.cy - 9, align: "left" },
        { lx: p.cx - 9, ly: p.cy - 9, align: "right" },
        { lx: p.cx + 9, ly: p.cy + 11, align: "left" },
        { lx: p.cx - 9, ly: p.cy + 11, align: "right" },
        { lx: p.cx, ly: p.cy - 14, align: "center" },
      ];
      for (const cand of candidates) {
        const boxX =
          cand.align === "left"
            ? cand.lx
            : cand.align === "right"
            ? cand.lx - tw
            : cand.lx - tw / 2;
        const boxY = cand.ly - th / 2;
        if (
          boxX < padL + 3 ||
          boxX + tw > padL + plotW - 3 ||
          boxY < padT + 3 ||
          boxY + th > padT + plotH - 3
        ) {
          continue;
        }
        const overlaps = placedLabels.some(
          (b) => !(boxX + tw + 2 < b.x || boxX - 2 > b.x + b.w || boxY + th + 2 < b.y || boxY - 2 > b.y + b.h)
        );
        if (!overlaps) {
          placedLabels.push({
            x: boxX,
            y: boxY,
            w: tw,
            h: th,
            txt,
            pt: p,
            color: baitColorMap.get(p.row.Bait) || "#2563eb",
          });
          break;
        }
      }
    }

    for (const p of pts) {
      if (!p.isSig) {
        hitRegions.push({ cx: p.cx, cy: p.cy, r: 6, row: p.row, isSig: false });
      }
    }
    for (const p of sigPts) {
      hitRegions.push({ cx: p.cx, cy: p.cy, r: 8, row: p.row, isSig: true });
    }

    let hoveredPt = null;

    function drawVolcano() {
      ctx.fillStyle = "#ffffff";
      ctx.fillRect(0, 0, width, height);

      // Highlight upper-right significant enrichment zone
      if (xCutPx < padL + plotW && yCutPx > padT) {
        ctx.fillStyle = "rgba(16, 185, 129, 0.06)";
        ctx.fillRect(xCutPx, padT, padL + plotW - xCutPx, yCutPx - padT);

        // Subtle enrichment zone badge
        ctx.fillStyle = "#047857";
        ctx.font = "600 10.5px Inter, system-ui, sans-serif";
        ctx.textAlign = "right";
        ctx.textBaseline = "top";
        ctx.fillText("HIGH-CONFIDENCE ENRICHMENT ZONE", padL + plotW - 10, padT + 8);
      }

      // Gridlines & ticks
      ctx.strokeStyle = "#f1f5f9";
      ctx.lineWidth = 1;
      const nTicksX = 6;
      const nTicksY = 5;
      ctx.font = "500 11px 'JetBrains Mono', Inter, monospace";
      ctx.fillStyle = "#64748b";

      for (let i = 0; i <= nTicksX; i++) {
        const val = minX + (i / nTicksX) * (maxX - minX);
        const px = toPxX(val);
        ctx.beginPath();
        ctx.moveTo(px, padT);
        ctx.lineTo(px, padT + plotH);
        ctx.stroke();
        ctx.textAlign = "center";
        ctx.textBaseline = "top";
        ctx.fillText(val.toFixed(1), px, padT + plotH + 8);
      }

      for (let i = 0; i <= nTicksY; i++) {
        const val = minY + (i / nTicksY) * (maxY - minY);
        if (val < 0 || ((yAxisMode === "SaintScore" || yAxisMode === "MaxP") && val > 1.001)) continue;
        const py = toPxY(val);
        ctx.beginPath();
        ctx.moveTo(padL, py);
        ctx.lineTo(padL + plotW, py);
        ctx.stroke();
        ctx.textAlign = "right";
        ctx.textBaseline = "middle";
        ctx.fillText(val.toFixed(2), padL - 8, py);
      }

      // Threshold dashed lines
      ctx.save();
      ctx.setLineDash([5, 4]);
      ctx.strokeStyle = "#94a3b8";
      ctx.lineWidth = 1.4;
      if (xCutVal >= minX && xCutVal <= maxX) {
        ctx.beginPath();
        ctx.moveTo(xCutPx, padT);
        ctx.lineTo(xCutPx, padT + plotH);
        ctx.stroke();
      }
      if (yCutVal >= minY && yCutVal <= maxY) {
        ctx.beginPath();
        ctx.moveTo(padL, yCutPx);
        ctx.lineTo(padL + plotW, yCutPx);
        ctx.stroke();
      }
      ctx.restore();

      // Plot frame
      ctx.strokeStyle = "#cbd5e1";
      ctx.lineWidth = 1.3;
      ctx.strokeRect(padL, padT, plotW, plotH);

      // Axis titles
      ctx.fillStyle = "#0f172a";
      ctx.font = "600 12px Inter, system-ui, sans-serif";
      ctx.textAlign = "center";
      ctx.textBaseline = "bottom";
      ctx.fillText(xLabels[xAxisMode] || xAxisMode, padL + plotW / 2, height - 10);

      ctx.save();
      ctx.translate(18, padT + plotH / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.textAlign = "center";
      ctx.textBaseline = "top";
      ctx.fillText(yLabels[yAxisMode] || yAxisMode, 0, 0);
      ctx.restore();

      // Background non-significant points
      for (const p of pts) {
        if (!p.isSig) {
          ctx.beginPath();
          ctx.arc(p.cx, p.cy, 4.0, 0, Math.PI * 2);
          ctx.fillStyle = "rgba(148, 163, 184, 0.42)";
          ctx.fill();
        }
      }

      // Significant points
      for (const p of sigPts) {
        const col = baitColorMap.get(p.row.Bait) || "#2563eb";
        const isSel =
          selectedRow && selectedRow.Bait === p.row.Bait && selectedRow.Prey === p.row.Prey;
        const isHov = hoveredPt && hoveredPt.row === p.row;

        if (isSel || isHov) {
          ctx.beginPath();
          ctx.arc(p.cx, p.cy, 10, 0, Math.PI * 2);
          ctx.fillStyle = "rgba(37, 99, 235, 0.18)";
          ctx.fill();
        }

        ctx.beginPath();
        ctx.arc(p.cx, p.cy, isSel || isHov ? 6.8 : 5.6, 0, Math.PI * 2);
        ctx.fillStyle = col;
        ctx.fill();
        ctx.strokeStyle = "#ffffff";
        ctx.lineWidth = 1.4;
        ctx.stroke();
      }

      // Gene Callout Pills
      ctx.font = "600 10.5px Inter, system-ui, sans-serif";
      ctx.textBaseline = "middle";
      for (const lb of placedLabels) {
        ctx.strokeStyle = "rgba(100,116,139,0.55)";
        ctx.lineWidth = 0.9;
        ctx.beginPath();
        ctx.moveTo(lb.pt.cx, lb.pt.cy);
        ctx.lineTo(lb.x + lb.w / 2, lb.y + lb.h / 2);
        ctx.stroke();

        drawRoundedRect(ctx, lb.x, lb.y, lb.w, lb.h, 4);
        ctx.fillStyle = "rgba(255, 255, 255, 0.92)";
        ctx.fill();
        ctx.strokeStyle = "#cbd5e1";
        ctx.lineWidth = 0.8;
        ctx.stroke();

        ctx.fillStyle = "#0f172a";
        ctx.textAlign = "center";
        ctx.fillText(lb.txt, lb.x + lb.w / 2, lb.y + lb.h / 2 + 0.5);
      }

      // Top Header Summary + Bait Color Legend
      ctx.fillStyle = "#0f172a";
      ctx.font = "700 12.5px Inter, system-ui, sans-serif";
      ctx.textAlign = "left";
      ctx.textBaseline = "middle";
      ctx.fillText(`Significant Hits: ${sigPts.length} / ${pts.length} pairs`, padL, 22);

      // Draw Bait Legend pills on top right if <= 8 baits
      const shownBaits =
        baitFilter === "ALL" ? uniqueBaits.slice(0, 6) : [baitFilter];
      let bx = padL + 230;
      ctx.font = "600 11px Inter, system-ui, sans-serif";
      for (const b of shownBaits) {
        if (bx + 80 > width - padR) break;
        const c = baitColorMap.get(b) || "#2563eb";
        ctx.beginPath();
        ctx.arc(bx, 22, 4.5, 0, Math.PI * 2);
        ctx.fillStyle = c;
        ctx.fill();
        ctx.fillStyle = "#334155";
        ctx.fillText(b, bx + 8, 22);
        bx += ctx.measureText(b).width + 24;
      }
    }

    drawVolcano();

    canvas.onmousemove = (ev) => {
      const rect = canvas.getBoundingClientRect();
      const mx = ev.clientX - rect.left;
      const my = ev.clientY - rect.top;
      let hit = null;
      for (let i = hitRegions.length - 1; i >= 0; i--) {
        const h = hitRegions[i];
        const dx = mx - h.cx;
        const dy = my - h.cy;
        if (dx * dx + dy * dy <= h.r * h.r) {
          hit = h;
          break;
        }
      }
      if (hit !== hoveredPt) {
        hoveredPt = hit;
        drawVolcano();
      }
      canvas.style.cursor = hit ? "pointer" : "default";
      if (hit && callbacks.onHoverRow) {
        callbacks.onHoverRow(hit.row, ev.clientX, ev.clientY);
      } else if (!hit && callbacks.onLeave) {
        callbacks.onLeave();
      }
    };

    canvas.onclick = (ev) => {
      const rect = canvas.getBoundingClientRect();
      const mx = ev.clientX - rect.left;
      const my = ev.clientY - rect.top;
      for (let i = hitRegions.length - 1; i >= 0; i--) {
        const h = hitRegions[i];
        const dx = mx - h.cx;
        const dy = my - h.cy;
        if (dx * dx + dy * dy <= h.r * h.r) {
          if (callbacks.onSelectRow) callbacks.onSelectRow(h.row);
          break;
        }
      }
    };

    canvas.onmouseleave = () => {
      if (hoveredPt) {
        hoveredPt = null;
        drawVolcano();
      }
      canvas.style.cursor = "default";
      if (callbacks.onLeave) callbacks.onLeave();
    };

    // Complete Publication Vector SVG Generator (with axes, ticks, gene pills & legend!)
    function toSVG() {
      const svg = [
        `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" font-family="Inter, -apple-system, sans-serif">`,
        `<rect width="${width}" height="${height}" fill="#ffffff"/>`,
        `<text x="${padL}" y="26" font-size="12.5" font-weight="700" fill="#0f172a">Significant Hits: ${sigPts.length} / ${pts.length} pairs</text>`,
      ];

      if (xCutPx < padL + plotW && yCutPx > padT) {
        svg.push(
          `<rect x="${xCutPx}" y="${padT}" width="${padL + plotW - xCutPx}" height="${
            yCutPx - padT
          }" fill="rgba(16,185,129,0.06)"/>`,
          `<text x="${padL + plotW - 10}" y="${
            padT + 18
          }" text-anchor="end" font-size="10.5" font-weight="600" fill="#047857">HIGH-CONFIDENCE ENRICHMENT ZONE</text>`
        );
      }

      // Ticks & gridlines
      for (let i = 0; i <= 6; i++) {
        const val = minX + (i / 6) * (maxX - minX);
        const px = toPxX(val);
        svg.push(
          `<line x1="${px.toFixed(1)}" y1="${padT}" x2="${px.toFixed(1)}" y2="${
            padT + plotH
          }" stroke="#f1f5f9"/>`,
          `<text x="${px.toFixed(1)}" y="${
            padT + plotH + 18
          }" text-anchor="middle" font-size="11" fill="#64748b">${val.toFixed(1)}</text>`
        );
      }
      for (let i = 0; i <= 5; i++) {
        const val = minY + (i / 5) * (maxY - minY);
        if (val < 0 || ((yAxisMode === "SaintScore" || yAxisMode === "MaxP") && val > 1.001)) continue;
        const py = toPxY(val);
        svg.push(
          `<line x1="${padL}" y1="${py.toFixed(1)}" x2="${padL + plotW}" y2="${py.toFixed(
            1
          )}" stroke="#f1f5f9"/>`,
          `<text x="${padL - 8}" y="${(py + 4).toFixed(
            1
          )}" text-anchor="end" font-size="11" fill="#64748b">${val.toFixed(2)}</text>`
        );
      }

      svg.push(
        `<line x1="${xCutPx.toFixed(1)}" y1="${padT}" x2="${xCutPx.toFixed(1)}" y2="${
          padT + plotH
        }" stroke="#94a3b8" stroke-width="1.4" stroke-dasharray="5,4"/>`,
        `<line x1="${padL}" y1="${yCutPx.toFixed(1)}" x2="${padL + plotW}" y2="${yCutPx.toFixed(
          1
        )}" stroke="#94a3b8" stroke-width="1.4" stroke-dasharray="5,4"/>`,
        `<rect x="${padL}" y="${padT}" width="${plotW}" height="${plotH}" fill="none" stroke="#cbd5e1" stroke-width="1.3"/>`
      );

      for (const p of pts) {
        if (!p.isSig) {
          svg.push(
            `<circle cx="${p.cx.toFixed(1)}" cy="${p.cy.toFixed(
              1
            )}" r="4.0" fill="rgba(148,163,184,0.42)"/>`
          );
        }
      }
      for (const p of sigPts) {
        const col = baitColorMap.get(p.row.Bait) || "#2563eb";
        svg.push(
          `<circle cx="${p.cx.toFixed(1)}" cy="${p.cy.toFixed(
            1
          )}" r="5.6" fill="${col}" stroke="#ffffff" stroke-width="1.4"/>`
        );
      }

      for (const lb of placedLabels) {
        svg.push(
          `<line x1="${lb.pt.cx.toFixed(1)}" y1="${lb.pt.cy.toFixed(1)}" x2="${(
            lb.x + lb.w / 2
          ).toFixed(1)}" y2="${(lb.y + lb.h / 2).toFixed(
            1
          )}" stroke="rgba(100,116,139,0.55)" stroke-width="0.9"/>`,
          `<rect x="${lb.x.toFixed(1)}" y="${lb.y.toFixed(1)}" width="${lb.w.toFixed(
            1
          )}" height="${lb.h.toFixed(
            1
          )}" rx="4" fill="rgba(255,255,255,0.92)" stroke="#cbd5e1" stroke-width="0.8"/>`,
          `<text x="${(lb.x + lb.w / 2).toFixed(1)}" y="${(lb.y + lb.h / 2 + 3.5).toFixed(
            1
          )}" text-anchor="middle" font-size="10.5" font-weight="600" fill="#0f172a">${escapeXml(
            lb.txt
          )}</text>`
        );
      }

      svg.push(
        `<text x="${padL + plotW / 2}" y="${
          height - 12
        }" text-anchor="middle" font-size="12" font-weight="600" fill="#0f172a">${escapeXml(
          xLabels[xAxisMode] || xAxisMode
        )}</text>`,
        `<text x="18" y="${padT + plotH / 2}" transform="rotate(-90 18 ${
          padT + plotH / 2
        })" text-anchor="middle" font-size="12" font-weight="600" fill="#0f172a">${escapeXml(
          yLabels[yAxisMode] || yAxisMode
        )}</text>`,
        `</svg>`
      );
      return svg.join("\n");
    }

    return { hitRegions, toSVG };
  }

  // ============================================================================
  // 3. FORCE-DIRECTED BAIT-PREY INTERACTOME NETWORK GRAPH
  // ============================================================================

  function renderNetworkGraph(canvas, rows, options = {}, callbacks = {}) {
    const {
      minScore = 0.8,
      maxBfdr = 0.05,
      maxEdges = 120,
      baitFilter = "ALL",
      onlySharedPreys = false,
      selectedRow = null,
    } = options;

    const width = Math.max(600, canvas.parentElement ? canvas.parentElement.clientWidth - 24 : 700);
    const height = 470;
    const ctx = setupHiDPICanvas(canvas, width, height);

    let edges = rows.filter(
      (r) =>
        (baitFilter === "ALL" || r.Bait === baitFilter) &&
        r.SaintScore >= minScore &&
        (Number.isNaN(r.computedBFDR) || r.computedBFDR <= maxBfdr)
    );

    if (edges.length === 0 && rows.length > 0) {
      edges = [...rows]
        .filter((r) => baitFilter === "ALL" || r.Bait === baitFilter)
        .sort((a, b) => b.SaintScore - a.SaintScore)
        .slice(0, 35);
    }

    const preyBaitCount = new Map();
    for (const e of edges) {
      if (!preyBaitCount.has(e.Prey)) preyBaitCount.set(e.Prey, new Set());
      preyBaitCount.get(e.Prey).add(e.Bait);
    }

    if (onlySharedPreys) {
      const shared = edges.filter((e) => (preyBaitCount.get(e.Prey)?.size || 0) >= 2);
      if (shared.length > 0) edges = shared;
    }

    edges.sort((a, b) => b.SaintScore - a.SaintScore || b.AvgSpec - a.AvgSpec);
    if (edges.length > maxEdges) {
      edges = edges.slice(0, maxEdges);
    }

    const nodesMap = new Map();
    const uniqueBaits = Array.from(new Set(rows.map((r) => r.Bait)));
    const baitColorMap = new Map(
      uniqueBaits.map((b, idx) => [b, BAIT_PALETTE[idx % BAIT_PALETTE.length]])
    );

    const activeBaits = Array.from(new Set(edges.map((e) => e.Bait)));
    activeBaits.forEach((b, idx) => {
      const angle = (2 * Math.PI * idx) / Math.max(1, activeBaits.length) - Math.PI / 2;
      const radius = activeBaits.length === 1 ? 0 : Math.min(width, height) * 0.27;
      nodesMap.set(`BAIT:${b}`, {
        id: `BAIT:${b}`,
        label: b,
        type: "bait",
        color: baitColorMap.get(b) || "#2563eb",
        x: width / 2 + Math.cos(angle) * radius,
        y: height / 2 + 12 + Math.sin(angle) * radius,
        vx: 0,
        vy: 0,
        degree: 0,
        r: 18,
      });
    });

    edges.forEach((e, idx) => {
      const bNode = nodesMap.get(`BAIT:${e.Bait}`);
      bNode.degree++;
      const pKey = `PREY:${e.Prey}`;
      if (!nodesMap.has(pKey)) {
        const sharedCount = preyBaitCount.get(e.Prey)?.size || 1;
        const jitterAngle = (idx * 2.39996) % (2 * Math.PI);
        const dist = sharedCount >= 2 ? 42 : 88 + (idx % 38);
        nodesMap.set(pKey, {
          id: pKey,
          prey_id: e.Prey,
          label: e.PreyGene || e.Prey,
          type: "prey",
          sharedCount,
          primaryBaitColor: bNode.color,
          x: (sharedCount >= 2 ? width / 2 : bNode.x) + Math.cos(jitterAngle) * dist,
          y: (sharedCount >= 2 ? height / 2 + 12 : bNode.y) + Math.sin(jitterAngle) * dist,
          vx: 0,
          vy: 0,
          degree: 1,
          maxSpec: e.AvgSpec,
          bestRow: e,
          r: sharedCount >= 2 ? 12 : 9,
        });
      } else {
        const pNode = nodesMap.get(pKey);
        pNode.degree++;
        if (e.AvgSpec > pNode.maxSpec) pNode.maxSpec = e.AvgSpec;
        if (e.SaintScore > pNode.bestRow.SaintScore) pNode.bestRow = e;
      }
    });

    const nodes = Array.from(nodesMap.values());
    const graphLinks = edges.map((e) => ({
      source: nodesMap.get(`BAIT:${e.Bait}`),
      target: nodesMap.get(`PREY:${e.Prey}`),
      row: e,
    }));

    // Deterministic force-directed layout (150 iterations)
    for (let iter = 0; iter < 150; iter++) {
      const alpha = 1 - iter / 160;
      for (let i = 0; i < nodes.length; i++) {
        for (let j = i + 1; j < nodes.length; j++) {
          const n1 = nodes[i];
          const n2 = nodes[j];
          let dx = n2.x - n1.x;
          let dy = n2.y - n1.y;
          let d2 = dx * dx + dy * dy;
          if (d2 < 1) {
            dx = 0.5;
            dy = 0.5;
            d2 = 0.5;
          }
          const dist = Math.sqrt(d2);
          const minDist = n1.r + n2.r + 20;
          const repulse = ((n1.type === "bait" && n2.type === "bait" ? 3200 : 920) / d2) * alpha;
          const fx = (dx / dist) * repulse;
          const fy = (dy / dist) * repulse;
          n1.vx -= fx;
          n1.vy -= fy;
          n2.vx += fx;
          n2.vy += fy;

          if (dist < minDist) {
            const push = (minDist - dist) * 0.28;
            n1.vx -= (dx / dist) * push;
            n1.vy -= (dy / dist) * push;
            n2.vx += (dx / dist) * push;
            n2.vy += (dy / dist) * push;
          }
        }
      }

      for (const link of graphLinks) {
        const s = link.source;
        const t = link.target;
        const dx = t.x - s.x;
        const dy = t.y - s.y;
        const dist = Math.sqrt(dx * dx + dy * dy) || 1;
        const targetLen = t.sharedCount >= 2 ? 110 : 84;
        const force = (dist - targetLen) * 0.045 * alpha;
        const fx = (dx / dist) * force;
        const fy = (dy / dist) * force;
        s.vx += fx * 0.35;
        s.vy += fy * 0.35;
        t.vx -= fx;
        t.vy -= fy;
      }

      for (const n of nodes) {
        n.vx += (width / 2 - n.x) * 0.012 * alpha;
        n.vy += (height / 2 + 12 - n.y) * 0.012 * alpha;
        n.x += Math.max(-12, Math.min(12, n.vx));
        n.y += Math.max(-12, Math.min(12, n.vy));
        n.vx *= 0.65;
        n.vy *= 0.65;
        n.x = Math.max(44, Math.min(width - 44, n.x));
        n.y = Math.max(54, Math.min(height - 34, n.y));
      }
    }

    let hoveredNode = null;
    let draggedNode = null;

    function draw() {
      ctx.fillStyle = "#ffffff";
      ctx.fillRect(0, 0, width, height);

      // Subtle radial background grid for topology feel
      ctx.strokeStyle = "#f8fafc";
      ctx.lineWidth = 1;
      for (const rad of [90, 170, 250]) {
        ctx.beginPath();
        ctx.arc(width / 2, height / 2 + 12, rad, 0, Math.PI * 2);
        ctx.stroke();
      }

      if (nodes.length === 0) {
        ctx.fillStyle = "#64748b";
        ctx.font = "500 14px Inter, system-ui, sans-serif";
        ctx.textAlign = "center";
        ctx.fillText("No interactions meet the network threshold.", width / 2, height / 2);
        return;
      }

      // Draw edges
      for (const link of graphLinks) {
        const isHighlighted =
          !hoveredNode || link.source === hoveredNode || link.target === hoveredNode;
        ctx.beginPath();
        ctx.moveTo(link.source.x, link.source.y);
        ctx.lineTo(link.target.x, link.target.y);
        const w = 1.2 + Math.max(0, link.row.SaintScore - 0.5) * 3.8;
        ctx.lineWidth = isHighlighted ? w : 0.8;
        ctx.strokeStyle = isHighlighted
          ? link.target.sharedCount >= 2
            ? "rgba(217, 119, 6, 0.55)"
            : "rgba(100, 116, 139, 0.45)"
          : "rgba(226, 232, 240, 0.45)";
        ctx.stroke();
      }

      // Draw nodes
      for (const n of nodes) {
        const isHovered = hoveredNode === n;
        const isSelected =
          selectedRow &&
          n.type === "prey" &&
          n.prey_id === selectedRow.Prey;

        if (n.type === "bait") {
          if (isHovered) {
            ctx.beginPath();
            ctx.arc(n.x, n.y, n.r + 6, 0, Math.PI * 2);
            ctx.fillStyle = "rgba(37, 99, 235, 0.15)";
            ctx.fill();
          }
          // Draw Hexagon for Bait Hub
          ctx.beginPath();
          for (let k = 0; k < 6; k++) {
            const ang = (Math.PI / 3) * k - Math.PI / 6;
            const hx = n.x + n.r * Math.cos(ang);
            const hy = n.y + n.r * Math.sin(ang);
            if (k === 0) ctx.moveTo(hx, hy);
            else ctx.lineTo(hx, hy);
          }
          ctx.closePath();
          ctx.fillStyle = n.color;
          ctx.fill();
          ctx.strokeStyle = "#0f172a";
          ctx.lineWidth = 2.2;
          ctx.stroke();

          // Label badge below bait
          const lblW = ctx.measureText(n.label).width + 12;
          drawRoundedRect(ctx, n.x - lblW / 2, n.y + n.r + 3, lblW, 17, 4);
          ctx.fillStyle = "rgba(15, 23, 42, 0.88)";
          ctx.fill();
          ctx.fillStyle = "#ffffff";
          ctx.font = "700 11px Inter, system-ui, sans-serif";
          ctx.textAlign = "center";
          ctx.textBaseline = "middle";
          ctx.fillText(n.label, n.x, n.y + n.r + 11.5);
        } else {
          // Prey node
          if (isHovered || isSelected) {
            ctx.beginPath();
            ctx.arc(n.x, n.y, n.r + 5, 0, Math.PI * 2);
            ctx.fillStyle = isSelected ? "rgba(37, 99, 235, 0.22)" : "rgba(15, 23, 42, 0.12)";
            ctx.fill();
          }

          ctx.beginPath();
          ctx.arc(n.x, n.y, n.r, 0, Math.PI * 2);
          ctx.fillStyle = n.sharedCount >= 2 ? "#fef3c7" : "#ffffff";
          ctx.fill();
          ctx.strokeStyle = n.sharedCount >= 2 ? "#d97706" : n.primaryBaitColor;
          ctx.lineWidth = n.sharedCount >= 2 ? 2.5 : 2.0;
          ctx.stroke();

          ctx.fillStyle = "#0f172a";
          ctx.font = n.sharedCount >= 2 ? "600 10.5px Inter, sans-serif" : "500 10px Inter, sans-serif";
          ctx.textAlign = "center";
          ctx.textBaseline = "bottom";
          ctx.fillText(n.label, n.x, n.y - n.r - 3);
        }
      }

      // Top Header Legend Pill Bar
      drawRoundedRect(ctx, 12, 10, width - 24, 30, 6);
      ctx.fillStyle = "rgba(248, 250, 252, 0.94)";
      ctx.fill();
      ctx.strokeStyle = "#e2e8f0";
      ctx.lineWidth = 1;
      ctx.stroke();

      ctx.font = "600 11.5px Inter, system-ui, sans-serif";
      ctx.fillStyle = "#0f172a";
      ctx.textAlign = "left";
      ctx.textBaseline = "middle";
      ctx.fillText(
        `⬡ Baits: ${activeBaits.length}   ·   ○ Preys: ${nodes.length - activeBaits.length}   ·   ◎ Shared Hub Preys: ${
          nodes.filter((n) => n.sharedCount >= 2).length
        }   ·   Edges: ${graphLinks.length}`,
        24,
        25
      );

      ctx.fillStyle = "#64748b";
      ctx.font = "500 11px Inter, system-ui, sans-serif";
      ctx.textAlign = "right";
      ctx.fillText("Drag any node to rearrange · Click prey to inspect replicates", width - 24, 25);
    }

    draw();

    canvas.onmousedown = (ev) => {
      const rect = canvas.getBoundingClientRect();
      const mx = ev.clientX - rect.left;
      const my = ev.clientY - rect.top;
      for (let i = nodes.length - 1; i >= 0; i--) {
        const n = nodes[i];
        const dx = mx - n.x;
        const dy = my - n.y;
        if (dx * dx + dy * dy <= (n.r + 5) * (n.r + 5)) {
          draggedNode = n;
          if (n.type === "prey" && n.bestRow && callbacks.onSelectRow) {
            callbacks.onSelectRow(n.bestRow);
          }
          break;
        }
      }
    };

    canvas.onmousemove = (ev) => {
      const rect = canvas.getBoundingClientRect();
      const mx = ev.clientX - rect.left;
      const my = ev.clientY - rect.top;
      if (draggedNode) {
        draggedNode.x = Math.max(28, Math.min(width - 28, mx));
        draggedNode.y = Math.max(50, Math.min(height - 28, my));
        draw();
        return;
      }
      let found = null;
      for (let i = nodes.length - 1; i >= 0; i--) {
        const n = nodes[i];
        const dx = mx - n.x;
        const dy = my - n.y;
        if (dx * dx + dy * dy <= (n.r + 5) * (n.r + 5)) {
          found = n;
          break;
        }
      }
      if (found !== hoveredNode) {
        hoveredNode = found;
        canvas.style.cursor = found ? "grab" : "default";
        draw();
      }
      if (found && found.bestRow && callbacks.onHoverRow) {
        callbacks.onHoverRow(found.bestRow, ev.clientX, ev.clientY);
      } else if (!found && callbacks.onLeave) {
        callbacks.onLeave();
      }
    };

    canvas.onmouseup = () => {
      draggedNode = null;
    };
    canvas.onmouseleave = () => {
      draggedNode = null;
      hoveredNode = null;
      draw();
      if (callbacks.onLeave) callbacks.onLeave();
    };

    function toSVG() {
      const svg = [
        `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" font-family="Inter, -apple-system, sans-serif">`,
        `<rect width="${width}" height="${height}" fill="#ffffff"/>`,
        `<rect x="12" y="10" width="${width - 24}" height="30" rx="6" fill="#f8fafc" stroke="#e2e8f0"/>`,
        `<text x="24" y="29" font-size="11.5" font-weight="600" fill="#0f172a">${escapeXml(
          `Baits: ${activeBaits.length} · Preys: ${nodes.length - activeBaits.length} · Shared Hub Preys: ${
            nodes.filter((n) => n.sharedCount >= 2).length
          } · Edges: ${graphLinks.length}`
        )}</text>`,
      ];
      for (const link of graphLinks) {
        const w = (1.2 + Math.max(0, link.row.SaintScore - 0.5) * 3.8).toFixed(2);
        const stroke =
          link.target.sharedCount >= 2 ? "rgba(217,119,6,0.55)" : "rgba(100,116,139,0.45)";
        svg.push(
          `<line x1="${link.source.x.toFixed(1)}" y1="${link.source.y.toFixed(
            1
          )}" x2="${link.target.x.toFixed(1)}" y2="${link.target.y.toFixed(
            1
          )}" stroke="${stroke}" stroke-width="${w}"/>`
        );
      }
      for (const n of nodes) {
        if (n.type === "bait") {
          const ptsHex = [];
          for (let k = 0; k < 6; k++) {
            const ang = (Math.PI / 3) * k - Math.PI / 6;
            ptsHex.push(
              `${(n.x + n.r * Math.cos(ang)).toFixed(1)},${(n.y + n.r * Math.sin(ang)).toFixed(1)}`
            );
          }
          svg.push(
            `<polygon points="${ptsHex.join(" ")}" fill="${n.color}" stroke="#0f172a" stroke-width="2.2"/>`,
            `<text x="${n.x.toFixed(1)}" y="${(n.y + n.r + 15).toFixed(
              1
            )}" text-anchor="middle" font-size="11.5" font-weight="700" fill="#0f172a">${escapeXml(
              n.label
            )}</text>`
          );
        } else {
          svg.push(
            `<circle cx="${n.x.toFixed(1)}" cy="${n.y.toFixed(1)}" r="${n.r}" fill="${
              n.sharedCount >= 2 ? "#fef3c7" : "#ffffff"
            }" stroke="${n.sharedCount >= 2 ? "#d97706" : n.primaryBaitColor}" stroke-width="${
              n.sharedCount >= 2 ? 2.5 : 2.0
            }"/>`,
            `<text x="${n.x.toFixed(1)}" y="${(n.y - n.r - 4).toFixed(
              1
            )}" text-anchor="middle" font-size="10.5" font-weight="${
              n.sharedCount >= 2 ? 600 : 500
            }" fill="#0f172a">${escapeXml(n.label)}</text>`
          );
        }
      }
      svg.push(`</svg>`);
      return svg.join("\n");
    }

    return { toSVG };
  }

  // ============================================================================
  // 4. REPLICATE BREAKDOWN & SAINT MODEL FIT INSPECTOR
  // ============================================================================

  function renderReplicateInspector(canvas, row, quant_type = "spc") {
    const width = Math.max(420, canvas.parentElement ? canvas.parentElement.clientWidth - 20 : 560);
    const height = 260;
    const ctx = setupHiDPICanvas(canvas, width, height);

    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, width, height);

    if (!row) {
      ctx.fillStyle = "#64748b";
      ctx.font = "500 13px Inter, system-ui, sans-serif";
      ctx.textAlign = "center";
      ctx.fillText("Click any row, dot, or network node to inspect replicates.", width / 2, height / 2);
      return { toSVG: () => "" };
    }

    const reps = row._replicatesRaw || [];
    const repNames = row._replicateNames || reps.map((_, i) => `Rep ${i + 1}`);
    const ctrls = row._controlsRaw || [];
    const ctrlNames = row._controlNames || ctrls.map((_, i) => `Ctrl ${i + 1}`);
    const zBits = row._zBits || [];
    const repScores = row._repScores || [];

    const bars = [];
    for (let i = 0; i < reps.length; i++) {
      bars.push({
        label: repNames[i] || `T${i + 1}`,
        val: Number.isNaN(reps[i]) ? 0 : reps[i],
        isMissing: Number.isNaN(reps[i]),
        isCtrl: false,
        z: Boolean(zBits[i]),
        pScore: repScores[i] !== undefined ? repScores[i] : null,
      });
    }
    const maxCtrlsToShow = Math.min(12, ctrls.length);
    for (let i = 0; i < maxCtrlsToShow; i++) {
      bars.push({
        label: ctrlNames[i] || `C${i + 1}`,
        val: Number.isNaN(ctrls[i]) ? 0 : ctrls[i],
        isMissing: Number.isNaN(ctrls[i]),
        isCtrl: true,
        z: false,
        pScore: null,
      });
    }

    const padL = 56;
    const padR = 24;
    const padT = 36;
    const padB = 60;
    const plotW = width - padL - padR;
    const plotH = height - padT - padB;

    const eta = row._eta ?? 0;
    const muTrue = (row._eta ?? 0) + (row._d ?? 0);

    let minVal = 0;
    let maxVal = Math.max(muTrue * 1.15, 1);
    for (const b of bars) {
      if (!b.isMissing) {
        if (b.val < minVal) minVal = b.val;
        if (b.val > maxVal) maxVal = b.val;
      }
    }
    if (quant_type === "int") {
      minVal = Math.min(minVal, eta - 1.0, 0);
      maxVal = Math.max(maxVal + 0.6, muTrue + 0.6);
    } else {
      minVal = 0;
      maxVal = maxVal * 1.16;
    }

    const toPxY = (v) => padT + plotH - ((v - minVal) / (maxVal - minVal || 1)) * plotH;
    const zeroY = toPxY(Math.max(minVal, 0));

    const slotW = plotW / Math.max(1, bars.length);
    const barW = Math.max(12, Math.min(38, slotW * 0.66));

    // Shaded background region for Test Replicates vs Control IPs
    const testSplitX = padL + reps.length * slotW;
    ctx.fillStyle = "rgba(239, 246, 255, 0.45)";
    ctx.fillRect(padL, padT, reps.length * slotW, plotH);
    ctx.fillStyle = "rgba(248, 250, 252, 0.75)";
    ctx.fillRect(testSplitX, padT, plotW - reps.length * slotW, plotH);

    // Divider line between Test and Controls
    if (reps.length > 0 && maxCtrlsToShow > 0) {
      ctx.save();
      ctx.setLineDash([3, 3]);
      ctx.strokeStyle = "#94a3b8";
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(testSplitX, padT);
      ctx.lineTo(testSplitX, padT + plotH);
      ctx.stroke();
      ctx.restore();
    }

    // Plot border
    ctx.strokeStyle = "#cbd5e1";
    ctx.lineWidth = 1.2;
    ctx.strokeRect(padL, padT, plotW, plotH);

    // Draw Bars
    for (let i = 0; i < bars.length; i++) {
      const b = bars[i];
      const cx = padL + i * slotW + slotW / 2;
      const topY = b.isMissing ? zeroY : toPxY(b.val);
      const y0 = Math.min(zeroY, topY);
      const bh = Math.max(2, Math.abs(zeroY - topY));

      if (b.isCtrl) {
        ctx.fillStyle = "#94a3b8";
      } else if (b.z) {
        ctx.fillStyle = "#10b981"; // Z = 1 (True interaction state)
      } else {
        ctx.fillStyle = "#f59e0b"; // Z = 0 in test replicate
      }

      if (!b.isMissing) {
        drawRoundedRect(ctx, cx - barW / 2, y0, barW, bh, 3);
        ctx.fill();
      }

      // Value label above bar
      ctx.fillStyle = "#0f172a";
      ctx.font = "600 10px 'JetBrains Mono', Inter, monospace";
      ctx.textAlign = "center";
      ctx.textBaseline = "bottom";
      const valStr = b.isMissing
        ? "NA"
        : quant_type === "spc"
        ? String(Math.round(b.val))
        : b.val.toFixed(2);
      ctx.fillText(valStr, cx, Math.min(y0 - 3, padT + plotH - 4));

      // Sample name rotated
      ctx.save();
      ctx.translate(cx, padT + plotH + 7);
      ctx.rotate(-Math.PI / 5);
      ctx.textAlign = "right";
      ctx.textBaseline = "top";
      ctx.fillStyle = b.isCtrl ? "#64748b" : "#0f172a";
      ctx.font = b.isCtrl ? "500 9.5px Inter, sans-serif" : "600 10px Inter, sans-serif";
      ctx.fillText(b.label.slice(0, 14), 0, 0);
      ctx.restore();
    }

    // Horizontal reference lines for SAINT background (eta_i) and True signal (eta_i + d_i)
    if (eta >= minVal && eta <= maxVal) {
      const py = toPxY(eta);
      ctx.save();
      ctx.setLineDash([4, 3]);
      ctx.strokeStyle = "#dc2626";
      ctx.lineWidth = 1.4;
      ctx.beginPath();
      ctx.moveTo(padL, py);
      ctx.lineTo(padL + plotW, py);
      ctx.stroke();
      ctx.restore();
      ctx.fillStyle = "#dc2626";
      ctx.font = "600 10px 'JetBrains Mono', Inter, sans-serif";
      ctx.textAlign = "right";
      ctx.textBaseline = "bottom";
      ctx.fillText(`ηᵢ (bg) = ${eta.toFixed(2)}`, padL + plotW - 6, py - 3);
    }

    if (muTrue >= minVal && muTrue <= maxVal) {
      const py = toPxY(muTrue);
      ctx.save();
      ctx.setLineDash([2, 2]);
      ctx.strokeStyle = "#059669";
      ctx.lineWidth = 1.4;
      ctx.beginPath();
      ctx.moveTo(padL, py);
      ctx.lineTo(padL + plotW, py);
      ctx.stroke();
      ctx.restore();
      ctx.fillStyle = "#059669";
      ctx.font = "600 10px 'JetBrains Mono', Inter, sans-serif";
      ctx.textAlign = "left";
      ctx.textBaseline = "bottom";
      ctx.fillText(`ηᵢ + dᵢ (true) = ${muTrue.toFixed(2)}`, padL + 6, py - 3);
    }

    // Header title & legend
    ctx.fillStyle = "#0f172a";
    ctx.font = "700 12px Inter, system-ui, sans-serif";
    ctx.textAlign = "left";
    ctx.textBaseline = "top";
    ctx.fillText(
      `${row.Bait} → ${row.PreyGene}   ·   SaintScore: ${row.SaintScore.toFixed(4)}   ·   BFDR: ${
        Number.isNaN(row.computedBFDR) ? "N/A" : row.computedBFDR.toFixed(4)
      }`,
      padL,
      10
    );

    // Mini key at top-right
    ctx.font = "600 10.5px Inter, system-ui, sans-serif";
    ctx.textAlign = "right";
    ctx.fillStyle = "#059669";
    ctx.fillText("■ Test (Z=1)   ■ Test (Z=0)   ■ Controls", width - padR, 10);
  }

  // ============================================================================
  // 5. SCORE DISTRIBUTION HISTOGRAM & ICM CONVERGENCE TRACE
  // ============================================================================

  function renderQcCharts(canvas, rows, summary) {
    const width = Math.max(520, canvas.parentElement ? canvas.parentElement.clientWidth - 24 : 660);
    const height = 280;
    const ctx = setupHiDPICanvas(canvas, width, height);

    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, width, height);

    const halfW = Math.floor((width - 48) / 2);
    const padT = 42;
    const padB = 46;
    const plotH = height - padT - padB;

    // Left Chart: SaintScore Bimodal Distribution (20 bins in [0, 1])
    const leftX = 48;
    const leftPlotW = halfW - 24;
    const bins = new Int32Array(20);
    for (const r of rows) {
      const b = Math.min(19, Math.max(0, Math.floor((r.SaintScore || 0) * 20)));
      bins[b]++;
    }
    const maxBin = Math.max(1, ...bins);

    ctx.fillStyle = "#f8fafc";
    ctx.fillRect(leftX, padT, leftPlotW, plotH);
    ctx.strokeStyle = "#cbd5e1";
    ctx.lineWidth = 1.2;
    ctx.strokeRect(leftX, padT, leftPlotW, plotH);

    const bw = leftPlotW / 20;
    for (let i = 0; i < 20; i++) {
      const bh = (bins[i] / maxBin) * (plotH - 16);
      const bx = leftX + i * bw;
      const by = padT + plotH - bh;
      ctx.fillStyle = i >= 16 ? "#059669" : i >= 10 ? "#d97706" : "#94a3b8";
      if (bh > 0) {
        drawRoundedRect(ctx, bx + 1.5, by, Math.max(2, bw - 3), bh, 2);
        ctx.fill();
      }
    }

    ctx.fillStyle = "#0f172a";
    ctx.font = "700 12px Inter, sans-serif";
    ctx.textAlign = "left";
    ctx.fillText("SaintScore (AvgP) Bimodal Distribution", leftX, 16);

    ctx.fillStyle = "#64748b";
    ctx.font = "500 10.5px Inter, sans-serif";
    ctx.fillText("0.0 (Background)", leftX, padT + plotH + 16);
    ctx.textAlign = "right";
    ctx.fillText("1.0 (True Interactors)", leftX + leftPlotW, padT + plotH + 16);

    // Right Chart: ICM Log-Likelihood Convergence Trace
    const rightX = leftX + halfW + 24;
    const rightW = width - rightX - 24;
    ctx.fillStyle = "#f8fafc";
    ctx.fillRect(rightX, padT, rightW, plotH);
    ctx.strokeStyle = "#cbd5e1";
    ctx.lineWidth = 1.2;
    ctx.strokeRect(rightX, padT, rightW, plotH);

    ctx.fillStyle = "#0f172a";
    ctx.font = "700 12px Inter, sans-serif";
    ctx.textAlign = "left";
    ctx.fillText(
      `ICM Convergence (${summary?.icmIterations || 0} iters · β₁ = ${(summary?.beta1 ?? 0).toFixed(3)})`,
      rightX,
      16
    );

    const trace = summary?.icmTrace || [];
    const lls = trace.length > 0 ? [trace[0].old_ll, ...trace.map((t) => t.new_ll)] : [];
    if (lls.length > 0) {
      const minLL = Math.min(...lls);
      const maxLL = Math.max(...lls);
      const spanLL = Math.max(1e-6, maxLL - minLL);

      ctx.beginPath();
      for (let i = 0; i < lls.length; i++) {
        const px = rightX + 18 + (i / Math.max(1, lls.length - 1)) * (rightW - 36);
        const py = padT + plotH - 16 - ((lls[i] - minLL) / spanLL) * (plotH - 32);
        if (i === 0) ctx.moveTo(px, py);
        else ctx.lineTo(px, py);
      }
      ctx.strokeStyle = "#2563eb";
      ctx.lineWidth = 2.4;
      ctx.stroke();

      for (let i = 0; i < lls.length; i++) {
        const px = rightX + 18 + (i / Math.max(1, lls.length - 1)) * (rightW - 36);
        const py = padT + plotH - 16 - ((lls[i] - minLL) / spanLL) * (plotH - 32);
        ctx.beginPath();
        ctx.arc(px, py, 4, 0, Math.PI * 2);
        ctx.fillStyle = "#1d4ed8";
        ctx.fill();
        ctx.strokeStyle = "#ffffff";
        ctx.lineWidth = 1.2;
        ctx.stroke();
      }

      ctx.fillStyle = "#64748b";
      ctx.font = "500 10.5px 'JetBrains Mono', Inter, monospace";
      ctx.textAlign = "left";
      ctx.fillText(`Start LL: ${lls[0].toFixed(1)}`, rightX + 6, padT + plotH + 16);
      ctx.textAlign = "right";
      ctx.fillText(`Final LL: ${lls[lls.length - 1].toFixed(1)}`, rightX + rightW - 6, padT + plotH + 16);
    }

    function toSVG() {
      const svg = [
        `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" font-family="Inter, -apple-system, sans-serif">`,
        `<rect width="${width}" height="${height}" fill="#ffffff"/>`,
        `<text x="${leftX}" y="24" font-size="12" font-weight="700" fill="#0f172a">SaintScore (AvgP) Bimodal Distribution</text>`,
        `<rect x="${leftX}" y="${padT}" width="${leftPlotW}" height="${plotH}" fill="#f8fafc" stroke="#cbd5e1"/>`,
      ];
      for (let i = 0; i < 20; i++) {
        const bh = (bins[i] / maxBin) * (plotH - 16);
        const bx = leftX + i * bw;
        const by = padT + plotH - bh;
        const col = i >= 16 ? "#059669" : i >= 10 ? "#d97706" : "#94a3b8";
        if (bh > 0) {
          svg.push(
            `<rect x="${(bx + 1.5).toFixed(1)}" y="${by.toFixed(1)}" width="${Math.max(
              2,
              bw - 3
            ).toFixed(1)}" height="${bh.toFixed(1)}" rx="2" fill="${col}"/>`
          );
        }
      }
      svg.push(
        `<text x="${rightX}" y="24" font-size="12" font-weight="700" fill="#0f172a">${escapeXml(
          `ICM Convergence (${summary?.icmIterations || 0} iters · β₁ = ${(summary?.beta1 ?? 0).toFixed(3)})`
        )}</text>`,
        `<rect x="${rightX}" y="${padT}" width="${rightW}" height="${plotH}" fill="#f8fafc" stroke="#cbd5e1"/>`
      );
      if (lls.length > 0) {
        const minLL = Math.min(...lls);
        const maxLL = Math.max(...lls);
        const spanLL = Math.max(1e-6, maxLL - minLL);
        const polyPts = [];
        for (let i = 0; i < lls.length; i++) {
          const px = rightX + 18 + (i / Math.max(1, lls.length - 1)) * (rightW - 36);
          const py = padT + plotH - 16 - ((lls[i] - minLL) / spanLL) * (plotH - 32);
          polyPts.push(`${px.toFixed(1)},${py.toFixed(1)}`);
        }
        svg.push(
          `<polyline points="${polyPts.join(" ")}" fill="none" stroke="#2563eb" stroke-width="2.4"/>`
        );
      }
      svg.push(`</svg>`);
      return svg.join("\n");
    }

    return { toSVG };
  }

  return {
    renderDotPlot,
    renderVolcanoPlot,
    renderNetworkGraph,
    renderReplicateInspector,
    renderQcCharts,
    exportCanvasPNG,
    downloadBlob,
  };
});
