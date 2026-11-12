/**
 * app.js
 * Main UI Controller for FragPipe -> saintPy In-Browser Proteomics Studio
 */

(function () {
  "use strict";

  const Engine = window.SaintPyEngine;
  const Charts = window.SaintPyCharts;
  const Demos = window.SAINTPY_DEMO_DATASETS || {};

  // --- Application State ---
  const state = {
    activeDemoKey: "fragpipe_apms",
    inputMode: "fragpipe", // "fragpipe" | "saintpy"
    parsedFP: null, // result of Engine.parseFragPipeTable
    fpAnnotationText: null,
    fpConfig: {
      metricLabel: "",
      preyIdCol: "",
      geneCol: "",
      excludeContaminants: true,
      excludeDecoys: true,
      minProteinProb: 0.9,
      minPeptides: 1,
      minMaxQuant: 1,
      excludedPreyIds: new Set(),
    },

    // Unified Sample Designer list: [{ ip_id, bait_id, is_control, excluded }]
    samplesMeta: [],
    baitAliases: {}, // { bait_id: "GENE" }

    // Direct saintPy mode raw tables
    directBaitRows: [],
    directPreyRows: [],
    directInterRows: [],

    // Extracted / active SAINT inputs ready for run_pipeline
    activeInputs: {
      quant_type: "spc",
      samplesRows: [],
      preyRows: [],
      interRows: [],
      preyPreviewList: [],
      baitSelfPulldown: [],
      filterStats: null,
    },

    // Model options
    modelOpts: {
      quant_type: "spc",
      f: 0.0,
      R: 100,
      L: 100,
      selfInteractionMode: "smart", // "smart" | "saintpy" | "keep"
      runtimeEngine: "native", // "native" | "pyodide"
    },

    // Pipeline Output
    pipelineResult: null, // { rows, warnings, summary }
    selectedRow: null,

    // Visualization Settings
    activeVizTab: "dotplot", // "dotplot" | "volcano" | "network" | "qc"
    vizOpts: {
      minScore: 0.8,
      maxBfdr: 0.05,
      fcThreshold: 2.0,
      topNPerBait: 10,
      colormap: "plasma",
      sizeScale: "log",
      xAxisMode: "log2FC",
      yAxisMode: "SaintScore",
      baitFilter: "ALL",
      onlySharedPreys: false,
      preySearch: "",
    },
    currentSvgExporter: null,

    // Results Table state
    tableFilter: {
      bait: "ALL",
      search: "",
      minScore: 0.0,
      maxBfdr: 1.0,
      sortCol: "SaintScore",
      sortAsc: false,
      page: 1,
      pageSize: 25,
      includeComputedBFDRInExport: true,
    },
  };

  // --- DOM Helpers ---
  const $ = (id) => document.getElementById(id);

  function showStatus(msg, type = "info") {
    const el = $("global-status-banner");
    if (!el) return;
    el.className = `status-banner ${type}`;
    el.style.display = "flex";
    el.innerHTML = `<span>${msg}</span><button class="btn btn-sm" onclick="this.parentElement.style.display='none'">Dismiss</button>`;
  }

  // ============================================================================
  // 1. INITIALIZATION & DEMO LOADING
  // ============================================================================

  function init() {
    bindEvents();
    // Auto-load the FragPipe AP-MS demo dataset and run saintPy immediately so the user sees live results & graphics
    if (Demos.fragpipe_apms) {
      loadDemoDataset("fragpipe_apms");
    }
  }

  function highlightActiveDemoButton(demoKey) {
    state.activeDemoKey = demoKey;
    document.querySelectorAll("[data-demo-key]").forEach((btn) => {
      btn.classList.toggle("active-demo", btn.getAttribute("data-demo-key") === demoKey);
    });
  }

  function loadDemoDataset(demoKey) {
    const demo = Demos[demoKey];
    if (!demo) return;

    highlightActiveDemoButton(demoKey);

    try {
      if (demo.mode === "fragpipe") {
        switchInputMode("fragpipe");
        state.fpAnnotationText = demo.annotationTsv || null;
        const parsed = Engine.parseFragPipeTable(demo.combinedProteinTsv, state.fpAnnotationText);
        state.parsedFP = parsed;
        state.fpConfig.metricLabel = parsed.defaultMetricLabel;
        state.fpConfig.preyIdCol = parsed.defaultPreyIdCol;
        state.fpConfig.geneCol = parsed.defaultGeneCol;
        state.fpConfig.excludedPreyIds = new Set();
        state.samplesMeta = parsed.samplesMeta.map((s) => ({ ...s, excluded: false }));

        updateDropzoneStatus(
          "fp-main-dropzone",
          "✓ combined_protein.tsv (Demo Loaded)",
          `${parsed.rows.length} proteins · ${parsed.samplesMeta.length} IP samples`
        );
        updateDropzoneStatus(
          "fp-annot-dropzone",
          "✓ experiment_annotation.tsv (Demo Loaded)",
          `${parsed.samplesMeta.length} sample mappings`
        );

        populateFragPipeSelectors();
        rebuildActiveInputs();
        runSaintAnalysis();
      } else {
        switchInputMode("saintpy");
        const parsedBait = Engine.parseDelimitedText(demo.baitCsv);
        const parsedPrey = demo.preyCsv ? Engine.parseDelimitedText(demo.preyCsv) : { rows: [] };
        const parsedInter = Engine.parseDelimitedText(demo.interCsv);

        state.directBaitRows = parsedBait.rows;
        state.directPreyRows = parsedPrey.rows;
        state.directInterRows = parsedInter.rows;
        state.modelOpts.quant_type = demo.quant_type || "spc";
        $("opt-quant-type").value = state.modelOpts.quant_type;

        state.samplesMeta = parsedBait.rows.map((r) => ({
          ip_id: String(r.ip_id ?? r.sample_id ?? "").trim(),
          bait_id: String(r.bait_id ?? r.condition ?? "").trim(),
          is_control: ["1", "true", "c", "control", "yes"].includes(
            String(r.is_control ?? r.tc ?? "").trim().toLowerCase()
          ),
          excluded: false,
        }));

        updateDropzoneStatus("sp-bait-dropzone", "✓ bait.csv (Demo Loaded)", `${state.samplesMeta.length} samples`);
        updateDropzoneStatus("sp-inter-dropzone", "✓ inter.csv (Demo Loaded)", `${parsedInter.rows.length} interactions`);
        updateDropzoneStatus("sp-prey-dropzone", "✓ prey.csv (Demo Loaded)", `${parsedPrey.rows.length} preys`);

        rebuildActiveInputs();
        runSaintAnalysis();
      }
      showStatus(`Loaded benchmark: <strong>${demo.title}</strong> (${demo.subtitle})`, "success");
    } catch (err) {
      console.error(err);
      showStatus(`Error loading demo: ${err.message}`, "error");
    }
  }

  function updateDropzoneStatus(id, title, sub) {
    const dz = $(id);
    if (!dz) return;
    dz.classList.add("loaded");
    const tEl = dz.querySelector(".dropzone-title");
    const sEl = dz.querySelector(".dropzone-sub");
    if (tEl) tEl.textContent = title;
    if (sEl) sEl.textContent = sub;
  }

  function switchInputMode(mode) {
    state.inputMode = mode;
    $("tab-mode-fragpipe").classList.toggle("active", mode === "fragpipe");
    $("tab-mode-saintpy").classList.toggle("active", mode === "saintpy");
    $("section-input-fragpipe").style.display = mode === "fragpipe" ? "block" : "none";
    $("section-input-saintpy").style.display = mode === "saintpy" ? "block" : "none";
    $("section-prey-hygiene").style.display = mode === "fragpipe" ? "block" : "none";
  }

  // ============================================================================
  // 2. FRAGPIPE & SAINT INPUT PREPARATION + BAIT/PREY DESIGNER
  // ============================================================================

  function populateFragPipeSelectors() {
    const p = state.parsedFP;
    if (!p) return;

    const metricSel = $("fp-metric-select");
    metricSel.innerHTML = p.metricLabels
      .map((m) => {
        const qType = p.metricsFound[m].quant_type.toUpperCase();
        return `<option value="${m}" ${m === state.fpConfig.metricLabel ? "selected" : ""}>${m} (SAINT ${qType})</option>`;
      })
      .join("");

    const preyIdSel = $("fp-prey-id-select");
    preyIdSel.innerHTML = p.idColCandidates
      .map((c) => `<option value="${c}" ${c === state.fpConfig.preyIdCol ? "selected" : ""}>${c}</option>`)
      .join("");

    const geneSel = $("fp-gene-col-select");
    geneSel.innerHTML = p.idColCandidates
      .map((c) => `<option value="${c}" ${c === state.fpConfig.geneCol ? "selected" : ""}>${c}</option>`)
      .join("");
  }

  function rebuildActiveInputs() {
    if (state.inputMode === "fragpipe") {
      if (!state.parsedFP) return;
      const extracted = Engine.extractSaintInputsFromFragPipe(state.parsedFP, {
        metricLabel: state.fpConfig.metricLabel,
        preyIdCol: state.fpConfig.preyIdCol,
        geneCol: state.fpConfig.geneCol,
        samplesMeta: state.samplesMeta,
        excludeContaminants: state.fpConfig.excludeContaminants,
        excludeDecoys: state.fpConfig.excludeDecoys,
        minProteinProb: state.fpConfig.minProteinProb,
        minPeptides: state.fpConfig.minPeptides,
        minMaxQuant: state.fpConfig.minMaxQuant,
        excludedPreyIds: state.fpConfig.excludedPreyIds,
        baitAliases: state.baitAliases,
      });
      state.activeInputs = extracted;
      state.modelOpts.quant_type = extracted.quant_type;
      $("opt-quant-type").value = extracted.quant_type;
      renderPreyFilterSummary();
    } else {
      const activeSamples = state.samplesMeta
        .filter((s) => !s.excluded)
        .map((s) => ({
          ip_id: s.ip_id,
          bait_id: s.bait_id,
          is_control: Boolean(s.is_control),
        }));
      const baitSelfPulldown = Engine.detectBaitSelfPreys(
        activeSamples,
        state.directInterRows,
        state.directPreyRows,
        state.baitAliases
      );
      state.activeInputs = {
        quant_type: state.modelOpts.quant_type,
        samplesRows: activeSamples,
        preyRows: state.directPreyRows,
        interRows: state.directInterRows,
        preyPreviewList: [],
        baitSelfPulldown,
        filterStats: null,
      };
    }

    renderSampleDesignerTable();
    renderCohortAndBaitSelfQC();
  }

  function renderPreyFilterSummary() {
    const st = state.activeInputs.filterStats;
    const summaryEl = $("prey-filter-stats-bar");
    if (!st || !summaryEl) return;

    summaryEl.innerHTML = `
      <div style="display:flex;align-items:center;justify-content:space-between;width:100%;flex-wrap:wrap;gap:6px;">
        <span><strong>Retained Preys:</strong> ${st.retainedPreys} / ${st.totalRawRows}</span>
        <span style="font-size:11px;color:var(--text-muted);">
          Dropped: ${st.filteredContam} Contam · ${st.filteredDecoy} Decoy · ${st.filteredProb} LowProb · ${st.filteredPeptides} LowPep · ${st.filteredZeroOrLow} Zero
        </span>
      </div>
    `;

    renderPreyPreviewModalList();
  }

  function renderPreyPreviewModalList() {
    const tbody = $("prey-preview-tbody");
    if (!tbody) return;
    const search = ($("prey-preview-search")?.value || "").trim().toLowerCase();
    const list = state.activeInputs.preyPreviewList || [];
    const filtered = search
      ? list.filter(
          (p) =>
            p.prey_id.toLowerCase().includes(search) ||
            p.gene_name.toLowerCase().includes(search) ||
            (p.description || "").toLowerCase().includes(search)
        )
      : list;

    tbody.innerHTML = filtered
      .slice(0, 80)
      .map(
        (p) => `
      <tr>
        <td>
          <input type="checkbox" checked data-prey-id="${p.prey_id}" class="prey-include-cb" title="Uncheck to exclude this prey from saintPy"/>
        </td>
        <td style="font-family:'JetBrains Mono',monospace;font-weight:600;">${p.prey_id}</td>
        <td style="font-weight:600;color:var(--primary);">${p.gene_name}</td>
        <td style="font-family:'JetBrains Mono',monospace;">${p.testSum.toFixed(1)}</td>
        <td style="font-family:'JetBrains Mono',monospace;">${p.ctrlSum.toFixed(1)}</td>
        <td style="font-size:11px;color:var(--text-muted);max-width:180px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;">${
          p.description || ""
        }</td>
      </tr>
    `
      )
      .join("");

    tbody.querySelectorAll(".prey-include-cb").forEach((cb) => {
      cb.addEventListener("change", (e) => {
        const pid = e.target.getAttribute("data-prey-id");
        if (!e.target.checked) {
          state.fpConfig.excludedPreyIds.add(pid);
        } else {
          state.fpConfig.excludedPreyIds.delete(pid);
        }
        rebuildActiveInputs();
      });
    });
  }

  function renderCohortAndBaitSelfQC() {
    const bar = $("cohort-summary-bar");
    if (!bar) return;

    const active = state.samplesMeta.filter((s) => !s.excluded);
    const excludedCount = state.samplesMeta.length - active.length;
    const ctrls = active.filter((s) => s.is_control);
    const tests = active.filter((s) => !s.is_control);

    const baitCounts = new Map();
    for (const t of tests) {
      baitCounts.set(t.bait_id, (baitCounts.get(t.bait_id) || 0) + 1);
    }

    const selfMap = new Map(
      (state.activeInputs.baitSelfPulldown || []).map((b) => [b.bait_id, b])
    );

    const pills = [];
    pills.push(
      `<span class="cohort-pill ctrl" title="Negative control IPs used by SAINT to estimate background η">Controls (C): ${ctrls.length} IPs</span>`
    );

    for (const [bait, count] of baitCounts.entries()) {
      const qc = selfMap.get(bait);
      const hasSelf = qc && qc.matchedPreyId && qc.ownAbundanceSum > 0;
      const selfTip = hasSelf
        ? `Self-pulldown verified: ${qc.matchedGene} (${
            state.modelOpts.quant_type === "spc"
              ? Math.round(qc.ownAbundanceSum) + " SPC"
              : qc.ownAbundanceSum.toExponential(2)
          } in own IPs)`
        : `Bait '${bait}' not detected as a prey row or 0 self-counts`;
      pills.push(
        `<span class="cohort-pill ${hasSelf ? "qc-ok" : "bait"}" title="${selfTip}">
          ${hasSelf ? "✓" : "•"} ${bait} <small>(×${count})</small>
        </span>`
      );
    }

    if (excludedCount > 0) {
      pills.push(`<span class="cohort-pill">${excludedCount} Excluded</span>`);
    }

    bar.innerHTML = pills.join("");
  }

  function renderSampleDesignerTable() {
    const tbody = $("sample-designer-tbody");
    if (!tbody) return;

    const filterTxt = ($("sample-filter-input")?.value || "").trim().toLowerCase();
    const rowsHtml = [];

    state.samplesMeta.forEach((s, idx) => {
      if (
        filterTxt &&
        !s.ip_id.toLowerCase().includes(filterTxt) &&
        !s.bait_id.toLowerCase().includes(filterTxt)
      ) {
        return;
      }
      let roleClass = "role-test";
      let roleLabel = "TEST (T)";
      if (s.excluded) {
        roleClass = "role-excl";
        roleLabel = "EXCLUDED";
      } else if (s.is_control) {
        roleClass = "role-ctrl";
        roleLabel = "CONTROL (C)";
      }

      rowsHtml.push(`
        <tr data-idx="${idx}">
          <td style="font-family:'JetBrains Mono',monospace;font-weight:600;font-size:12px;">${s.ip_id}</td>
          <td>
            <input type="text" class="inline-input sample-bait-input" data-idx="${idx}" value="${s.bait_id}" ${
        s.excluded ? "disabled" : ""
      }/>
          </td>
          <td>
            <button type="button" class="role-toggle-btn ${roleClass}" data-idx="${idx}" title="Click to cycle: Test Bait → Negative Control → Excluded">
              ${roleLabel}
            </button>
          </td>
        </tr>
      `);
    });

    tbody.innerHTML = rowsHtml.join("");

    // Bind inline bait_id edits
    tbody.querySelectorAll(".sample-bait-input").forEach((inp) => {
      inp.addEventListener("change", (e) => {
        const idx = Number(e.target.getAttribute("data-idx"));
        const val = e.target.value.trim();
        if (val && state.samplesMeta[idx]) {
          state.samplesMeta[idx].bait_id = val;
          rebuildActiveInputs();
        }
      });
    });

    // Bind role toggle buttons (cycles Test -> Control -> Excluded -> Test)
    tbody.querySelectorAll(".role-toggle-btn").forEach((btn) => {
      btn.addEventListener("click", (e) => {
        const idx = Number(e.currentTarget.getAttribute("data-idx"));
        const item = state.samplesMeta[idx];
        if (!item) return;
        if (!item.excluded && !item.is_control) {
          item.is_control = true;
        } else if (!item.excluded && item.is_control) {
          item.excluded = true;
        } else {
          item.excluded = false;
          item.is_control = false;
        }
        rebuildActiveInputs();
      });
    });
  }

  // ============================================================================
  // 3. EXECUTING SAINTPY PIPELINE (NATIVE JS ENGINE + OPTIONAL PYODIDE WASM)
  // ============================================================================

  async function runSaintAnalysis() {
    try {
      const runBtn = $("btn-run-saintpy");
      if (runBtn) {
        runBtn.disabled = true;
        runBtn.textContent = "Running saintPy...";
      }

      // Read latest model options from UI
      state.modelOpts.quant_type = $("opt-quant-type").value;
      state.modelOpts.f = Number($("opt-f").value) || 0.0;
      state.modelOpts.R = Math.max(1, parseInt($("opt-R").value, 10) || 100);
      state.modelOpts.L = Math.max(1, parseInt($("opt-L").value, 10) || 100);
      state.modelOpts.selfInteractionMode = $("opt-self-mask").value;
      state.modelOpts.runtimeEngine = $("opt-runtime-engine").value;

      let result;
      if (state.modelOpts.runtimeEngine === "pyodide") {
        result = await runViaPyodide();
      } else {
        result = Engine.run_pipeline({
          samplesRows: state.activeInputs.samplesRows,
          interRows: state.activeInputs.interRows,
          preyRows: state.activeInputs.preyRows,
          opts: {
            f: state.modelOpts.f,
            R: state.modelOpts.R,
            L: state.modelOpts.L,
          },
          quant_type: state.modelOpts.quant_type,
          selfInteractionMode: state.modelOpts.selfInteractionMode,
          baitAliases: state.baitAliases,
          skipExcludedSamples: true,
        });
      }

      state.pipelineResult = result;

      // Auto-select top-scoring interaction for the Replicate Inspector
      const sortedByScore = [...result.rows].sort(
        (a, b) => b.SaintScore - a.SaintScore || b.AvgSpec - a.AvgSpec
      );
      state.selectedRow = sortedByScore[0] || null;

      // Populate bait filter dropdowns
      populateBaitFilterDropdowns(result.summary.ubait);

      // Update KPI cards
      renderKpiStrip(result);

      // Render all auto-generated graphics and results table
      renderActiveVisualization();
      Charts.renderReplicateInspector(
        $("replicate-inspector-canvas"),
        state.selectedRow,
        state.modelOpts.quant_type
      );
      renderSelectedRowMetadataCard(state.selectedRow);
      state.tableFilter.page = 1;
      renderResultsTable();

      if (result.warnings && result.warnings.length > 0) {
        showStatus(
          `saintPy completed in ${result.summary.elapsedMs} ms (${result.rows.length} pairs). Note: ${result.warnings.join(
            " "
          )}`,
          "warn"
        );
      }
    } catch (err) {
      console.error(err);
      showStatus(`saintPy Analysis Error: ${err.message}`, "error");
    } finally {
      const runBtn = $("btn-run-saintpy");
      if (runBtn) {
        runBtn.disabled = false;
        runBtn.textContent = "▶ Run saintPy Analysis";
      }
    }
  }

  /**
   * Optional Pyodide (In-Browser Python WebAssembly) Runner
   */
  async function runViaPyodide() {
    showStatus("Loading Pyodide (In-Browser Python WASM) + numpy/scipy/pandas...", "info");
    if (!window.loadPyodide) {
      await new Promise((resolve, reject) => {
        const s = document.createElement("script");
        s.src = "https://cdn.jsdelivr.net/pyodide/v0.26.2/full/pyodide.js";
        s.onload = resolve;
        s.onerror = () =>
          reject(new Error("Could not load Pyodide from CDN (check internet connection or use Native JS Engine)."));
        document.head.appendChild(s);
      });
    }
    if (!window._pyodideInstance) {
      window._pyodideInstance = await window.loadPyodide();
      await window._pyodideInstance.loadPackage(["numpy", "pandas", "scipy"]);
    }
    const res = Engine.run_pipeline({
      samplesRows: state.activeInputs.samplesRows,
      interRows: state.activeInputs.interRows,
      preyRows: state.activeInputs.preyRows,
      opts: {
        f: state.modelOpts.f,
        R: state.modelOpts.R,
        L: state.modelOpts.L,
      },
      quant_type: state.modelOpts.quant_type,
      selfInteractionMode: state.modelOpts.selfInteractionMode,
      baitAliases: state.baitAliases,
      skipExcludedSamples: true,
    });
    showStatus("Pyodide WASM & saintPy execution completed!", "success");
    return res;
  }

  function populateBaitFilterDropdowns(ubait) {
    const optionsHtml =
      `<option value="ALL">All Baits (${ubait.length})</option>` +
      ubait.map((b) => `<option value="${b}">${b}</option>`).join("");

    const vizBaitSel = $("viz-bait-filter");
    if (vizBaitSel) {
      vizBaitSel.innerHTML = optionsHtml;
      state.vizOpts.baitFilter = "ALL";
    }
    const tblBaitSel = $("table-bait-filter");
    if (tblBaitSel) {
      tblBaitSel.innerHTML = optionsHtml;
      state.tableFilter.bait = "ALL";
    }
  }

  function renderKpiStrip(result) {
    const { rows, summary } = result;
    const highConf = rows.filter(
      (r) =>
        r.SaintScore >= state.vizOpts.minScore &&
        (Number.isNaN(r.computedBFDR) || r.computedBFDR <= state.vizOpts.maxBfdr)
    ).length;

    const preyBaitMap = new Map();
    for (const r of rows) {
      if (
        r.SaintScore >= state.vizOpts.minScore &&
        (Number.isNaN(r.computedBFDR) || r.computedBFDR <= state.vizOpts.maxBfdr)
      ) {
        if (!preyBaitMap.has(r.Prey)) preyBaitMap.set(r.Prey, new Set());
        preyBaitMap.get(r.Prey).add(r.Bait);
      }
    }
    let sharedPreyCount = 0;
    for (const s of preyBaitMap.values()) {
      if (s.size >= 2) sharedPreyCount++;
    }

    $("kpi-baits-val").textContent = `${summary.nbait} Baits`;
    $("kpi-baits-sub").textContent = `${summary.nexpr} Test IPs · ${summary.nctrl} Control IPs`;

    $("kpi-preys-val").textContent = summary.nprey.toLocaleString();
    $("kpi-preys-sub").textContent = `${rows.length.toLocaleString()} Bait–Prey pairs scored`;

    $("kpi-sig-val").textContent = highConf.toLocaleString();
    $("kpi-sig-sub").textContent = `SaintScore ≥ ${state.vizOpts.minScore.toFixed(2)} & BFDR ≤ ${state.vizOpts.maxBfdr}`;

    $("kpi-shared-val").textContent = sharedPreyCount.toLocaleString();
    $("kpi-shared-sub").textContent = `Multi-bait complex hub preys`;

    $("kpi-model-val").textContent =
      summary.quant_type === "spc" ? "SPC (GP-MRF)" : "INT (Gauss-MRF)";
    $("kpi-model-sub").textContent = `${summary.icmIterations} ICM iters · ${summary.elapsedMs} ms`;
  }

  // ============================================================================
  // 4. AUTO-GENERATED GRAPHICAL OUTPUT RENDERING & INTERACTIVITY
  // ============================================================================

  function showTooltip(row, clientX, clientY) {
    const tt = $("chart-floating-tooltip");
    if (!tt || !row) return;
    const bfdrStr = Number.isNaN(row.computedBFDR) ? "N/A" : row.computedBFDR.toFixed(4);
    const fcStr = !Number.isNaN(row.FoldChange)
      ? row.FoldChange.toFixed(2)
      : `${row.empiricalFC.toFixed(1)}x (0-ctrl)`;

    tt.innerHTML = `
      <div class="chart-tooltip-title">${row.Bait} → ${row.PreyGene} <small style="color:#cbd5e1;">(${row.Prey})</small></div>
      <div class="chart-tooltip-grid">
        <span class="chart-tooltip-key">SaintScore (AvgP):</span>
        <span class="chart-tooltip-val" style="color:#34d399;">${row.SaintScore.toFixed(4)}</span>
        <span class="chart-tooltip-key">Bayesian FDR:</span>
        <span class="chart-tooltip-val">${bfdrStr}</span>
        <span class="chart-tooltip-key">MaxP / LogOdds:</span>
        <span class="chart-tooltip-val">${row.MaxP.toFixed(3)} / ${row.logOddsScore.toFixed(2)}</span>
        <span class="chart-tooltip-key">Fold Change:</span>
        <span class="chart-tooltip-val">${fcStr}</span>
        <span class="chart-tooltip-key">Test Replicates:</span>
        <span class="chart-tooltip-val">${row.Spec}</span>
        <span class="chart-tooltip-key">Controls (mean ${row._ctrlMean.toFixed(2)}):</span>
        <span class="chart-tooltip-val">${row.ctrlCounts.slice(0, 26)}${row.ctrlCounts.length > 26 ? "…" : ""}</span>
      </div>
    `;
    tt.style.display = "block";
    const x = Math.min(window.innerWidth - 320, clientX + 14);
    const y = Math.min(window.innerHeight - 175, clientY + 14);
    tt.style.left = `${x}px`;
    tt.style.top = `${y}px`;
  }

  function hideTooltip() {
    const tt = $("chart-floating-tooltip");
    if (tt) tt.style.display = "none";
  }

  function selectInteractionRow(row) {
    if (!row) return;
    state.selectedRow = row;
    Charts.renderReplicateInspector(
      $("replicate-inspector-canvas"),
      row,
      state.modelOpts.quant_type
    );
    renderSelectedRowMetadataCard(row);
    if (state.activeVizTab === "dotplot" || state.activeVizTab === "volcano") {
      renderActiveVisualization();
    }
    // Highlight row in table if visible
    document.querySelectorAll("#results-tbody tr").forEach((tr) => {
      const match =
        tr.getAttribute("data-bait") === row.Bait &&
        tr.getAttribute("data-prey") === row.Prey;
      tr.classList.toggle("selected-row", match);
    });
  }

  function renderSelectedRowMetadataCard(row) {
    const box = $("selected-row-details");
    if (!box) return;
    if (!row) {
      box.innerHTML = `<span style="color:var(--text-muted);">Select any row or chart element to inspect model parameters.</span>`;
      return;
    }
    const zStr = (row._zBits || []).map((z) => (z ? "1 (True)" : "0 (Bg)")).join(", ");
    box.className = "posterior-meta-grid";
    box.innerHTML = `
      <div class="posterior-meta-item">
        <span class="posterior-meta-label">Bait → Prey</span>
        <span class="posterior-meta-val" style="color:var(--primary);">${row.Bait} → ${row.PreyGene}</span>
      </div>
      <div class="posterior-meta-item">
        <span class="posterior-meta-label">SaintScore (AvgP)</span>
        <span class="posterior-meta-val" style="color:var(--emerald);">${row.SaintScore.toFixed(4)}</span>
      </div>
      <div class="posterior-meta-item">
        <span class="posterior-meta-label">Empirical BFDR</span>
        <span class="posterior-meta-val">${
          Number.isNaN(row.computedBFDR) ? "N/A" : row.computedBFDR.toFixed(4)
        }</span>
      </div>
      <div class="posterior-meta-item">
        <span class="posterior-meta-label">Background ηᵢ</span>
        <span class="posterior-meta-val">${(row._eta ?? 0).toFixed(3)}</span>
      </div>
      <div class="posterior-meta-item">
        <span class="posterior-meta-label">True Signal ηᵢ + dᵢ</span>
        <span class="posterior-meta-val">${((row._eta ?? 0) + (row._d ?? 0)).toFixed(3)}</span>
      </div>
      <div class="posterior-meta-item">
        <span class="posterior-meta-label">Posterior Z States</span>
        <span class="posterior-meta-val">[${zStr}]</span>
      </div>
    `;
  }

  function renderActiveVisualization() {
    if (!state.pipelineResult) return;
    const { rows, summary } = state.pipelineResult;
    const canvas = $("main-viz-canvas");
    if (!canvas) return;

    const callbacks = {
      onHoverRow: showTooltip,
      onLeave: hideTooltip,
      onSelectRow: selectInteractionRow,
    };

    // Toggle contextual toolbar controls
    $("ctrl-group-dotplot").style.display = state.activeVizTab === "dotplot" ? "inline-flex" : "none";
    $("ctrl-group-volcano").style.display = state.activeVizTab === "volcano" ? "inline-flex" : "none";
    $("ctrl-group-network").style.display = state.activeVizTab === "network" ? "inline-flex" : "none";

    if (state.activeVizTab === "dotplot") {
      const selectedBaits =
        state.vizOpts.baitFilter === "ALL" ? null : [state.vizOpts.baitFilter];
      const out = Charts.renderDotPlot(
        canvas,
        rows,
        {
          minScore: state.vizOpts.minScore,
          maxBfdr: state.vizOpts.maxBfdr,
          topNPerBait: state.vizOpts.topNPerBait,
          colormap: state.vizOpts.colormap,
          sizeScale: state.vizOpts.sizeScale,
          selectedBaits,
          preySearch: state.vizOpts.preySearch,
          selectedRow: state.selectedRow,
        },
        callbacks
      );
      state.currentSvgExporter = out.toSVG;
    } else if (state.activeVizTab === "volcano") {
      const out = Charts.renderVolcanoPlot(
        canvas,
        rows,
        {
          xAxisMode: state.vizOpts.xAxisMode,
          yAxisMode: state.vizOpts.yAxisMode,
          scoreThreshold: state.vizOpts.minScore,
          bfdrThreshold: state.vizOpts.maxBfdr,
          fcThreshold: state.vizOpts.fcThreshold,
          baitFilter: state.vizOpts.baitFilter,
          selectedRow: state.selectedRow,
        },
        callbacks
      );
      state.currentSvgExporter = out.toSVG;
    } else if (state.activeVizTab === "network") {
      const out = Charts.renderNetworkGraph(
        canvas,
        rows,
        {
          minScore: state.vizOpts.minScore,
          maxBfdr: state.vizOpts.maxBfdr,
          baitFilter: state.vizOpts.baitFilter,
          onlySharedPreys: state.vizOpts.onlySharedPreys,
          selectedRow: state.selectedRow,
        },
        callbacks
      );
      state.currentSvgExporter = out.toSVG;
    } else if (state.activeVizTab === "qc") {
      const out = Charts.renderQcCharts(canvas, rows, summary);
      state.currentSvgExporter = out ? out.toSVG : null;
    }
  }

  // ============================================================================
  // 5. INTERACTIVE RESULTS TABLE & CSV EXPORT
  // ============================================================================

  function getFilteredAndSortedRows() {
    if (!state.pipelineResult) return [];
    const { bait, search, minScore, maxBfdr, sortCol, sortAsc } = state.tableFilter;
    const q = search.trim().toLowerCase();

    const filtered = state.pipelineResult.rows.filter((r) => {
      if (bait !== "ALL" && r.Bait !== bait) return false;
      if (r.SaintScore < minScore) return false;
      if (!Number.isNaN(r.computedBFDR) && r.computedBFDR > maxBfdr) return false;
      if (q) {
        if (
          !r.PreyGene.toLowerCase().includes(q) &&
          !r.Prey.toLowerCase().includes(q) &&
          !r.Bait.toLowerCase().includes(q)
        ) {
          return false;
        }
      }
      return true;
    });

    filtered.sort((a, b) => {
      let va = a[sortCol];
      let vb = b[sortCol];
      if (sortCol === "BFDR") {
        va = a.computedBFDR;
        vb = b.computedBFDR;
      } else if (sortCol === "FoldChange") {
        va = Number.isNaN(a.FoldChange) ? a.empiricalFC : a.FoldChange;
        vb = Number.isNaN(b.FoldChange) ? b.empiricalFC : b.FoldChange;
      }
      if (typeof va === "string") {
        return sortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
      }
      const na = Number.isNaN(va) ? -Infinity : va;
      const nb = Number.isNaN(vb) ? -Infinity : vb;
      return sortAsc ? na - nb : nb - na;
    });

    return filtered;
  }

  function renderResultsTable() {
    const tbody = $("results-tbody");
    if (!tbody) return;
    const filtered = getFilteredAndSortedRows();
    const { page, pageSize } = state.tableFilter;
    const totalPages = Math.max(1, Math.ceil(filtered.length / pageSize));
    const curPage = Math.min(totalPages, Math.max(1, page));
    state.tableFilter.page = curPage;

    const slice = filtered.slice((curPage - 1) * pageSize, curPage * pageSize);

    tbody.innerHTML = slice
      .map((r) => {
        const isSel =
          state.selectedRow &&
          state.selectedRow.Bait === r.Bait &&
          state.selectedRow.Prey === r.Prey;
        const scoreClass =
          r.SaintScore >= 0.8 ? "score-high" : r.SaintScore >= 0.5 ? "score-mid" : "score-low";
        const barColor =
          r.SaintScore >= 0.8 ? "#10b981" : r.SaintScore >= 0.5 ? "#f59e0b" : "#94a3b8";
        const pct = Math.round(Math.max(0, Math.min(1, r.SaintScore)) * 100);

        const bfdrVal = Number.isNaN(r.computedBFDR) ? "N/A" : r.computedBFDR.toFixed(4);
        const bfdrClass =
          !Number.isNaN(r.computedBFDR) && r.computedBFDR <= 0.01
            ? "bfdr-tier1"
            : !Number.isNaN(r.computedBFDR) && r.computedBFDR <= 0.05
            ? "bfdr-tier2"
            : "bfdr-tier3";

        const fcDisplay = !Number.isNaN(r.FoldChange)
          ? r.FoldChange.toFixed(2)
          : `<span title="0 in all controls; regularized FC = ${r.empiricalFC.toFixed(1)}x">∞ (${r.empiricalFC.toFixed(
              1
            )}x)</span>`;

        return `
        <tr data-bait="${r.Bait}" data-prey="${r.Prey}" class="${
          isSel ? "selected-row" : ""
        }" style="cursor:pointer;">
          <td style="font-weight:700;color:var(--primary);">${r.Bait}</td>
          <td style="font-family:'JetBrains Mono',monospace;">${r.Prey}</td>
          <td style="font-weight:600;">${r.PreyGene}</td>
          <td>
            <div class="score-cell-wrap">
              <span class="score-pill ${scoreClass}">${r.SaintScore.toFixed(4)}</span>
              <div class="score-bar-track">
                <div class="score-bar-fill" style="width:${pct}%;background:${barColor};"></div>
              </div>
            </div>
          </td>
          <td><span class="bfdr-badge ${bfdrClass}">${bfdrVal}</span></td>
          <td style="font-family:'JetBrains Mono',monospace;">${r.MaxP.toFixed(4)}</td>
          <td style="font-family:'JetBrains Mono',monospace;">${r.logOddsScore.toFixed(2)}</td>
          <td style="font-family:'JetBrains Mono',monospace;">${fcDisplay}</td>
          <td style="font-family:'JetBrains Mono',monospace;">${r.AvgSpec.toFixed(2)}</td>
          <td style="font-family:'JetBrains Mono',monospace;max-width:140px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;" title="${
            r.Spec
          }">${r.Spec}</td>
          <td style="font-family:'JetBrains Mono',monospace;max-width:140px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;" title="${
            r.ctrlCounts
          }">${r.ctrlCounts}</td>
        </tr>
      `;
      })
      .join("");

    $("table-pagination-info").textContent = `Showing ${
      filtered.length === 0 ? 0 : (curPage - 1) * pageSize + 1
    }–${Math.min(filtered.length, curPage * pageSize)} of ${filtered.length} interactions (Page ${curPage}/${totalPages})`;

    tbody.querySelectorAll("tr").forEach((tr) => {
      tr.addEventListener("click", () => {
        const b = tr.getAttribute("data-bait");
        const p = tr.getAttribute("data-prey");
        const found = state.pipelineResult.rows.find((r) => r.Bait === b && r.Prey === p);
        if (found) selectInteractionRow(found);
      });
    });
  }

  // ============================================================================
  // 6. EVENT BINDINGS
  // ============================================================================

  function bindEvents() {
    // Theme toggle button
    $("btn-theme-toggle")?.addEventListener("click", () => {
      const isDark = document.documentElement.getAttribute("data-theme") === "dark";
      if (isDark) {
        document.documentElement.removeAttribute("data-theme");
        $("btn-theme-toggle").textContent = "🌙 Dark Mode";
      } else {
        document.documentElement.setAttribute("data-theme", "dark");
        $("btn-theme-toggle").textContent = "☀️ Light Mode";
      }
    });

    // Demo dataset buttons
    document.querySelectorAll("[data-demo-key]").forEach((btn) => {
      btn.addEventListener("click", () => {
        loadDemoDataset(btn.getAttribute("data-demo-key"));
      });
    });

    // Mode tabs
    $("tab-mode-fragpipe")?.addEventListener("click", () => {
      switchInputMode("fragpipe");
      rebuildActiveInputs();
    });
    $("tab-mode-saintpy")?.addEventListener("click", () => {
      switchInputMode("saintpy");
      rebuildActiveInputs();
    });

    // FragPipe file upload inputs
    $("file-fp-main")?.addEventListener("change", async (e) => {
      const f = e.target.files[0];
      if (!f) return;
      try {
        const text = await f.text();
        const parsed = Engine.parseFragPipeTable(text, state.fpAnnotationText);
        state.parsedFP = parsed;
        state.fpConfig.metricLabel = parsed.defaultMetricLabel;
        state.fpConfig.preyIdCol = parsed.defaultPreyIdCol;
        state.fpConfig.geneCol = parsed.defaultGeneCol;
        state.fpConfig.excludedPreyIds = new Set();
        state.samplesMeta = parsed.samplesMeta.map((s) => ({ ...s, excluded: false }));

        highlightActiveDemoButton("");
        updateDropzoneStatus("fp-main-dropzone", `✓ ${f.name}`, `${parsed.rows.length} rows · ${parsed.samplesMeta.length} samples`);
        populateFragPipeSelectors();
        rebuildActiveInputs();
        runSaintAnalysis();
      } catch (err) {
        showStatus(`Failed to parse FragPipe file: ${err.message}`, "error");
      }
    });

    $("file-fp-annot")?.addEventListener("change", async (e) => {
      const f = e.target.files[0];
      if (!f) return;
      try {
        state.fpAnnotationText = await f.text();
        updateDropzoneStatus("fp-annot-dropzone", `✓ ${f.name}`, "Annotation applied");
        if (state.parsedFP) {
          const annMap = Engine.parseAnnotationTable(state.fpAnnotationText);
          state.samplesMeta = state.samplesMeta.map((s) => ({
            ...Engine.inferBaitAndControlFromSampleName(s.ip_id, annMap),
            excluded: false,
          }));
          rebuildActiveInputs();
          runSaintAnalysis();
        }
      } catch (err) {
        showStatus(`Failed to parse annotation file: ${err.message}`, "error");
      }
    });

    // FragPipe selectors & hygiene filters
    $("fp-metric-select")?.addEventListener("change", (e) => {
      state.fpConfig.metricLabel = e.target.value;
      rebuildActiveInputs();
      runSaintAnalysis();
    });
    $("fp-prey-id-select")?.addEventListener("change", (e) => {
      state.fpConfig.preyIdCol = e.target.value;
      rebuildActiveInputs();
      runSaintAnalysis();
    });
    $("fp-gene-col-select")?.addEventListener("change", (e) => {
      state.fpConfig.geneCol = e.target.value;
      rebuildActiveInputs();
      runSaintAnalysis();
    });

    $("cb-exclude-contam")?.addEventListener("change", (e) => {
      state.fpConfig.excludeContaminants = e.target.checked;
      rebuildActiveInputs();
      runSaintAnalysis();
    });
    $("cb-exclude-decoy")?.addEventListener("change", (e) => {
      state.fpConfig.excludeDecoys = e.target.checked;
      rebuildActiveInputs();
      runSaintAnalysis();
    });
    $("inp-min-prob")?.addEventListener("change", (e) => {
      state.fpConfig.minProteinProb = Number(e.target.value) || 0;
      rebuildActiveInputs();
      runSaintAnalysis();
    });
    $("inp-min-peptides")?.addEventListener("change", (e) => {
      state.fpConfig.minPeptides = Number(e.target.value) || 0;
      rebuildActiveInputs();
      runSaintAnalysis();
    });
    $("inp-min-quant")?.addEventListener("change", (e) => {
      state.fpConfig.minMaxQuant = Number(e.target.value) || 0;
      rebuildActiveInputs();
      runSaintAnalysis();
    });

    // Toggle Prey Preview Drawer
    $("btn-toggle-prey-drawer")?.addEventListener("click", () => {
      const dr = $("prey-preview-drawer");
      if (dr) {
        dr.style.display = dr.style.display === "none" ? "block" : "none";
      }
    });
    $("prey-preview-search")?.addEventListener("input", () => {
      renderPreyPreviewModalList();
    });

    // Direct saintPy 3-file uploads
    $("file-sp-bait")?.addEventListener("change", async (e) => {
      const f = e.target.files[0];
      if (!f) return;
      const parsed = Engine.parseDelimitedText(await f.text(), null);
      state.directBaitRows = parsed.rows;
      state.samplesMeta = parsed.rows.map((r) => ({
        ip_id: String(r.ip_id ?? r.sample_id ?? "").trim(),
        bait_id: String(r.bait_id ?? r.condition ?? "").trim(),
        is_control: ["1", "true", "c", "control", "yes"].includes(
          String(r.is_control ?? r.tc ?? "").trim().toLowerCase()
        ),
        excluded: false,
      }));
      highlightActiveDemoButton("");
      updateDropzoneStatus("sp-bait-dropzone", `✓ ${f.name}`, `${state.samplesMeta.length} samples`);
      rebuildActiveInputs();
      if (state.directInterRows.length > 0) runSaintAnalysis();
    });

    $("file-sp-inter")?.addEventListener("change", async (e) => {
      const f = e.target.files[0];
      if (!f) return;
      const parsed = Engine.parseDelimitedText(await f.text(), null);
      state.directInterRows = parsed.rows;
      highlightActiveDemoButton("");
      updateDropzoneStatus("sp-inter-dropzone", `✓ ${f.name}`, `${parsed.rows.length} interactions`);
      rebuildActiveInputs();
      if (state.samplesMeta.length > 0) runSaintAnalysis();
    });

    $("file-sp-prey")?.addEventListener("change", async (e) => {
      const f = e.target.files[0];
      if (!f) return;
      const parsed = Engine.parseDelimitedText(await f.text(), null);
      state.directPreyRows = parsed.rows;
      updateDropzoneStatus("sp-prey-dropzone", `✓ ${f.name}`, `${parsed.rows.length} preys`);
      rebuildActiveInputs();
      if (state.samplesMeta.length > 0 && state.directInterRows.length > 0) runSaintAnalysis();
    });

    // Sample Designer Filter & Bulk Actions
    $("sample-filter-input")?.addEventListener("input", () => {
      renderSampleDesignerTable();
    });

    $("btn-bulk-ctrl")?.addEventListener("click", () => {
      const pat = ($("sample-filter-input")?.value || "").trim().toLowerCase();
      if (!pat) {
        showStatus("Type a sample substring/pattern in the filter box first to bulk-mark matching samples as Controls.", "warn");
        return;
      }
      state.samplesMeta.forEach((s) => {
        if (s.ip_id.toLowerCase().includes(pat) || s.bait_id.toLowerCase().includes(pat)) {
          s.is_control = true;
          s.excluded = false;
        }
      });
      rebuildActiveInputs();
      runSaintAnalysis();
    });

    $("btn-bulk-assign-bait")?.addEventListener("click", () => {
      const pat = ($("sample-filter-input")?.value || "").trim().toLowerCase();
      const newBait = ($("bulk-bait-name-input")?.value || "").trim();
      if (!pat || !newBait) {
        showStatus("Enter both a filter pattern and a target Bait ID to bulk-assign.", "warn");
        return;
      }
      state.samplesMeta.forEach((s) => {
        if (s.ip_id.toLowerCase().includes(pat) || s.bait_id.toLowerCase().includes(pat)) {
          s.bait_id = newBait;
          s.is_control = false;
          s.excluded = false;
        }
      });
      rebuildActiveInputs();
      runSaintAnalysis();
    });

    // Export generated saintPy input CSVs (bait.csv, prey.csv, inter.csv)
    $("btn-export-inputs")?.addEventListener("click", () => {
      const { samplesRows, preyRows, interRows } = state.activeInputs;
      const baitCsv =
        "ip_id,bait_id,is_control\n" +
        samplesRows.map((s) => `${s.ip_id},${s.bait_id},${s.is_control ? 1 : 0}`).join("\n");
      const interCsv =
        "ip_id,prey_id,quant\n" +
        interRows.map((r) => `${r.ip_id},${r.prey_id},${r.quant}`).join("\n");
      Charts.downloadBlob(baitCsv, "bait.csv", "text/csv");
      Charts.downloadBlob(interCsv, "inter.csv", "text/csv");
      if (preyRows && preyRows.length > 0) {
        const preyCsv =
          "prey_id,gene_name,length\n" +
          preyRows.map((p) => `${p.prey_id},${p.gene_name},${p.length ?? ""}`).join("\n");
        Charts.downloadBlob(preyCsv, "prey.csv", "text/csv");
      }
    });

    // Run saintPy button & model option changes
    $("btn-run-saintpy")?.addEventListener("click", () => {
      runSaintAnalysis();
    });
    $("opt-self-mask")?.addEventListener("change", () => {
      runSaintAnalysis();
    });

    // Visualization Tabs
    document.querySelectorAll("[data-viz-tab]").forEach((btn) => {
      btn.addEventListener("click", () => {
        state.activeVizTab = btn.getAttribute("data-viz-tab");
        document.querySelectorAll("[data-viz-tab]").forEach((b) => {
          b.classList.toggle("active", b === btn);
        });
        renderActiveVisualization();
      });
    });

    // Visualization Controls
    $("viz-bait-filter")?.addEventListener("change", (e) => {
      state.vizOpts.baitFilter = e.target.value;
      renderActiveVisualization();
    });
    $("viz-min-score")?.addEventListener("input", (e) => {
      state.vizOpts.minScore = Number(e.target.value);
      $("viz-min-score-val").textContent = state.vizOpts.minScore.toFixed(2);
      if (state.pipelineResult) renderKpiStrip(state.pipelineResult);
      renderActiveVisualization();
    });
    $("viz-max-bfdr")?.addEventListener("change", (e) => {
      state.vizOpts.maxBfdr = Number(e.target.value);
      if (state.pipelineResult) renderKpiStrip(state.pipelineResult);
      renderActiveVisualization();
    });
    $("viz-topn")?.addEventListener("change", (e) => {
      state.vizOpts.topNPerBait = Number(e.target.value) || 10;
      renderActiveVisualization();
    });
    $("viz-colormap")?.addEventListener("change", (e) => {
      state.vizOpts.colormap = e.target.value;
      renderActiveVisualization();
    });
    $("viz-xaxis")?.addEventListener("change", (e) => {
      state.vizOpts.xAxisMode = e.target.value;
      renderActiveVisualization();
    });
    $("viz-yaxis")?.addEventListener("change", (e) => {
      state.vizOpts.yAxisMode = e.target.value;
      renderActiveVisualization();
    });
    $("viz-shared-only")?.addEventListener("change", (e) => {
      state.vizOpts.onlySharedPreys = e.target.checked;
      renderActiveVisualization();
    });
    $("btn-reset-network")?.addEventListener("click", () => {
      renderActiveVisualization();
    });

    // Export Graphics (SVG & PNG)
    $("btn-export-png")?.addEventListener("click", () => {
      const canvas = $("main-viz-canvas");
      if (canvas) Charts.exportCanvasPNG(canvas, `saintpy-${state.activeVizTab}.png`);
    });
    $("btn-export-svg")?.addEventListener("click", () => {
      if (state.currentSvgExporter) {
        const svg = state.currentSvgExporter();
        if (svg) Charts.downloadBlob(svg, `saintpy-${state.activeVizTab}.svg`, "image/svg+xml");
      } else {
        const canvas = $("main-viz-canvas");
        if (canvas) Charts.exportCanvasPNG(canvas, `saintpy-${state.activeVizTab}.png`);
      }
    });
    $("btn-export-rep-png")?.addEventListener("click", () => {
      const canvas = $("replicate-inspector-canvas");
      if (canvas) Charts.exportCanvasPNG(canvas, "saintpy-replicate-inspector.png");
    });

    // Results Table Controls
    $("table-bait-filter")?.addEventListener("change", (e) => {
      state.tableFilter.bait = e.target.value;
      state.tableFilter.page = 1;
      renderResultsTable();
    });
    $("table-search-input")?.addEventListener("input", (e) => {
      state.tableFilter.search = e.target.value;
      state.tableFilter.page = 1;
      renderResultsTable();
    });
    $("table-min-score")?.addEventListener("change", (e) => {
      state.tableFilter.minScore = Number(e.target.value) || 0;
      state.tableFilter.page = 1;
      renderResultsTable();
    });

    document.querySelectorAll("#results-table th[data-sort]").forEach((th) => {
      th.addEventListener("click", () => {
        const col = th.getAttribute("data-sort");
        if (state.tableFilter.sortCol === col) {
          state.tableFilter.sortAsc = !state.tableFilter.sortAsc;
        } else {
          state.tableFilter.sortCol = col;
          state.tableFilter.sortAsc = col === "BFDR" || col === "Bait" || col === "PreyGene";
        }
        renderResultsTable();
      });
    });

    $("btn-prev-page")?.addEventListener("click", () => {
      if (state.tableFilter.page > 1) {
        state.tableFilter.page--;
        renderResultsTable();
      }
    });
    $("btn-next-page")?.addEventListener("click", () => {
      state.tableFilter.page++;
      renderResultsTable();
    });

    $("btn-export-results-csv")?.addEventListener("click", () => {
      if (!state.pipelineResult) return;
      const includeBfdr = $("cb-export-bfdr")?.checked ?? true;
      const csv = Engine.formatOutputCSV(state.pipelineResult.rows, includeBfdr);
      Charts.downloadBlob(csv, "saintpy_results.csv", "text/csv");
    });

    $("btn-export-filtered-csv")?.addEventListener("click", () => {
      if (!state.pipelineResult) return;
      const includeBfdr = $("cb-export-bfdr")?.checked ?? true;
      const filtered = getFilteredAndSortedRows();
      const csv = Engine.formatOutputCSV(filtered, includeBfdr);
      Charts.downloadBlob(csv, "saintpy_filtered_results.csv", "text/csv");
    });

    // Responsive window resize handler for Canvases
    let resizeTimer = null;
    window.addEventListener("resize", () => {
      clearTimeout(resizeTimer);
      resizeTimer = setTimeout(() => {
        if (!state.pipelineResult) return;
        renderActiveVisualization();
        Charts.renderReplicateInspector(
          $("replicate-inspector-canvas"),
          state.selectedRow,
          state.modelOpts.quant_type
        );
      }, 150);
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
