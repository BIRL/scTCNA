/* Plotly helpers shared by the pages: one colour per compartment everywhere, a common
 * layout, and the plot types the dashboards use (grouped bars, paired patient lines,
 * violins, the cross-cancer heatmap). Exposed on window.Plots. */
(function (global) {
  'use strict';

  const COLORS = {
    malignant: '#e4572e',
    non_malignant: '#1b998b',
    significant: '#1f5fa8',
    neutral: '#b8c2cc',
    accent: '#1f5fa8',
    line: '#9aa7b4',
  };
  const LABELS = { malignant: 'Malignant', non_malignant: 'Non-malignant' };
  const FONT = { family: 'system-ui, -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif', size: 12, color: '#22303c' };
  const GRID = '#eef1f5';

  function withAlpha(hex, alpha) {
    const n = parseInt(hex.slice(1), 16);
    return 'rgba(' + ((n >> 16) & 255) + ',' + ((n >> 8) & 255) + ',' + (n & 255) + ',' + alpha + ')';
  }

  function node(div) {
    return typeof div === 'string' ? document.getElementById(div) : div;
  }

  function layout(over) {
    over = over || {};
    const base = {
      font: FONT,
      paper_bgcolor: '#ffffff',
      plot_bgcolor: '#ffffff',
      margin: { l: 60, r: 16, t: 36, b: 56 },
      legend: { orientation: 'h', x: 0, y: 1.02, yanchor: 'bottom' },
      hovermode: 'closest',
      hoverlabel: { bgcolor: '#ffffff', bordercolor: '#c5cfd9', font: FONT },
    };
    const out = Object.assign({}, base, over);
    out.xaxis = Object.assign({ automargin: true, zeroline: false, gridcolor: GRID }, over.xaxis || {});
    out.yaxis = Object.assign({ automargin: true, zeroline: false, gridcolor: GRID }, over.yaxis || {});
    return out;
  }

  function config(filename) {
    return {
      responsive: true,
      displaylogo: false,
      modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d', 'toggleSpikelines'],
      toImageButtonOptions: { format: 'png', filename: filename || 'sctcna_plot', scale: 2 },
    };
  }

  function draw(div, data, layoutOver, filename) {
    const target = node(div);
    target.classList.remove('plot-empty');
    target.textContent = '';
    return Plotly.react(target, data, layout(layoutOver), config(filename));
  }

  /* Clear a plot and show a short message in its place. */
  function empty(div, message) {
    const target = node(div);
    if (target.data) Plotly.purge(target);
    target.classList.add('plot-empty');
    target.textContent = message || 'No data.';
  }

  function hoverNumbers(values, name) {
    return values.map((v) => name + ': ' + App.fmt(v));
  }

  /* Grouped bars, one trace per series: series = [{key, name, values, error, color, hover}]. */
  function groupedBars(div, categories, series, opts) {
    opts = opts || {};
    const data = series.map((s) => {
      const trace = {
        type: 'bar',
        name: s.name || LABELS[s.key] || s.key,
        x: categories,
        y: s.values,
        marker: { color: s.color || COLORS[s.key] || COLORS.accent },
        hovertext: s.hover || null,
        hoverinfo: s.hover ? 'text' : 'x+y+name',
      };
      if (s.error) trace.error_y = { type: 'data', array: s.error, visible: true, color: '#55606b', thickness: 1.2, width: 4 };
      return trace;
    });
    return draw(
      div,
      data,
      {
        barmode: 'group',
        bargap: 0.25,
        xaxis: { title: opts.xtitle || '', type: 'category', tickangle: categories.length > 8 ? -45 : 0 },
        yaxis: { title: opts.ytitle || '', rangemode: opts.rangemode || 'tozero' },
        showlegend: series.length > 1,
      },
      opts.filename
    );
  }

  /* Single bar series with an explicit colour per bar (e.g. significant vs not). */
  function bars(div, categories, values, opts) {
    opts = opts || {};
    const data = [
      {
        type: 'bar',
        x: opts.horizontal ? values : categories,
        y: opts.horizontal ? categories : values,
        orientation: opts.horizontal ? 'h' : 'v',
        marker: { color: opts.colors || COLORS.accent },
        text: opts.text || null,
        textposition: opts.text ? 'outside' : 'none',
        cliponaxis: false,
        hovertext: opts.hover || null,
        hoverinfo: opts.hover ? 'text' : 'x+y',
      },
    ];
    const lay = {
      showlegend: false,
      bargap: 0.3,
      xaxis: { title: opts.xtitle || '', zeroline: !!opts.horizontal, zerolinecolor: '#8a97a4' },
      yaxis: { title: opts.ytitle || '', zeroline: !opts.horizontal, zerolinecolor: '#8a97a4' },
    };
    if (opts.horizontal) {
      lay.yaxis.type = 'category';
      lay.yaxis.autorange = 'reversed';
      lay.margin = { l: 90, r: 24, t: 24, b: 48 };
    } else {
      lay.xaxis.type = 'category';
      lay.xaxis.tickangle = categories.length > 8 ? -45 : 0;
    }
    if (opts.height) lay.height = opts.height;
    return draw(div, data, lay, opts.filename);
  }

  /* One line per patient joining its Malignant and Non-malignant values (TCNA-style
   * paired view). pairs = [{label, a, b}] with a = malignant, b = non-malignant; an optional
   * second set (opts.overlay = {label, pairs}) is drawn dashed, e.g. the downsampled CV. */
  function pairedLines(div, pairs, opts) {
    opts = opts || {};
    const x = [LABELS.malignant, LABELS.non_malignant];
    const data = [];
    function lines(set, dashed, name) {
      set.forEach((p, i) => {
        data.push({
          type: 'scatter',
          mode: 'lines+markers',
          x: x,
          y: [p.a, p.b],
          name: name,
          legendgroup: name,
          showlegend: i === 0,
          line: { color: dashed ? '#c3cbd3' : COLORS.line, width: dashed ? 1.2 : 1.5, dash: dashed ? 'dot' : 'solid' },
          marker: { color: [COLORS.malignant, COLORS.non_malignant], size: dashed ? 7 : 9, symbol: dashed ? 'circle-open' : 'circle', line: { width: dashed ? 1.5 : 0, color: [COLORS.malignant, COLORS.non_malignant] } },
          hovertext: [p.label + ' (' + LABELS.malignant + '): ' + App.fmt(p.a), p.label + ' (' + LABELS.non_malignant + '): ' + App.fmt(p.b)],
          hoverinfo: 'text',
        });
      });
    }
    lines(pairs, false, opts.name || 'per patient');
    if (opts.overlay && opts.overlay.pairs && opts.overlay.pairs.length) lines(opts.overlay.pairs, true, opts.overlay.label || 'downsampled');
    // medians of each compartment as wide dashes
    const med = (vals) => {
      const v = vals.filter(App.isNum).sort((p, q) => p - q);
      if (!v.length) return null;
      const m = Math.floor(v.length / 2);
      return v.length % 2 ? v[m] : (v[m - 1] + v[m]) / 2;
    };
    const medA = med(pairs.map((p) => p.a));
    const medB = med(pairs.map((p) => p.b));
    if (medA !== null && medB !== null) {
      data.push({
        type: 'scatter',
        mode: 'markers',
        x: x,
        y: [medA, medB],
        name: 'median across patients',
        marker: { symbol: 'line-ew', size: 26, color: [COLORS.malignant, COLORS.non_malignant], line: { width: 3, color: [COLORS.malignant, COLORS.non_malignant] } },
        hovertext: ['Median ' + LABELS.malignant + ': ' + App.fmt(medA), 'Median ' + LABELS.non_malignant + ': ' + App.fmt(medB)],
        hoverinfo: 'text',
      });
    }
    return draw(
      div,
      data,
      {
        xaxis: { type: 'category', range: [-0.5, 1.5], tickfont: { size: 13 } },
        yaxis: { title: opts.ytitle || '', rangemode: opts.rangemode || 'normal' },
        showlegend: true,
      },
      opts.filename
    );
  }

  /* Violins of the per-cell values of each compartment with per-patient medians as dots. */
  function violins(div, groups, patients, opts) {
    opts = opts || {};
    const data = [];
    ['malignant', 'non_malignant'].forEach((key) => {
      const g = groups[key];
      if (!g || !g.values || !g.values.length) return;
      data.push({
        type: 'violin',
        name: LABELS[key],
        x0: LABELS[key],
        y: g.values,
        box: { visible: true, width: 0.15, fillcolor: '#ffffff', line: { color: '#55606b' } },
        meanline: { visible: false },
        points: false,
        spanmode: 'hard',
        line: { color: COLORS[key], width: 1.2 },
        fillcolor: withAlpha(COLORS[key], 0.35),
        hoverinfo: 'y+name',
        showlegend: false,
      });
    });
    const pts = (patients || []).filter((p) => App.isNum(p.median));
    if (pts.length) {
      data.push({
        type: 'scatter',
        mode: 'markers',
        name: 'patient medians',
        x: pts.map((p) => LABELS[p.compartment] || p.compartment),
        y: pts.map((p) => p.median),
        marker: { color: pts.map((p) => COLORS[p.compartment] || COLORS.accent), size: 9, line: { color: '#22303c', width: 1 } },
        hovertext: pts.map((p) => p.patient + ' (' + (LABELS[p.compartment] || p.compartment) + '): median ' + App.fmt(p.median) + ', ' + App.fmtInt(p.n_cells) + ' cells'),
        hoverinfo: 'text',
      });
    }
    return draw(
      div,
      data,
      {
        violinmode: 'group',
        xaxis: { type: 'category', categoryorder: 'array', categoryarray: [LABELS.malignant, LABELS.non_malignant], tickfont: { size: 13 } },
        yaxis: { title: opts.ytitle || '' },
        showlegend: pts.length > 0,
      },
      opts.filename
    );
  }

  /* z-score every row over its non-null entries (rows with < 2 values become all null). */
  function zscoreRows(values) {
    return values.map((row) => {
      const v = row.filter(App.isNum);
      if (v.length < 2) return row.map(() => null);
      const mean = v.reduce((a, b) => a + b, 0) / v.length;
      const sd = Math.sqrt(v.reduce((a, b) => a + (b - mean) * (b - mean), 0) / (v.length - 1));
      return row.map((x) => (App.isNum(x) ? (sd > 0 ? (x - mean) / sd : 0) : null));
    });
  }

  /* Genes x cancer types heatmap; nulls stay blank. opts: {diverging, zscore, qValues, title}. */
  function heatmap(div, genes, cancerTypes, values, opts) {
    opts = opts || {};
    const z = opts.zscore ? zscoreRows(values) : values;
    const diverging = opts.zscore || opts.diverging;
    let zmax = 0;
    z.forEach((row) => row.forEach((x) => {
      if (App.isNum(x)) zmax = Math.max(zmax, Math.abs(x));
    }));
    const text = z.map((row, i) =>
      row.map((x, j) => {
        let s = genes[i] + ' · ' + cancerTypes[j] + '<br>' + (opts.title || 'value') + ': ' + App.fmt(values[i][j]);
        if (opts.zscore) s += '<br>row z-score: ' + App.fmt(x);
        if (opts.qValues && opts.qValues[i]) s += '<br>q: ' + App.fmtP(opts.qValues[i][j]);
        return s;
      })
    );
    const trace = {
      type: 'heatmap',
      x: cancerTypes,
      y: genes,
      z: z,
      text: text,
      hoverinfo: 'text',
      hoverongaps: false,
      xgap: 2,
      ygap: 2,
      colorscale: diverging ? 'RdBu' : 'YlOrRd',
      reversescale: false, // RdBu as drawn by Plotly: high values red, low values blue
      colorbar: { title: { text: opts.zscore ? 'row z-score' : opts.colorbar || '', side: 'right' }, thickness: 14, len: 0.8 },
    };
    if (diverging) {
      trace.zmid = 0;
      trace.zmin = -zmax || -1;
      trace.zmax = zmax || 1;
    }
    const height = Math.max(340, 24 * genes.length + 140);
    return draw(
      div,
      [trace],
      {
        height: height,
        margin: { l: 90, r: 16, t: 24, b: 70 },
        xaxis: { type: 'category', side: 'bottom', tickangle: cancerTypes.length > 8 ? -45 : 0, showgrid: false },
        yaxis: { type: 'category', autorange: 'reversed', showgrid: false, tickfont: { size: genes.length > 60 ? 9 : 11 } },
      },
      opts.filename
    );
  }

  global.Plots = {
    COLORS: COLORS,
    LABELS: LABELS,
    withAlpha: withAlpha,
    draw: draw,
    empty: empty,
    groupedBars: groupedBars,
    bars: bars,
    pairedLines: pairedLines,
    violins: violins,
    heatmap: heatmap,
    zscoreRows: zscoreRows,
    hoverNumbers: hoverNumbers,
  };
})(window);
