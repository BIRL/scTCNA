/* About page: fills the data-source table and the study citations from /api/datasets. */
(function () {
  'use strict';
  const $ = App.$;
  const el = App.el;

  async function init() {
    let rows;
    try {
      rows = await App.apiGet('/datasets');
    } catch (e) {
      return;
    }
    const container = $('datasets-table');
    container.innerHTML = '';
    if (!rows.length) {
      container.appendChild(el('div', { class: 'placeholder', text: 'No datasets loaded.' }));
      return;
    }
    $('datasets-note').textContent = rows.length + (rows.length === 1 ? ' dataset' : ' datasets') + '; cell counts are after QC and doublet removal.';
    const head = ['Dataset', 'Cancer', 'Study', 'Source', 'Technology', 'Site', 'Patients (paired)', 'Malignant cells', 'Non-malignant cells', 'Genes (expressed)', 'Downsampling depth', 'Download date', 'Data'];
    const table = el('table', { class: 'data-table' }, [el('thead', {}, [el('tr', {}, head.map((h) => el('th', { text: h })))])]);
    const tbody = el('tbody');
    rows.forEach((r) => {
      const cells = [
        el('td', { text: r.dataset }),
        el('td', { class: 'text', text: r.cancer_type + (r.cancer_name ? ' – ' + r.cancer_name : '') + (r.tcna_site ? ' (TCNA site: ' + r.tcna_site + ')' : '') }),
        el('td', { class: 'text', text: r.study || '–' }),
        el('td', { class: 'text', text: r.source || '–' }),
        el('td', { class: 'text', text: r.technology || '–' }),
        el('td', { class: 'text', text: r.site || '–' }),
        el('td', { text: App.fmtInt(r.n_patients) + ' (' + App.fmtInt(r.n_patients_paired) + ')' }),
        el('td', { text: App.fmtInt(r.n_cells_malignant) }),
        el('td', { text: App.fmtInt(r.n_cells_non_malignant) }),
        el('td', { text: App.fmtInt(r.n_genes) + ' (' + App.fmtInt(r.n_genes_expressed) + ')' }),
        el('td', { text: App.fmtInt(r.downsampling_target_depth) + ' UMIs' }),
        el('td', { class: 'text', text: r.download_date ? String(r.download_date).slice(0, 10) : '–' }),
        el('td', { class: 'text' }, r.data_url ? [el('a', { href: r.data_url, target: '_blank', rel: 'noopener', text: 'bundle' })] : ['–']),
      ];
      tbody.appendChild(el('tr', {}, cells));
    });
    table.appendChild(tbody);
    container.appendChild(el('div', { class: 'table-wrap' }, [table]));

    const refs = $('datasets-citations');
    refs.innerHTML = '';
    rows.forEach((r) => {
      refs.appendChild(el('li', { text: (r.cancer_type ? r.cancer_type + ': ' : '') + (r.citation || r.study || r.dataset) + (r.notes ? ' – ' + r.notes : '') }));
    });
  }

  document.addEventListener('DOMContentLoaded', init);
})();
