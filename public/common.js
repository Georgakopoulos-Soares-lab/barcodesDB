function renderResults(container, data) {
  container.innerHTML = '';

  function compOf(seq) {
    const r = { A: 0, C: 0, G: 0, T: 0 };
    const s = String(seq || '').toUpperCase();
    for (let i = 0; i < s.length; i++) {
      const c = s[i];
      if (r[c] !== undefined) r[c]++;
    }
    return r;
  }
  function gcPctFromComp(comp, len) {
    if (!len) return 0;
    const gc = (comp.G || 0) + (comp.C || 0);
    return (gc * 100.0) / len;
  }
  function downloadText(filename, text) {
    const blob = new Blob([text], { type: 'text/plain;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }
  function toTSV(rows) {
    return rows.map(r => r.map(x => {
      const s = String(x ?? '');
      return s.includes('\t') || s.includes('\n') ? s.replace(/\t/g, ' ').replace(/\n/g, ' ') : s;
    }).join('\t')).join('\n') + '\n';
  }

  const header = document.createElement('div');
  header.className = 'row';

  const s1 = document.createElement('div'); s1.className='col card';
  s1.innerHTML = `<div class="stat"><strong>Total</strong><span class="pill">${data.total}</span></div>`;

  const s2 = document.createElement('div'); s2.className='col card';
  s2.innerHTML = `<div class="stat"><strong>Found</strong><span class="pill">${data.found}</span><span class="pill">${data.foundPct.toFixed(2)}%</span></div>`;

  container.appendChild(header);
  header.appendChild(s1);
  header.appendChild(s2);

  // Download block
  const tools = document.createElement('div');
  tools.className = 'card';
  tools.style.marginTop = '14px';
  tools.innerHTML = `
    <div style="display:flex;gap:12px;align-items:flex-end;flex-wrap:wrap;justify-content:space-between">
      <div>
        <div style="font-weight:800;margin-bottom:4px">Download</div>
        <div class="help">Export k-mers with composition (A/C/G/T) and GC%.</div>
      </div>
      <div style="display:flex;gap:10px;align-items:flex-end;flex-wrap:wrap">
        <div style="min-width:240px">
          <label class="help">Include</label>
          <select id="dlMode" class="input">
            <option value="present">Present only</option>
            <option value="absent">Not present only</option>
            <option value="both" selected>Both</option>
          </select>
        </div>
        <button class="btn secondary" id="dlBtn">Download TSV</button>
      </div>
    </div>
  `;
  container.appendChild(tools);

  const card = document.createElement('div');
  card.className = 'card';
  const table = document.createElement('table');
  table.className = 'table';
  table.innerHTML = `
    <thead>
      <tr>
        <th>#</th>
        <th>K-mer</th>
        <th>Present</th>
        <th>GC %</th>
        <th>A</th><th>C</th><th>G</th><th>T</th>
      </tr>
    </thead>
    <tbody></tbody>`;

  const tbody = table.querySelector('tbody');
  data.results.forEach((r, i) => {
    const tr = document.createElement('tr');
    tr.innerHTML = `
      <td>${i+1}</td>
      <td class="kmer">${r.kmer}</td>
      <td>${r.present ? '✓' : '—'}</td>
      <td>${r.gc.toFixed(2)}</td>
      <td>${r.comp.A||0}</td>
      <td>${r.comp.C||0}</td>
      <td>${r.comp.G||0}</td>
      <td>${r.comp.T||0}</td>`;
    tbody.appendChild(tr);
  });

  card.appendChild(table);
  container.appendChild(card);

  // wire download
  const dlMode = tools.querySelector('#dlMode');
  const dlBtn = tools.querySelector('#dlBtn');
  dlBtn.onclick = () => {
    const mode = dlMode.value;
    const filtered = data.results.filter(r => {
      if (mode === 'present') return !!r.present;
      if (mode === 'absent') return !r.present;
      return true;
    });
    const rows = [
      ['kmer', 'present', 'gc_pct', 'A', 'C', 'G', 'T'],
      ...filtered.map(r => {
        // Be robust if comp/gc were missing for any reason
        const comp = r.comp || compOf(r.kmer);
        const gc = (typeof r.gc === 'number') ? r.gc : gcPctFromComp(comp, String(r.kmer||'').length);
        return [r.kmer, r.present ? 1 : 0, gc.toFixed(2), comp.A||0, comp.C||0, comp.G||0, comp.T||0];
      })
    ];
    const stamp = new Date().toISOString().replace(/[:.]/g, '-');
    downloadText(`barcodes_kmer_lookup_${mode}_${stamp}.tsv`, toTSV(rows));
  };
}