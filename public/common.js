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

  // Hamming distance summary card (only when check was applied)
  const hammingApplied = data.hamming_distance_check_applied;
  if (hammingApplied) {
    const hdr = document.createElement('div');
    hdr.className = 'card';
    hdr.style.marginTop = '14px';
    hdr.innerHTML = `
      <div style="font-weight:800;margin-bottom:8px">Hamming Distance Check</div>
      <div class="row" style="gap:16px">
        <div class="col" style="flex:1"><span class="pill">Min requested: ${data.min_hamming_distance_requested}</span></div>
        <div class="col" style="flex:1"><span class="pill">Nearest observed: ${data.nearest_hamming_distance_observed !== null && data.nearest_hamming_distance_observed !== undefined ? data.nearest_hamming_distance_observed : 'N/A'}</span></div>
        <div class="col" style="flex:1"><span class="pill" style="${data.passes_min_hamming_distance ? 'background:#1b5e20;color:#fff' : 'background:#b71c1c;color:#fff'}">${data.passes_min_hamming_distance ? '✓ Passes' : '✗ Fails'}</span></div>
      </div>
    `;
    container.appendChild(hdr);
  }

  // Motif filter summary card
  const motifApplied = data.motif_filter_applied;
  if (motifApplied) {
    const mhdr = document.createElement('div');
    mhdr.className = 'card';
    mhdr.style.marginTop = '14px';
    let optsHtml = '';
    const mo = data.motif_options || {};
    if (mo.filter_homopolymers) optsHtml += `<span class="pill">Homopolymers (max ${mo.max_homopolymer})</span>`;
    if (mo.filter_low_complexity) optsHtml += `<span class="pill">Low complexity (min entropy ${mo.min_shannon_entropy})</span>`;
    if (mo.filter_dinucleotide_repeats) optsHtml += `<span class="pill">Dinucleotide repeats</span>`;
    if (mo.filter_trinucleotide_repeats) optsHtml += `<span class="pill">Trinucleotide repeats</span>`;
    if (mo.filter_tetranucleotide_repeats) optsHtml += `<span class="pill">Tetranucleotide repeats</span>`;
    if (mo.filter_restriction_sites) optsHtml += `<span class="pill">Restriction sites</span>`;
    if (mo.filter_functional_motifs) optsHtml += `<span class="pill">Functional motifs</span>`;
    mhdr.innerHTML = `
      <div style="font-weight:800;margin-bottom:8px">Sequence Motif Filters</div>
      <div class="row" style="gap:16px">
        <div class="col" style="flex:1"><span class="pill">Mode: ${data.motif_mode || 'flag'}</span></div>
        <div class="col" style="flex:1"><span class="pill">Failing: ${data.motif_fail_count ?? 0}</span></div>
      </div>
      <div style="margin-top:8px;display:flex;gap:6px;flex-wrap:wrap">${optsHtml}</div>
      <div class="help" style="margin-top:6px">These filters only annotate the queried sequence. They do not alter the exact existence lookup result.</div>
    `;
    container.appendChild(mhdr);
  }

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

  const hasHamming = hammingApplied;
  const hasMotif = motifApplied;
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
        ${hasHamming ? '<th>Nearest Hamming</th><th>Passes Min</th>' : ''}
        ${hasMotif ? '<th>Motif passes</th><th>Motif hits</th>' : ''}
      </tr>
    </thead>
    <tbody></tbody>`;

  const tbody = table.querySelector('tbody');
  data.results.forEach((r, i) => {
    const tr = document.createElement('tr');
    let hammingCells = '';
    if (hasHamming) {
      const nearest = r.nearest_hamming_distance_observed !== undefined ? r.nearest_hamming_distance_observed : '—';
      const passes = r.passes_min_hamming_distance !== undefined ? (r.passes_min_hamming_distance ? '✓' : '✗') : '—';
      hammingCells = `<td>${nearest}</td><td>${passes}</td>`;
    }
    let motifCells = '';
    if (hasMotif) {
      const passes = r.motif_passes !== undefined ? (r.motif_passes ? '✓' : '✗') : '—';
      const hits = r.motif_hits || r.motif_hit_count ? (r.motif_hits || r.motif_hit_count) : '—';
      motifCells = `<td>${passes}</td><td style="max-width:300px;word-break:break-all;font-size:0.85em">${hits}</td>`;
    }
    tr.innerHTML = `
      <td>${i+1}</td>
      <td class="kmer">${r.kmer}</td>
      <td>${r.present ? '✓' : '—'}</td>
      <td>${r.gc.toFixed(2)}</td>
      <td>${r.comp.A||0}</td>
      <td>${r.comp.C||0}</td>
      <td>${r.comp.G||0}</td>
      <td>${r.comp.T||0}</td>
      ${hammingCells}
      ${motifCells}`;
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
    const headers = ['kmer', 'present', 'gc_pct', 'A', 'C', 'G', 'T'];
    if (hasHamming) headers.push('nearest_hamming', 'passes_min_hamming');
    if (hasMotif) headers.push('motif_passes', 'motif_hits');
    const rows = [
      headers,
      ...filtered.map(r => {
        // Be robust if comp/gc were missing for any reason
        const comp = r.comp || compOf(r.kmer);
        const gc = (typeof r.gc === 'number') ? r.gc : gcPctFromComp(comp, String(r.kmer||'').length);
        const row = [r.kmer, r.present ? 1 : 0, gc.toFixed(2), comp.A||0, comp.C||0, comp.G||0, comp.T||0];
        if (hasHamming) {
          row.push(r.nearest_hamming_distance_observed !== undefined ? r.nearest_hamming_distance_observed : '');
          row.push(r.passes_min_hamming_distance !== undefined ? (r.passes_min_hamming_distance ? 1 : 0) : '');
        }
        if (hasMotif) {
          row.push(r.motif_passes !== undefined ? (r.motif_passes ? 1 : 0) : '');
          row.push(r.motif_hits || '');
        }
        return row;
      })
    ];
    const stamp = new Date().toISOString().replace(/[:.]/g, '-');
    downloadText(`barcodes_kmer_lookup_${mode}_${stamp}.tsv`, toTSV(rows));
  };
}