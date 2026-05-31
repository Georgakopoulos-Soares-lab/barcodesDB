// barcodesDB server
// Express API wrapping local binaries: query_kmer_bitmap and query_substring_bitmap_stream (sharded)

const express = require('express');
const multer = require('multer');
const { spawn } = require('child_process');
const path = require('path');

const fs = require('fs');

const app = express();
const uploadsDir = path.join(__dirname, 'uploads');
if (!fs.existsSync(uploadsDir)) fs.mkdirSync(uploadsDir, { recursive: true });
const upload = multer({ dest: uploadsDir });

const ROOT = path.resolve(__dirname, '..');

// Binaries
const BIN_QUERY_KMER = path.join(ROOT, 'query_kmer_bitmap');
const BIN_QUERY_SUBSTR = path.join(ROOT, 'query_substring_bitmap_stream'); // sharded substring binary (updated)

// Data
const BITMAP_16 = path.join(ROOT, 'roar_barcodes_16.bin');
const BITMAP_17 = path.join(ROOT, 'roar_barcodes_17.bin');
const BITMAP_18 = path.join(ROOT, 'roar_barcodes_18.bin');

const SHARDS_16 = path.join(ROOT, 'shards_16');
const SHARDS_17 = path.join(ROOT, 'shards_17');
const SHARDS_18 = path.join(ROOT, 'shards_18');

const GC_HIST_16 = path.join(ROOT, 'gc_hist_shards_16.json');
const GC_HIST_17 = path.join(ROOT, 'gc_hist_shards_17.json');
const GC_HIST_18 = path.join(ROOT, 'gc_hist_shards_18.json'); // per-shard GC histogram

function getShardsForK(k) {
  if (k === 16) return { shards: SHARDS_16, gcHist: GC_HIST_16, baseK: 16 };
  if (k === 17) return { shards: SHARDS_17, gcHist: GC_HIST_17, baseK: 17 };
  // k==18 or expansion-from-18
  return { shards: SHARDS_18, gcHist: GC_HIST_18, baseK: 18 };
}

/**
 * Parse and sanitize motif filter options from a request body.
 * All filters default to off. Returns a clean options object.
 */
function parseMotifOptions(body) {
  const opts = {
    motif_mode: 'off',
    filter_homopolymers: false,
    max_homopolymer: 4,
    filter_low_complexity: false,
    min_shannon_entropy: 1.5,
    filter_dinucleotide_repeats: false,
    filter_trinucleotide_repeats: false,
    filter_tetranucleotide_repeats: false,
    filter_restriction_sites: false,
    filter_functional_motifs: false,
  };

  if (!body) return opts;

  // motif_mode: must be one of "off", "flag", "exclude"
  const modeRaw = String(body.motif_mode || 'off').trim().toLowerCase();
  if (['off', 'flag', 'exclude'].includes(modeRaw)) {
    opts.motif_mode = modeRaw;
  }

  // Booleans: accept true/false, "true"/"false", "1"/"0"
  function toBool(v) {
    if (v === true || v === 'true' || v === '1' || v === 1) return true;
    return false;
  }

  opts.filter_homopolymers = toBool(body.filter_homopolymers);
  opts.filter_low_complexity = toBool(body.filter_low_complexity);
  opts.filter_dinucleotide_repeats = toBool(body.filter_dinucleotide_repeats);
  opts.filter_trinucleotide_repeats = toBool(body.filter_trinucleotide_repeats);
  opts.filter_tetranucleotide_repeats = toBool(body.filter_tetranucleotide_repeats);
  opts.filter_restriction_sites = toBool(body.filter_restriction_sites);
  opts.filter_functional_motifs = toBool(body.filter_functional_motifs);

  // max_homopolymer: default 4, minimum 1
  if (body.max_homopolymer !== undefined && body.max_homopolymer !== null && body.max_homopolymer !== '') {
    const v = parseInt(String(body.max_homopolymer), 10);
    if (!isNaN(v) && v >= 1) opts.max_homopolymer = v;
  }

  // min_shannon_entropy: default 1.5, clamp between 0 and 2
  if (body.min_shannon_entropy !== undefined && body.min_shannon_entropy !== null && body.min_shannon_entropy !== '') {
    const v = parseFloat(String(body.min_shannon_entropy));
    if (!isNaN(v)) opts.min_shannon_entropy = Math.max(0, Math.min(2, v));
  }

  return opts;
}

/**
 * Append CLI args for motif filters to the given args array.
 * Only adds args when they differ from defaults or are enabled.
 */
function appendMotifArgs(args, opts) {
  if (!opts || opts.motif_mode === 'off') return;

  args.push('--motif-mode', opts.motif_mode);

  if (opts.filter_homopolymers) {
    args.push('--filter-homopolymers');
    if (opts.max_homopolymer !== 4) {
      args.push('--max-homopolymer', String(opts.max_homopolymer));
    }
  }
  if (opts.filter_low_complexity) {
    args.push('--filter-low-complexity');
    if (Math.abs(opts.min_shannon_entropy - 1.5) > 0.001) {
      args.push('--min-shannon-entropy', String(opts.min_shannon_entropy));
    }
  }
  if (opts.filter_dinucleotide_repeats) {
    args.push('--filter-dinucleotide-repeats');
  }
  if (opts.filter_trinucleotide_repeats) {
    args.push('--filter-trinucleotide-repeats');
  }
  if (opts.filter_tetranucleotide_repeats) {
    args.push('--filter-tetranucleotide-repeats');
  }
  if (opts.filter_restriction_sites) {
    args.push('--filter-restriction-sites');
  }
  if (opts.filter_functional_motifs) {
    args.push('--filter-functional-motifs');
  }
}

function isDNA(s) {
  return /^[ACGTacgt]+$/.test(s);
}

function gcContent(seq) {
  if (!seq || seq.length === 0) return 0;
  let gc = 0;
  for (const c of seq.toUpperCase()) if (c === 'G' || c === 'C') gc++;
  return (gc * 100.0) / seq.length;
}

// ── Motif filtering helpers for the random-generate endpoint ──────────────────
const GEN_RE_SITES = [
  ['EcoRI','GAATTC'],['BamHI','GGATCC'],['HindIII','AAGCTT'],['NotI','GCGGCCGC'],
  ['XhoI','CTCGAG'],['XbaI','TCTAGA'],['SpeI','ACTAGT'],['NheI','GCTAGC'],
  ['PstI','CTGCAG'],['KpnI','GGTACC'],['SacI','GAGCTC'],['SalI','GTCGAC'],
  ['SmaI','CCCGGG'],['MluI','ACGCGT'],['AgeI','ACCGGT'],['BglII','AGATCT'],
  ['AatII','GACGTC'],['AccI','GTATAC'],['AflII','CTTAAG'],['AscI','GGCGCGCC'],
  ['AvrII','CCTAGG'],['BclI','TGATCA'],['BsiWI','CGTACG'],['BspHI','TCATGA'],
  ['BsrGI','TGTACA'],['BssHII','GCGCGC'],['ClaI','ATCGAT'],['DraI','TTTAAA'],
  ['EcoRV','GATATC'],['FseI','GGCCGGCC'],['HpaI','GTTAAC'],['MfeI','CAATTG'],
  ['NcoI','CCATGG'],['NdeI','CATATG'],['NruI','TCGCGA'],['PacI','TTAATTAA'],
  ['PmeI','GTTTAAAC'],['PvuII','CAGCTG'],['SbfI','CCTGCAGG'],['ScaI','AGTACT'],
  ['SnaBI','TACGTA'],['SspI','AATATT'],['StuI','AGGCCT'],['SwaI','ATTTAAAT'],
  ['AluI','AGCT'],['HaeIII','GGCC'],['MseI','TTAA'],['MspI','CCGG'],['RsaI','GTAC'],['TaqI','TCGA']
];
const GEN_FUNC_MOTIFS = [
  ['TATA_box','TATAAA'],['GC_box','GGGCGG'],['CCAAT_box','CCAAT'],
  ['CRE','TGACGTCA'],['E_box_canonical','CACGTG'],['AP1_site','TGACTCA'],
  ['NFkB_site','GGGACTTTCC'],['STAT3_site','TTCCGGGAA'],['OCT4_site','ATGCAAAT'],
  ['E2F_site','TTTCCCGC'],['Ets_binding_core','CCGGAAGT'],
  ['polyA_signal','AATAAA'],['alt_polyA_signal','ATTAAA'],
  ['polyA_variant_AATACA','AATACA'],['polyA_variant_GATAAA','GATAAA'],
  ['polyA_variant_AATAGA','AATAGA'],['polyA_variant_AATGAA','AATGAA'],
  ['ARE_mRNA_stability','ATTTAT'],['splice_donor_like','GTAGT'],
  ['splice_acceptor_like','CAGG'],['Kozak_like_core','ACCATGG'],
  ['Shine_Dalgarno','AGGAGG'],['Pribnow_box','TATAAT'],
  ['minus35_sigma70','TTGACA'],['CpG','CG']
];

function revComp(seq) {
  const m = { A:'T', T:'A', C:'G', G:'C' };
  return seq.split('').reverse().map(b => m[b] || b).join('');
}

// Returns { passes: bool, hits: string } — hits is pipe-joined descriptor list
function evaluateMotifFilter(kmer, opts) {
  const up = kmer.toUpperCase();
  const rc = revComp(up);
  const hitList = [];

  if (opts.filter_homopolymers) {
    const maxRun = opts.max_homopolymer || 4;
    let run = 1, runStart = 0;
    for (let i = 1; i <= up.length; i++) {
      if (i < up.length && up[i] === up[i-1]) { run++; }
      else {
        if (run > maxRun) hitList.push(`homopolymer:${up.substring(runStart, i)}@${runStart}`);
        run = 1; runStart = i;
      }
    }
  }

  if (opts.filter_low_complexity) {
    const minH = opts.min_shannon_entropy || 1.5;
    const cnt = { A:0, C:0, G:0, T:0 };
    for (const b of up) cnt[b] = (cnt[b] || 0) + 1;
    let H = 0;
    for (const c of Object.values(cnt)) {
      if (c > 0) { const p = c / up.length; H -= p * Math.log2(p); }
    }
    if (H < minH) hitList.push(`low_complexity:H=${H.toFixed(2)}`);
  }

  if (opts.filter_dinucleotide_repeats) {
    for (let i = 0; i + 6 <= up.length; i++) {
      const di = up.substring(i, i + 2);
      if (di[0] !== di[1]) {
        let j = i + 2;
        while (j + 2 <= up.length && up.substring(j, j + 2) === di) j += 2;
        if (j - i >= 6) { hitList.push(`dinucl_repeat:${di}@${i}`); i = j - 1; }
      }
    }
  }

  if (opts.filter_trinucleotide_repeats) {
    for (let i = 0; i + 9 <= up.length; i++) {
      const tri = up.substring(i, i + 3);
      let j = i + 3;
      while (j + 3 <= up.length && up.substring(j, j + 3) === tri) j += 3;
      if (j - i >= 9) { hitList.push(`trinucl_repeat:${tri}@${i}`); i = j - 1; }
    }
  }

  if (opts.filter_tetranucleotide_repeats) {
    for (let i = 0; i + 12 <= up.length; i++) {
      const tetra = up.substring(i, i + 4);
      let j = i + 4;
      while (j + 4 <= up.length && up.substring(j, j + 4) === tetra) j += 4;
      if (j - i >= 12) { hitList.push(`tetranucl_repeat:${tetra}@${i}`); i = j - 1; }
    }
  }

  if (opts.filter_restriction_sites) {
    const sites = (opts.restriction_site_list && opts.restriction_site_list.length > 0)
      ? opts.restriction_site_list.map(s => ['custom', s])
      : GEN_RE_SITES;
    for (const [name, seq] of sites) {
      let idx = up.indexOf(seq);
      if (idx >= 0) hitList.push(`restriction_site:${name}:${seq}@${idx}`);
      else { idx = rc.indexOf(seq); if (idx >= 0) hitList.push(`restriction_site:${name}:${seq}@rc${idx}`); }
    }
  }

  if (opts.filter_functional_motifs) {
    const motifs = (opts.functional_motif_list && opts.functional_motif_list.length > 0)
      ? opts.functional_motif_list.map(s => ['custom', s])
      : GEN_FUNC_MOTIFS;
    for (const [name, seq] of motifs) {
      let idx = up.indexOf(seq);
      if (idx >= 0) hitList.push(`functional_motif:${name}:${seq}@${idx}`);
      else { idx = rc.indexOf(seq); if (idx >= 0) hitList.push(`functional_motif:${name}:${seq}@rc${idx}`); }
    }
  }

  return { passes: hitList.length === 0, hits: hitList.join('|') };
}

function passesMotifFilter(kmer, opts) {
  return evaluateMotifFilter(kmer, opts).passes;
}

function ntComp(seq) {
  const r = { A: 0, C: 0, G: 0, T: 0 };
  for (const c of seq.toUpperCase()) if (r[c] !== undefined) r[c]++;
  return r;
}

function runBinary(cmd, args, { stdinData, timeoutMs }) {
  return new Promise((resolve, reject) => {
    const p = spawn(cmd, args, { stdio: ['pipe', 'pipe', 'pipe'] });
    let stdout = '';
    let stderr = '';
    let timedOut = false;

    const to = timeoutMs ? setTimeout(() => {
      timedOut = true;
      p.kill('SIGKILL');
    }, timeoutMs) : null;

    p.stdout.on('data', (d) => (stdout += d.toString()));
    p.stderr.on('data', (d) => (stderr += d.toString()));

    p.on('error', (err) => {
      if (to) clearTimeout(to);
      reject(err);
    });

    p.on('close', (code) => {
      if (to) clearTimeout(to);
      if (timedOut) return reject(new Error('Process timed out'));
      if (code !== 0) {
        const e = new Error(`Process exited ${code}: ${stderr}`);
        e.code = code;
        e.stderr = stderr;
        e.stdout = stdout;
        return reject(e);
      }
      resolve({ stdout, stderr });
    });

    if (stdinData) p.stdin.write(stdinData);
    p.stdin.end();
  });
}

function parseSubstringStdout(stdout) {
  const lines = stdout.split(/\r?\n/).map((s) => s.trim()).filter(Boolean);
  if (!lines.length) {
    return { nextCursor: '', hasMore: false, returned: 0, kOut: null, kmers: [], motifMeta: null };
  }

  const meta = lines[0];
  if (!meta.startsWith('__META__')) {
    return { nextCursor: '', hasMore: false, returned: lines.length, kOut: null, kmers: lines, motifMeta: null };
  }

  const parts = meta.split('\t');
  const nextCursor = (parts[1] ?? '').trim();
  const hasMore = (parts[2] ?? '0').trim() === '1';
  const returned = Number(parts[3] ?? lines.length - 1) || (lines.length - 1);
  const kOut = Number(parts[4] ?? '') || null;

  // Parse extended __META__ fields (key=value pairs)
  const motifMeta = {};
  for (let i = 5; i < parts.length; i++) {
    const kv = parts[i].trim();
    const eq = kv.indexOf('=');
    if (eq > 0) {
      const key = kv.substring(0, eq);
      const val = kv.substring(eq + 1);
      motifMeta[key] = val;
    }
  }

  const kmers = lines.slice(1);
  return { nextCursor, hasMore, returned, kOut, kmers, motifMeta };
}

// app.use(express.static(path.join(__dirname, 'public')));
app.use('/', express.static(path.join(__dirname, 'public')));


app.use(express.json({ limit: '5mb' }));

// /api/query-kmer
app.post('/api/query-kmer', upload.single('kmersFile'), async (req, res) => {
  try {
    let kmers = [];
    if (req.file) {
      const text = fs.readFileSync(req.file.path, 'utf8');
      kmers = text.split(/\r?\n/).map((s) => s.trim()).filter((s) => s);
      fs.unlink(req.file.path, () => {});
    } else if (Array.isArray(req.body.kmers)) {
      kmers = req.body.kmers.map((s) => String(s).trim()).filter((s) => s);
    }
    if (!kmers.length) return res.status(400).json({ error: 'No kmers provided' });

    const uniq = Array.from(new Set(kmers));
    for (const k of uniq) {
      if (!isDNA(k)) return res.status(400).json({ error: `Invalid k-mer: ${k}` });
    }

    // Enforce: this endpoint accepts a single k per request (no mixed k-mer lengths).
    const ks = Array.from(new Set(uniq.map((s) => s.length)));
    if (ks.length !== 1) {
      return res.status(400).json({
        error: `Mixed k-mer lengths are not supported. Found lengths: ${ks.sort((a,b)=>a-b).join(', ')}`,
      });
    }
    const kReq = ks[0];
    if (![16, 17, 18].includes(kReq)) {
      return res.status(400).json({ error: `Only k=16, k=17, and k=18 are supported (got k=${kReq})` });
    }

    // Parse optional min_hamming_distance
    let minHamming = 0;
    if (req.body.min_hamming_distance !== undefined && req.body.min_hamming_distance !== null && req.body.min_hamming_distance !== '') {
      minHamming = parseInt(String(req.body.min_hamming_distance), 10);
      if (isNaN(minHamming) || minHamming < 0) minHamming = 0;
      if (minHamming > kReq) minHamming = kReq;
    }

    // Parse optional motif filter parameters
    const motifOpts = parseMotifOptions(req.body);

    const tmpFile = path.join(__dirname, 'uploads', `kmers_${Date.now()}.txt`);
    fs.writeFileSync(tmpFile, uniq.join('\n'));

    const { shards: shardsDir } = getShardsForK(kReq);
    const args = ['--shards', shardsDir, '--k', String(kReq), '--kmers', tmpFile];
    if (minHamming > 0) {
      args.push('--min-hamming-distance', String(minHamming));
    }
    // Add motif filter args
    appendMotifArgs(args, motifOpts);

    const { stdout } = await runBinary(BIN_QUERY_KMER, args, { timeoutMs: 120000 });
    fs.unlink(tmpFile, () => {});

    const results = [];
    let foundCount = 0;
    let hammingCheckApplied = (minHamming > 0);
    let overallPassesMinHamming = true;
    let overallNearestHamming = Infinity;
    let motifApplied = (motifOpts.motif_mode !== 'off');
    let motifFailCount = 0;
    for (const line of stdout.split(/\r?\n/)) {
      if (!line) continue;
      const parts = line.split('\t');
      const kmer = parts[0];
      const hit = parts[1];
      const present = hit === '1';
      if (present) foundCount++;

      const result = { kmer, present, gc: gcContent(kmer), comp: ntComp(kmer) };

      if (hammingCheckApplied && parts.length >= 4) {
        const nearestStr = parts[2];
        const passesStr = parts[3];
        const nearest = nearestStr ? parseInt(nearestStr, 10) : -1;
        const passes = passesStr === '1';
        if (nearest >= 0) {
          result.nearest_hamming_distance_observed = nearest;
          result.passes_min_hamming_distance = passes;
          if (nearest < overallNearestHamming) overallNearestHamming = nearest;
          if (!passes) overallPassesMinHamming = false;
        }
      }

      // Parse motif metadata (appended after hamming columns if both present)
      const motifOffset = hammingCheckApplied ? 4 : 2; // after kmer + present + optional hamming cols
      if (motifApplied && parts.length >= motifOffset + 3) {
        const motifPasses = parts[motifOffset] === '1';
        const motifHitCount = parseInt(parts[motifOffset + 1], 10) || 0;
        const motifHitsStr = parts[motifOffset + 2] || '';
        result.motif_passes = motifPasses;
        result.motif_hit_count = motifHitCount;
        result.motif_hits = motifHitsStr;
        if (!motifPasses) motifFailCount++;
      }

      results.push(result);
    }

    const response = {
      total: results.length,
      found: foundCount,
      foundPct: results.length ? (foundCount * 100.0) / results.length : 0,
      results,
    };

    if (hammingCheckApplied) {
      response.min_hamming_distance_requested = minHamming;
      response.hamming_distance_check_applied = true;
      response.passes_min_hamming_distance = overallPassesMinHamming;
      response.nearest_hamming_distance_observed = overallNearestHamming < Infinity ? overallNearestHamming : null;
    }

    if (motifApplied) {
      response.motif_filter_applied = true;
      response.motif_mode = motifOpts.motif_mode;
      response.motif_fail_count = motifFailCount;
      response.motif_options = {
        filter_homopolymers: motifOpts.filter_homopolymers || false,
        max_homopolymer: motifOpts.max_homopolymer,
        filter_low_complexity: motifOpts.filter_low_complexity || false,
        min_shannon_entropy: motifOpts.min_shannon_entropy,
        filter_dinucleotide_repeats: motifOpts.filter_dinucleotide_repeats || false,
        filter_restriction_sites: motifOpts.filter_restriction_sites || false,
        filter_functional_motifs: motifOpts.filter_functional_motifs || false,
      };
    }

    res.json(response);
  } catch (err) {
    console.error(err);
    res.status(500).json({ error: String(err.message || err) });
  }
});

// /api/query-substring (backend filters: substring optional + gc range required + optional constructK)
app.post('/api/query-substring', async (req, res) => {
  try {
    const body = req.body || {};
    const substringRaw = (typeof body.substring === 'string') ? body.substring.trim() : '';
    const substring = substringRaw; // optional

    const gcMin = Number.isFinite(Number(body.gcMin)) ? Math.floor(Number(body.gcMin)) : 0;
    const gcMax = Number.isFinite(Number(body.gcMax)) ? Math.floor(Number(body.gcMax)) : 100;
    if (gcMin < 0 || gcMax > 100 || gcMin > gcMax) {
      return res.status(400).json({ error: 'gcMin/gcMax must satisfy 0 <= gcMin <= gcMax <= 100' });
    }

    // ConstructK: optional; if empty/null => defaults to base-k (see below)
    const constructKRaw = (body.constructK === null || body.constructK === undefined) ? '' : String(body.constructK).trim();
    const constructK = constructKRaw ? Math.floor(Number(constructKRaw)) : null;
    if (constructKRaw && (!Number.isFinite(constructK) || constructK < 16 || constructK > 32)) {
      return res.status(400).json({ error: 'constructK must be an integer between 16 and 32' });
    }

    if (substring && !isDNA(substring)) {
      return res.status(400).json({ error: 'substring must be A/C/G/T only' });
    }
    const kOut = constructK ?? 18;
    const maxSubLen = kOut;
    if (substring && substring.length > maxSubLen) {
      return res.status(400).json({ error: `substring length must be <= ${maxSubLen}` });
    }

    const pageSize = Number.isFinite(Number(body.limit)) ? Math.floor(Number(body.limit)) : 200;
    if (pageSize < 1 || pageSize > 50000) {
      return res.status(400).json({ error: 'limit must be between 1 and 50000' });
    }

    // Optional: threads for backend; clamp to avoid overload
    const threadsReq = Number.isFinite(Number(body.threads)) ? Math.floor(Number(body.threads)) : 16;
    const threads = Math.max(1, Math.min(64, threadsReq));

    const cursorUsed = (typeof body.cursor === 'string' && body.cursor.trim()) ? body.cursor.trim() : '';

    // Parse optional motif filter parameters
    const motifOpts = parseMotifOptions(body);

    // Decide shard base.
    // Rules:
    // - For requested kOut in {16,17,18} => use that exact shard set.
    // - For kOut > 18 => expansion is only supported from 18-mer shards.
    const baseK = (kOut > 18) ? 18 : kOut;
    if (![16, 17, 18].includes(baseK)) {
      return res.status(400).json({ error: 'Only k=16,17,18 are supported as base lengths' });
    }
    const { shards: shardsDir, gcHist: gcHist } = getShardsForK(baseK);

    if (!fs.existsSync(gcHist)) {
      return res.status(500).json({ error: `GC histogram not found: ${gcHist}` });
    }

    const args = [
      '--shards', shardsDir,
      '--gc-hist', gcHist,
      '--limit', String(pageSize),
      '--threads', String(threads),
      '--gc-min', String(gcMin),
      '--gc-max', String(gcMax),
      '--random_access',
    ];
    if (kOut) args.push('--construct_k', String(kOut));
    if (substring) args.push('--substring', substring);
    if (cursorUsed) args.push('--cursor', cursorUsed);
    if (body.reverse_complement) args.push('--reverse_complement');

    // Add motif filter args
    appendMotifArgs(args, motifOpts);

    const { stdout } = await runBinary(BIN_QUERY_SUBSTR, args, { timeoutMs: 2 * 60 * 1000 });
    const parsed = parseSubstringStdout(stdout);

    const motifApplied = !!(motifOpts.motif_mode && motifOpts.motif_mode !== 'off');
    let motifFailCount = 0;
    const results = parsed.kmers.map((rawLine) => {
      const parts = rawLine.split('\t');
      const kmer = parts[0];
      const result = { kmer, gc: gcContent(kmer), comp: ntComp(kmer) };
      if (motifApplied && parts.length >= 3) {
        const motifPasses = parts[1] === '1';
        result.motif_passes = motifPasses;
        result.motif_hit_count = parseInt(parts[2], 10) || 0;
        result.motif_hits = parts[3] || '';
        if (!motifPasses) motifFailCount++;
      }
      return result;
    });

    res.json({
      cursorUsed,
      nextCursor: parsed.nextCursor || '',
      hasMore: !!parsed.hasMore,
      returned: parsed.returned ?? results.length,
      kOut: parsed.kOut ?? kOut,

      // echo backend filter state
      substring,
      gcMin,
      gcMax,
      threads,
      constructK: kOut,
      baseK,

      // top-level motif filter state (used by UI for hasMotif check)
      motif_filter_applied: motifApplied,
      motif_mode: motifApplied ? motifOpts.motif_mode : undefined,
      motif_fail_count: motifApplied ? motifFailCount : undefined,

      // motif filter metadata (legacy compat)
      motifMeta: parsed.motifMeta && Object.keys(parsed.motifMeta).length > 0 ? {
        motif_filter_applied: parsed.motifMeta.motif_filter_applied === '1',
        motif_mode: parsed.motifMeta.motif_mode || 'off',
        motif_filtered_count: parseInt(parsed.motifMeta.motif_filtered_count, 10) || 0,
        returned_count: parseInt(parsed.motifMeta.returned_count, 10) || results.length,
      } : undefined,

      results,
    });
  } catch (err) {
    console.error(err);
    res.status(500).json({ error: String(err.message || err) });
  }
});

// /api/generate-random-kmers — Generate random DNA k-mers matching criteria
// This does NOT query the bitmap. Users can generate barcodes independently.
app.post('/api/generate-random-kmers', async (req, res) => {
  try {
    const body = req.body || {};
    const k = parseInt(String(body.k || 18), 10);
    const n = parseInt(String(body.n ?? 10), 10);
    const gcMin = parseInt(String(body.gc_min || 0), 10);
    const gcMax = parseInt(String(body.gc_max || 100), 10);
    const substring = (typeof body.substring === 'string') ? body.substring.trim().toUpperCase() : '';

    // Validate
    if (isNaN(k) || k < 4 || k > 32) {
      return res.status(400).json({ error: 'k must be an integer between 4 and 32' });
    }
    if (isNaN(n) || n < 1 || n > 10000) {
      return res.status(400).json({ error: 'n must be between 1 and 10000' });
    }
    if (gcMin < 0 || gcMax > 100 || gcMin > gcMax) {
      return res.status(400).json({ error: 'gc_min/gc_max must satisfy 0 <= gc_min <= gc_max <= 100' });
    }
    if (substring && !isDNA(substring)) {
      return res.status(400).json({ error: 'substring must be A/C/G/T only' });
    }
    if (substring && substring.length > k) {
      return res.status(400).json({ error: 'substring length must be <= k' });
    }

    const bases = ['A', 'C', 'G', 'T'];
    const results = [];
    const MAX_ATTEMPTS = n * 500;  // more attempts since motif filter may reject many
    let attempts = 0;

    const motifOpts = parseMotifOptions(body);
    const motifActive = !!(motifOpts.motif_mode && motifOpts.motif_mode !== 'off');

    while (results.length < n && attempts < MAX_ATTEMPTS) {
      attempts++;
      // Generate a random k-mer
      let kmer = '';
      for (let i = 0; i < k; i++) {
        kmer += bases[Math.floor(Math.random() * 4)];
      }
      // Check GC%
      const gc = gcContent(kmer);
      if (gc < gcMin || gc > gcMax) continue;
      // Check substring (if specified)
      if (substring && !kmer.includes(substring)) continue;
      // Check motif filter (exclude mode rejects; flag mode accepts all but annotates)
      if (motifActive && motifOpts.motif_mode === 'exclude' && !passesMotifFilter(kmer, motifOpts)) continue;
      results.push(kmer);
    }

    res.json({
      k,
      n_requested: n,
      n_returned: results.length,
      gc_min: gcMin,
      gc_max: gcMax,
      substring: substring || null,
      attempts,
      motif_filter_applied: motifActive,
      motif_mode: motifActive ? motifOpts.motif_mode : undefined,
      results: results.map(kmer => {
        const r = { kmer, gc: gcContent(kmer), comp: ntComp(kmer) };
        if (motifActive && motifOpts.motif_mode === 'flag') {
          const ev = evaluateMotifFilter(kmer, motifOpts);
          r.motif_passes = ev.passes;
          r.motif_hits = ev.hits;
        }
        return r;
      }),
    });
  } catch (err) {
    console.error(err);
    res.status(500).json({ error: String(err.message || err) });
  }
});

// Pages
app.get('/', (req, res) => res.sendFile(path.join(__dirname, 'public', 'index.html')));
app.get('/kmer', (req, res) => res.sendFile(path.join(__dirname, 'public', 'kmer.html')));
app.get('/substring', (req, res) => res.sendFile(path.join(__dirname, 'public', 'substring.html')));
app.get('/generate', (req, res) => res.sendFile(path.join(__dirname, 'public', 'generate.html')));

const PORT = process.env.PORT || 8090;
const server = app.listen(PORT, () => {
  console.log(`barcodesDB listening on http://localhost:${PORT}`);
});
server.setTimeout(5 * 60 * 1000);