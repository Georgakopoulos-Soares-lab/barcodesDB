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
const BITMAP_18 = path.join(ROOT, 'roar_barcodes_18.bin');
const SHARDS_18 = path.join(ROOT, 'shards_18');
const GC_HIST_18 = path.join(ROOT, 'gc_hist_shards_18.json'); // per-shard GC histogram

function isDNA(s) {
  return /^[ACGTacgt]+$/.test(s);
}

function gcContent(seq) {
  if (!seq || seq.length === 0) return 0;
  let gc = 0;
  for (const c of seq.toUpperCase()) if (c === 'G' || c === 'C') gc++;
  return (gc * 100.0) / seq.length;
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
    return { nextCursor: '', hasMore: false, returned: 0, kOut: null, kmers: [] };
  }

  const meta = lines[0];
  if (!meta.startsWith('__META__')) {
    return { nextCursor: '', hasMore: false, returned: lines.length, kOut: null, kmers: lines };
  }

  const parts = meta.split('\t');
  const nextCursor = (parts[1] ?? '').trim();
  const hasMore = (parts[2] ?? '0').trim() === '1';
  const returned = Number(parts[3] ?? lines.length - 1) || (lines.length - 1);
  const kOut = Number(parts[4] ?? '') || null;

  const kmers = lines.slice(1);
  return { nextCursor, hasMore, returned, kOut, kmers };
}

// app.use(express.static(path.join(__dirname, 'public')));
app.use('/barcodesdb', express.static(path.join(__dirname, 'public')));

app.use(express.json({ limit: '5mb' }));

// /api/query-kmer (unchanged)
app.post('/barcodesdb/api/query-kmer', upload.single('kmersFile'), async (req, res) => {
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

    const tmpFile = path.join(__dirname, 'uploads', `kmers_${Date.now()}.txt`);
    fs.writeFileSync(tmpFile, uniq.join('\n'));

    const args = ['--bitmap', BITMAP_18, '--kmers', tmpFile];
    const { stdout } = await runBinary(BIN_QUERY_KMER, args, { timeoutMs: 120000 });
    fs.unlink(tmpFile, () => {});

    const results = [];
    let foundCount = 0;
    for (const line of stdout.split(/\r?\n/)) {
      if (!line) continue;
      const [kmer, hit] = line.split('\t');
      const present = hit === '1';
      if (present) foundCount++;
      results.push({ kmer, present, gc: gcContent(kmer), comp: ntComp(kmer) });
    }

    res.json({
      total: results.length,
      found: foundCount,
      foundPct: results.length ? (foundCount * 100.0) / results.length : 0,
      results,
    });
  } catch (err) {
    console.error(err);
    res.status(500).json({ error: String(err.message || err) });
  }
});

// /api/query-substring (backend filters: substring optional + gc range required + optional constructK)
app.post('/barcodesdb/api/query-substring', async (req, res) => {
  try {
    const body = req.body || {};
    const substringRaw = (typeof body.substring === 'string') ? body.substring.trim() : '';
    const substring = substringRaw; // optional

    const gcMin = Number.isFinite(Number(body.gcMin)) ? Math.floor(Number(body.gcMin)) : 0;
    const gcMax = Number.isFinite(Number(body.gcMax)) ? Math.floor(Number(body.gcMax)) : 100;
    if (gcMin < 0 || gcMax > 100 || gcMin > gcMax) {
      return res.status(400).json({ error: 'gcMin/gcMax must satisfy 0 <= gcMin <= gcMax <= 100' });
    }

    // ConstructK: optional; if empty/null => default (k=18 in binary)
    const constructKRaw = (body.constructK === null || body.constructK === undefined) ? '' : String(body.constructK).trim();
    const constructK = constructKRaw ? Math.floor(Number(constructKRaw)) : null;
    if (constructKRaw && (!Number.isFinite(constructK) || constructK < 18 || constructK > 32)) {
      return res.status(400).json({ error: 'constructK must be an integer between 18 and 32' });
    }

    if (substring && !isDNA(substring)) {
      return res.status(400).json({ error: 'substring must be A/C/G/T only' });
    }
    const maxSubLen = constructK || 18;
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

    if (!fs.existsSync(GC_HIST_18)) {
      return res.status(500).json({ error: `GC histogram not found: ${GC_HIST_18}` });
    }

    const args = [
      '--shards', SHARDS_18,
      '--gc-hist', GC_HIST_18,
      '--limit', String(pageSize),
      '--threads', String(threads),
      '--gc-min', String(gcMin),
      '--gc-max', String(gcMax),
      '--random_access',
    ];
    if (constructK) args.push('--construct_k', String(constructK));
    if (substring) args.push('--substring', substring);
    if (cursorUsed) args.push('--cursor', cursorUsed);
    if (body.reverse_complement) args.push('--reverse_complement');

    const { stdout } = await runBinary(BIN_QUERY_SUBSTR, args, { timeoutMs: 2 * 60 * 1000 });
    const parsed = parseSubstringStdout(stdout);

    const results = parsed.kmers.map((kmer) => ({
      kmer,
      gc: gcContent(kmer),
      comp: ntComp(kmer),
    }));

    res.json({
      cursorUsed,
      nextCursor: parsed.nextCursor || '',
      hasMore: !!parsed.hasMore,
      returned: parsed.returned ?? results.length,
      kOut: parsed.kOut ?? constructK ?? 18,

      // echo backend filter state
      substring,
      gcMin,
      gcMax,
      threads,
      constructK: constructK ?? 18,

      results,
    });
  } catch (err) {
    console.error(err);
    res.status(500).json({ error: String(err.message || err) });
  }
});

// Pages
app.get('/barcodesdb', (req, res) => res.sendFile(path.join(__dirname, 'public', 'index.html')));
app.get('/barcodesdb/kmer', (req, res) => res.sendFile(path.join(__dirname, 'public', 'kmer.html')));
app.get('/barcodesdb/substring', (req, res) => res.sendFile(path.join(__dirname, 'public', 'substring.html')));

const PORT = process.env.PORT || 8090;
const server = app.listen(PORT, () => {
  console.log(`barcodesDB listening on http://localhost:${PORT}`);
});
server.setTimeout(5 * 60 * 1000);