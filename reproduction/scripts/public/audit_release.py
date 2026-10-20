"""Read-only release integrity/privacy guard; not a substitute for human review."""
from pathlib import Path
import hashlib
import re
import sys

root = Path(__file__).resolve().parents[3]
expected = root / 'SHA256SUMS'
if not expected.is_file():
    raise SystemExit('SHA256SUMS is missing.')
listed = set()
for line in expected.read_text().splitlines():
    checksum, name = line.split('  ', 1)
    p = root / name
    if p.is_symlink() or not p.is_file() or root not in p.resolve().parents:
        raise SystemExit('Unsafe or missing release file: ' + name)
    if hashlib.sha256(p.read_bytes()).hexdigest() != checksum:
        raise SystemExit('Checksum mismatch: ' + name)
    listed.add(name)
actual = {str(p.relative_to(root)) for p in root.rglob('*')
          if p.is_file() and '.git' not in p.relative_to(root).parts}
extra = actual - listed - {'SHA256SUMS'}
if extra:
    raise SystemExit('Unreviewed files are present: ' + ', '.join(sorted(extra)[:10]))
patterns = [rb'-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----',
            rb'gh[pousr]_[A-Za-z0-9]{30,}', rb'github_pat_[A-Za-z0-9_]{40,}',
            rb'AKIA[0-9A-Z]{16}']
for name in listed:
    p = root / name
    if any(re.search(pattern,p.read_bytes()) for pattern in patterns):
        raise SystemExit('Potential credential in ' + name)
    if p.suffix.lower() in {'.rds','.rda','.rdata','.dta','.sav','.parquet','.sqlite','.pem','.key'}:
        raise SystemExit('Unexpected data/credential format: ' + name)
    if name.startswith('reproduction/data/') and name != 'reproduction/data/README.md':
        raise SystemExit('Unexpected participant input: ' + name)
print('Verified',len(listed),'release checksums; no selected secret signatures or participant input files.')
print('Generated local results are not release files. Never use git add . after running analyses.')
