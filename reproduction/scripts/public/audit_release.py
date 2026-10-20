"""Read-only release privacy guard; not a substitute for human review."""
from pathlib import Path
import re

root = Path(__file__).resolve().parents[3]
listed = set()
for p in root.rglob('*'):
    relative = p.relative_to(root)
    if '.git' in relative.parts:
        continue
    if p.is_symlink():
        raise SystemExit('Unsafe release symlink: ' + str(relative))
    if p.is_file():
        listed.add(str(relative))
patterns = [rb'-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----',
            rb'gh[pousr]_[A-Za-z0-9]{30,}', rb'github_pat_[A-Za-z0-9_]{40,}',
            rb'AKIA[0-9A-Z]{16}']
for name in sorted(listed):
    p = root / name
    if any(re.search(pattern,p.read_bytes()) for pattern in patterns):
        raise SystemExit('Potential credential in ' + name)
    if p.suffix.lower() in {'.rds','.rda','.rdata','.dta','.sav','.parquet','.sqlite','.pem','.key'}:
        raise SystemExit('Unexpected data/credential format: ' + name)
    if name.startswith('reproduction/data/') and name != 'reproduction/data/README.md':
        raise SystemExit('Unexpected participant input: ' + name)
print('Scanned',len(listed),'files; no selected secret signatures or participant input files found.')
print('Generated local results are not release files. Never use git add . after running analyses.')
