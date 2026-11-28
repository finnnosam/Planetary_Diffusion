"""Compare two ident.py fingerprints: the first differing step and keys, or IDENTICAL."""
import sys, json
a, b = (json.load(open(f)) for f in sys.argv[1:3])
ra, rb = a['rows'], b['rows']
for i, (x, y) in enumerate(zip(ra, rb)):
    bad = sorted(k for k in set(x) | set(y) if x.get(k) != y.get(k))
    if bad:
        print(f'DIFFERENT at step {i}: {bad[:12]}'); sys.exit(1)
if len(ra) != len(rb):
    print(f'DIFFERENT lengths {len(ra)} {len(rb)}'); sys.exit(1)
print(f'IDENTICAL over {len(ra)} steps ({a["seconds"]:.1f} s vs {b["seconds"]:.1f} s)')
