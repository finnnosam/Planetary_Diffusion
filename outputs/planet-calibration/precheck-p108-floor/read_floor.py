"""Pass 108 (journal 145): apply rule-floor.md to a floor trace. python read_floor.py <trace.jsonl> [> floor.txt]"""
import json
import sys

import numpy as np

EARTH = {'holes': 13., 'crossings_per_mkm2': 19.14, 'lowflat4_of_land': 0.1424}
MS = ('holes', 'crossings_per_mkm2', 'lowflat4_of_land')
SHORT = {'holes': 'H', 'crossings_per_mkm2': 'C', 'lowflat4_of_land': 'L'}


def main(f):
    recs = [json.loads(l) for l in open(f) if l.strip()]
    st = lambda name, m: float(np.mean([r['stages'][name][m] for r in recs if name in r['stages']]))
    F0 = {m: st('end', m) for m in MS}
    e = {m: (F0[m]-EARTH[m]) if m != 'lowflat4_of_land' else (EARTH[m]-F0[m]) for m in MS}
    sgn = {m: (1. if m != 'lowflat4_of_land' else -1.) for m in MS}

    def share(D, base='end'):
        return {m: sgn[m]*(st(base, m)-st(D, m))/e[m] for m in MS}

    print(f'{len(recs)} steps, t {recs[0]["t"]:+.0f} ... {recs[-1]["t"]:+.0f}')
    print('means: reading | H | C | L (% of land) | islands | rough_z')
    for D in ('end', 'end_local', 'F_T', 'F_S', 'F_A', 'F_TS', 'F_TSA', 'A_T', 'A_S', 'A_A', 'A_TS', 'A_TSA',
              'returned', 'ret_same', 'ret_next', 'world_next', 'world_next_local'):
        if D in recs[-1]['stages']:
            print(f"{D} | {st(D, 'holes'):.1f} | {st(D, 'crossings_per_mkm2'):.2f} | {100*st(D, 'lowflat4_of_land'):.2f} | "
                  f"{st(D, 'islands'):.1f} | {st(D, 'rough_z'):.1f}")
    print(f"Earth | {EARTH['holes']:.0f} | {EARTH['crossings_per_mkm2']:.2f} | {100*EARTH['lowflat4_of_land']:.2f}")
    print(f"excess e: H {e['holes']:.1f}, C {e['crossings_per_mkm2']:.2f}, L {100*e['lowflat4_of_land']:.2f} points")

    def fmt(s):
        return ', '.join(f'{SHORT[m]} {s[m]:+.2f}' for m in MS)

    print('\nR1 support: s_A0', fmt(share('end_local')))
    named = []
    print('R2 material under flexure (s_F; named if >= 0.25 on two of three):')
    sF = {}
    for X in ('T', 'S', 'A', 'TS', 'TSA'):
        s = share(f'F_{X}')
        sF[X] = s
        n_ok = sum(s[m] >= .25 for m in MS)
        n_neg = sum(s[m] < -.10 for m in MS)
        tag = ''
        if X in ('T', 'S', 'A') and n_ok >= 2:
            named.append(X); tag = '  NAMED'
        if X == 'TSA' and sum(s[m] >= .5 for m in MS) >= 2:
            tag = '  MAIN SOURCE (>= 0.50 on two of three)'
        if n_neg >= 2:
            tag += '  (rougher: compensates)'
        print(f'  F_{X}: {fmt(s)}{tag}')
    print('R3 channel (s\'_A from A0 on the same excess; flexure channel if s_F >= 2 s\'_A):')
    for X in ('T', 'S', 'A', 'TS', 'TSA'):
        s2 = {m: sgn[m]*(st('end_local', m)-st(f'A_{X}', m))/e[m] for m in MS}
        ch = [SHORT[m] for m in MS if sF[X][m] >= 2*s2[m] and sF[X][m] > 0]
        print(f"  A_{X}: {fmt(s2)}; flexure channel on {', '.join(ch) if ch else 'none'}")
    print('R4 resampling (jumps J = M(D) - M(end), H and C):')
    J = {D: {m: st(D, m)-F0[m] for m in MS[:2]} for D in ('ret_same', 'ret_next', 'world_next', 'world_next_local', 'returned')}
    for D, v in J.items():
        print(f"  J_{D}: H {v['holes']:+.1f}, C {v['crossings_per_mkm2']:+.2f}")
    ratio = {m: J['world_next'][m]/J['ret_next'][m] if J['ret_next'][m] else float('nan') for m in MS[:2]}
    verdict = ('world resampling makes most' if all(r >= .5 for r in ratio.values()) else
               'the rasters\' departure makes it' if all(r < .25 for r in ratio.values()) else 'both')
    print(f"  J_world_next / J_ret_next: H {ratio['holes']:.2f}, C {ratio['crossings_per_mkm2']:.2f} -> {verdict}")
    # validity: ret_next at step k against composed_olddyn at step k+1
    a = [r['stages']['ret_next'] for r in recs[:-1]]
    b = [r['stages']['composed_olddyn'] for r in recs[1:]]
    vh = np.mean([x['holes'] for x in a])/np.mean([x['holes'] for x in b])-1
    vc = np.mean([x['crossings_per_mkm2'] for x in a])/np.mean([x['crossings_per_mkm2'] for x in b])-1
    print(f"  validity: ret_next(k) vs composed_olddyn(k+1): holes {100*vh:+.1f}% (+-10), crossings {100*vc:+.1f}% (+-5) -> "
          f"{'valid' if abs(vh) <= .10 and abs(vc) <= .05 else 'NOT VALID'}")
    tk = ('world_next_fallback_share', 'world_next_fallback_cont_share', 'next_owner_change_share', 'unread_cont_share')
    print('  transfer:', ', '.join(f"{k} {np.mean([r['transfer'][k] for r in recs if k in r['transfer']]):.4f}" for k in tk))
    print(f"\nnamed by R2: {named or 'none'}")


if __name__ == '__main__':
    main(sys.argv[1])
