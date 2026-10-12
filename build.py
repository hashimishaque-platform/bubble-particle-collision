#!/usr/bin/env python3
"""Build the single-file page from template.html, flow.js, piv.js and checks.json.

    python3 build.py                       -> index.html, a complete page for GitHub Pages
    python3 build.py --artifact OUT URL    -> page body only, with the home link pointing to URL
"""
import json
import re
import sys
from pathlib import Path

here = Path(__file__).resolve().parent


def module(name):
    src = (here / name).read_text(encoding='utf-8')
    # the CommonJS export block is only for the Node scripts
    stripped = re.sub(r"\nif \(typeof module !== 'undefined'\) \{.*?\n\}\n", "\n", src, flags=re.S)
    assert stripped != src, f'no export block found in {name}'
    return stripped.strip()


def build(home, full):
    page = (here / 'template.html').read_text(encoding='utf-8')
    checks = json.loads((here / 'checks.json').read_text(encoding='utf-8'))
    page = page.replace('/*FLOW*/', module('flow.js'))
    page = page.replace('/*PIV*/', module('piv.js'))
    page = page.replace('/*CHECKS*/', json.dumps(checks, ensure_ascii=False))
    if not home.startswith('http'):
        page = page.replace('href="/*HOME*/" target="_blank" rel="noopener"', f'href="{home}"')
    page = page.replace('/*HOME*/', home)
    left = re.findall(r'/\*[A-Z]+\*/', page)
    assert not left, f'placeholders left: {left}'
    for ch in ('\u2013', '\u2014'):
        assert ch not in page, f'dash {ch!r} found'
    if full:
        i = page.index('</style>') + len('</style>')
        page = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
                '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
                '<meta name="description" content="Navier-Stokes flow around a rising bubble, particle collision efficiency and synthetic PIV, computed from theory.">\n'
                + page[:i].strip() + '\n</head>\n<body>\n' + page[i:].strip() + '\n</body>\n</html>\n')
    return page


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--artifact':
        out, home, full = Path(sys.argv[2]), sys.argv[3], False
    else:
        out, home, full = here / 'index.html', '/', True
    page = build(home, full)
    out.write_text(page, encoding='utf-8')
    print(f'wrote {out} ({len(page)} bytes)')
