#!/usr/bin/env python3
"""Build the single-file page from template.html, flow.js, piv.js and checks.json.

    python3 build.py                       -> index.html, a complete page for GitHub Pages
    python3 build.py --artifact OUT URL    -> page body only, with the home link pointing to URL
"""
import html
import json
import re
import sys
from pathlib import Path

here = Path(__file__).resolve().parent

# ---- search engines and link previews (full page only) ----
META = {
    'repo': 'bubble-particle-collision',
    'name': 'Bubble and particle collision',
    'title': 'Bubble-particle collision efficiency simulator | Hashim Ishaque',
    'headline': 'Bubble and particle collision: the flow around a rising bubble',
    'description': ('Navier-Stokes flow around a rising bubble, particle paths and collision efficiency for flotation, '
                    'and synthetic PIV of the same flow. A theoretical simulation.'),
    'image_alt': 'Streamlines around a 1 mm bubble with its wake, and the paths of 30 µm particles that hit or miss it',
    'keywords': ('flotation, bubble-particle collision, collision efficiency, Navier-Stokes, rising bubble, Maxey-Riley, '
                 'particle trajectories, Sutherland, Weber-Paddock, Yoon-Luttrell, PIV, cross-correlation'),
    'languages': ['JavaScript', 'Python'],
    'about': [
        ('Froth flotation', 'https://en.wikipedia.org/wiki/Froth_flotation'),
        ('Navier-Stokes equations', 'https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations'),
        ('Particle image velocimetry', 'https://en.wikipedia.org/wiki/Particle_image_velocimetry'),
        ('Stokes number', 'https://en.wikipedia.org/wiki/Stokes_number'),
    ],
}
SITE = 'https://hashimishaque.com'
PERSON = {
    '@type': 'Person', '@id': SITE + '/#person', 'name': 'Hashim Ishaque', 'url': SITE + '/',
    'sameAs': ['https://www.linkedin.com/in/hashimishaque/', 'https://github.com/hashimishaque-platform'],
}
WEBSITE = {'@type': 'WebSite', '@id': SITE + '/#website', 'url': SITE + '/', 'name': 'Hashim Ishaque'}
MIT = 'https://opensource.org/licenses/MIT'


def head(meta, extra=()):
    """Title, description, canonical address, icons, link preview and structured data (schema.org)."""
    url = f"{SITE}/{meta['repo']}/"
    image = url + 'social.png'
    page = {
        '@type': 'WebPage', '@id': url, 'url': url, 'name': meta['headline'], 'description': meta['description'],
        'inLanguage': 'en', 'isPartOf': {'@id': WEBSITE['@id']}, 'author': {'@id': PERSON['@id']},
        'primaryImageOfPage': {'@type': 'ImageObject', 'url': image, 'width': 1200, 'height': 630},
        'about': [{'@type': 'Thing', 'name': n, 'sameAs': u} for n, u in meta['about']],
        'keywords': meta['keywords'], 'breadcrumb': {'@id': url + '#breadcrumb'},
    }
    code = {
        '@type': 'SoftwareSourceCode', '@id': url + '#code', 'name': meta['name'], 'description': meta['description'],
        'url': url, 'codeRepository': f"https://github.com/hashimishaque-platform/{meta['repo']}",
        'programmingLanguage': meta['languages'], 'runtimePlatform': 'Web browser', 'license': MIT,
        'author': {'@id': PERSON['@id']},
    }
    crumbs = {'@type': 'BreadcrumbList', '@id': url + '#breadcrumb', 'itemListElement': [
        {'@type': 'ListItem', 'position': 1, 'name': 'Hashim Ishaque', 'item': SITE + '/'},
        {'@type': 'ListItem', 'position': 2, 'name': meta['name'], 'item': url}]}
    graph = [page, code, crumbs, WEBSITE, PERSON, *extra]
    ld = json.dumps({'@context': 'https://schema.org', '@graph': graph}, ensure_ascii=False, indent=1).replace('</', '<\\/')
    e = lambda s: html.escape(s, quote=True)
    return '\n'.join([
        f"<title>{e(meta['title'])}</title>",
        f'<meta name="description" content="{e(meta["description"])}">',
        '<meta name="author" content="Hashim Ishaque">',
        f'<link rel="canonical" href="{url}">',
        '<link rel="icon" href="/favicon.ico" sizes="32x32">',
        '<link rel="icon" href="/favicon.svg" type="image/svg+xml">',
        '<link rel="apple-touch-icon" href="/apple-touch-icon.png">',
        '<meta property="og:type" content="website">',
        '<meta property="og:site_name" content="Hashim Ishaque">',
        f'<meta property="og:title" content="{e(meta["headline"])}">',
        f'<meta property="og:description" content="{e(meta["description"])}">',
        f'<meta property="og:url" content="{url}">',
        f'<meta property="og:image" content="{image}">',
        '<meta property="og:image:width" content="1200">',
        '<meta property="og:image:height" content="630">',
        f'<meta property="og:image:alt" content="{e(meta["image_alt"])}">',
        '<meta name="twitter:card" content="summary_large_image">',
        f'<script type="application/ld+json">\n{ld}\n</script>',
    ])


EXTRA = ()


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
        page = re.sub(r'<title>.*?</title>\n', '', page, count=1)
        i = page.index('</style>') + len('</style>')
        page = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
                '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
                + head(META, EXTRA) + '\n'
                + page[:i].strip() + '\n</head>\n<body>\n' + page[i:].strip() + '\n</body>\n</html>\n')
        for ch in ('\u2013', '\u2014'):
            assert ch not in page, f'dash {ch!r} found'
    return page


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--artifact':
        out, home, full = Path(sys.argv[2]), sys.argv[3], False
    else:
        out, home, full = here / 'index.html', '/', True
    page = build(home, full)
    out.write_text(page, encoding='utf-8')
    print(f'wrote {out} ({len(page)} bytes)')
