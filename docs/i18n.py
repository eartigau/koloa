#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The French version of koloa's site.

make_page.py writes index.html in English; french() turns it into
index_fr.html with a catalogue of translations (docs/i18n/fr.json). Every
block of text of the page (a paragraph, an item of a list, a heading, a
caption, a cell of a table, a link of the sidebar, the title of a page) is
looked up by its English text, with its numbers and its equations replaced
by placeholders ({0}, {1}... and {m0}, {m1}...): the catalogue holds the
words, and the numbers of the demos, whatever they are, go back in. Code,
the output of the examples, the equations and the text inside the figures
stay as they are. A block that the catalogue lacks stays in English and is
listed in docs/i18n/missing_fr.json, to be translated. The links to
Wikipedia go to the French articles where there is one
(docs/i18n/wikipedia_fr.json, from Wikipedia's interlanguage links).

    python i18n.py extract     # every block of index.html, to translate
                               # (docs/i18n/source.json)
    python i18n.py prune       # drop the translations of blocks the page
                               # no longer has

Created on 2026-09-29

@author: artigau
"""
import json
import os
import re
import sys
import warnings
from html import unescape

from bs4 import BeautifulSoup, MarkupResemblesLocatorWarning, NavigableString, Tag

# a short translation (a word, a name) looks like a file name to
#   BeautifulSoup, which warns about it
warnings.filterwarnings('ignore', category=MarkupResemblesLocatorWarning)

HERE = os.path.dirname(os.path.abspath(__file__))
I18N = os.path.join(HERE, 'i18n')
#: the elements whose text is translated as a whole (their inline tags,
#: links, emphasis, code, go with them)
BLOCKS = ('title, h1, h2, h3, h4, p, li, figcaption, td, th, '
          '.uc-row > span, .uc-crit, .uc-measured, .uc-results, '
          '.ex-label > span, .nav-title, .home-link')
#: never translated: code, the output of the examples, scripts, styles
SKIP = ('script', 'style', 'pre', 'code')
SKIP_CLASSES = {'codeblock', 'ex-output', 'katex', 'mono'}
#: the equations of KaTeX, and the numbers (a decimal one may be followed
#: by a letter that is not ASCII, as in 0.7σ; an integer may not, so that
#: 2σ stays a word)
MATH = re.compile(r'\\\(.*?\\\)|\\\[.*?\\\]', re.S)
NUMBER = re.compile(r'(?<![\w.#&;/-])[-+]?(?:\d+\.\d+(?:[eE][-+]?\d+)?'
                    r'(?![A-Za-z0-9_])|\d+(?:[eE][-+]?\d+)?(?![\w]))')
#: a date (2020-07-03, with its time 04:59 or not) is one placeholder
DATE = re.compile(r'(?<![\w.-])\d{4}-\d{2}-\d{2}(?: \d{2}:\d{2})?(?![\w])')
NUMBER_OR_DATE = re.compile(DATE.pattern + '|' + NUMBER.pattern)
#: the badges of the tests
BADGES = {'PASS': 'RÉUSSI', 'FAIL': 'ÉCHEC', 'INFO': 'INFO'}


# =============================================================================
# Keys
# =============================================================================
def _skipped(tag: Tag) -> bool:
    """an element never translated, or inside one"""
    for node in [tag] + list(tag.parents):
        if not isinstance(node, Tag):
            continue
        if node.name in SKIP:
            return True
        if SKIP_CLASSES & set(node.get('class') or []):
            return True
    return False


def key_of(html: str):
    """
    The key of a block: its inner HTML, whitespace collapsed, the equations
    and the numbers of its text (never those of its tags) replaced by
    placeholders

    :param html: str, the inner HTML of the block

    :return: tuple, the key, the numbers and the equations it took out
    """
    html = ' '.join(html.split())
    maths, numbers = [], []

    def take_math(mt):
        maths.append(mt.group(0))
        return '{m%d}' % (len(maths) - 1)
    html = MATH.sub(take_math, html)
    parts = re.split(r'(<[^>]+>)', html)
    for it, part in enumerate(parts):
        if part.startswith('<'):
            continue

        def take_number(mt):
            numbers.append(mt.group(0))
            return '{%d}' % (len(numbers) - 1)
        parts[it] = NUMBER_OR_DATE.sub(take_number, part)
    return ''.join(parts).strip(), numbers, maths


def fill_key(text: str, numbers, maths) -> str:
    """a translation with the numbers and equations of the English block,
    in one pass (an equation may hold braces such as \\tfrac{1+k}{2}, which
    a second pass would take for a placeholder)"""
    return re.sub(r'\{(m?)(\d+)\}',
                  lambda mt: (maths if mt.group(1) else numbers)[
                      int(mt.group(2))], text)


def blocks(soup: BeautifulSoup):
    """
    The blocks of the page to translate: the elements of BLOCKS that hold
    no other one (a list item with a list in it gives its own link instead),
    outside code and equations

    :return: list of Tag
    """
    chosen = [tag for tag in soup.select(BLOCKS) if not _skipped(tag)]
    ids = {id(tag) for tag in chosen}
    out = []
    for tag in chosen:
        inner = [sub for sub in tag.find_all(True) if id(sub) in ids]
        if not inner:
            out.append(tag)
        elif tag.name == 'li':
            # a list item that holds a list: its link (the sidebar)
            out += [sub for sub in tag.find_all('a', recursive=False)]
    return out


def attributes(soup: BeautifulSoup):
    """the attributes that are text: the titles of the pages (the sidebar
    and the pager read them), the descriptions of the figures and links,
    and the description of the site"""
    for tag in soup.find_all(True):
        classes = tag.get('class') or []
        for name in ('data-title', 'alt', 'aria-label', 'title'):
            if not tag.has_attr(name) or _skipped(tag) \
                    or 'lang-link' in classes:
                continue
            # the pages of the modules are named after the modules
            if name == 'data-title' and 'ex-panel' in classes:
                continue
            yield tag, name
    for tag in soup.find_all('meta', attrs={'name': 'description'}):
        yield tag, 'content'


# =============================================================================
# Extract and translate
# =============================================================================
def extract(html: str) -> list:
    """every key of the page, in order, once"""
    soup = BeautifulSoup(html, 'html.parser')
    keys = []
    for tag in blocks(soup):
        key = key_of(tag.decode_contents())[0]
        if re.search(r'[A-Za-z]', re.sub(r'<[^>]+>|\{m?\d+\}', '', key)):
            keys.append(key)
    for tag, name in attributes(soup):
        key = key_of(tag[name])[0]
        if re.search(r'[A-Za-z]', key):
            keys.append(key)
    return list(dict.fromkeys(keys))


def load_catalogue(lang: str = 'fr') -> dict:
    """the catalogue of a language, or an empty one"""
    path = os.path.join(I18N, f'{lang}.json')
    if not os.path.exists(path):
        return {}
    with open(path) as handle:
        return json.load(handle)


def translate(html: str, lang: str = 'fr'):
    """
    The page in another language, from its catalogue

    :param html: str, the page in English
    :param lang: str, the language (its catalogue is docs/i18n/<lang>.json)

    :return: tuple, the page and the keys the catalogue lacks
    """
    cat = load_catalogue(lang)
    soup = BeautifulSoup(html, 'html.parser')
    missing = []
    for tag in blocks(soup):
        key, numbers, maths = key_of(tag.decode_contents())
        if not re.search(r'[A-Za-z]', re.sub(r'<[^>]+>|\{m?\d+\}', '', key)):
            continue
        if key not in cat:
            missing.append(key)
            continue
        new = BeautifulSoup(fill_key(cat[key], numbers, maths), 'html.parser')
        tag.clear()
        for node in list(new.contents):
            tag.append(node)
    for tag, name in attributes(soup):
        key, numbers, maths = key_of(tag[name])
        if not re.search(r'[A-Za-z]', key):
            continue
        if key in cat:
            # an attribute holds text, not HTML: its entities decoded
            tag[name] = unescape(fill_key(cat[key], numbers, maths))
        else:
            missing.append(key)
    # the badges of the tests
    for tag in soup.select('.uc-badge'):
        for word, new in BADGES.items():
            for node in tag.find_all(string=re.compile(word)):
                node.replace_with(NavigableString(str(node).replace(word, new)))
    soup.html['lang'] = lang
    # the definitions in the language of the page, where Wikipedia has one
    #   (docs/i18n/wikipedia_<lang>.json, from its interlanguage links)
    path = os.path.join(I18N, f'wikipedia_{lang}.json')
    if os.path.exists(path):
        with open(path) as handle:
            wiki = json.load(handle)
        for tag in soup.find_all('a', href=True):
            if tag['href'] in wiki:
                tag['href'] = wiki[tag['href']]
    # the screenshots of the GUI in the language of the page, where there
    #   are (figures/gui/<lang>/<name>, the GUI in that language)
    for tag in soup.find_all('img', src=True):
        src = tag['src']
        if src.startswith('figures/gui/') and src.count('/') == 2:
            local = src.replace('figures/gui/', f'figures/gui/{lang}/')
            if os.path.exists(os.path.join(HERE, local)):
                tag['src'] = local
    # the language link: back to the English page
    for tag in soup.select('a.lang-link'):
        tag['href'] = 'index.html'
        tag['hreflang'] = tag['lang'] = 'en'
        tag['title'] = 'English version'
        tag.string = 'EN'
    return str(soup), list(dict.fromkeys(missing))


def french(html: str) -> str:
    """
    index_fr.html from index.html: the translated page, its language link
    back to English, and the list of what the catalogue lacks

    :param html: str, the English page

    :return: str, the French page
    """
    page, missing = translate(html, 'fr')
    os.makedirs(I18N, exist_ok=True)
    with open(os.path.join(I18N, 'missing_fr.json'), 'w') as handle:
        json.dump(missing, handle, indent=1, ensure_ascii=False)
    return page, missing


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'extract':
        with open(os.path.join(HERE, 'index.html')) as handle:
            keys = extract(handle.read())
        os.makedirs(I18N, exist_ok=True)
        with open(os.path.join(I18N, 'source.json'), 'w') as handle:
            json.dump(keys, handle, indent=1, ensure_ascii=False)
        print(f'{len(keys)} blocks written to i18n/source.json')
    elif len(sys.argv) > 1 and sys.argv[1] == 'prune':
        with open(os.path.join(HERE, 'index.html')) as handle:
            keys = set(extract(handle.read()))
        catalogue = load_catalogue('fr')
        kept = {key: val for key, val in catalogue.items() if key in keys}
        with open(os.path.join(I18N, 'fr.json'), 'w') as handle:
            json.dump(kept, handle, indent=1, ensure_ascii=False)
        print(f'{len(catalogue) - len(kept)} translations the page no longer '
              f'uses removed; {len(kept)} kept')

# =============================================================================
# End of code
# =============================================================================
