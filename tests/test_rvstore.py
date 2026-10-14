#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The team's private repository of velocities (koloa.rvstore), with a
repository of GitHub played here: what is read into a folder, what is
written, a commit another one got ahead of, a star never to share

Created on 2026-10-14

@author: artigau
"""
import base64
import hashlib
import json
import os

from koloa import rvstore, shared


class Repo:
    """GitHub's API of a repository, as gh answers it"""

    def __init__(self):
        self.blobs, self.commits = {}, {}
        self.head = self._commit({}, None)
        self.late = 0

    def _commit(self, files, parent):
        sha = hashlib.sha1(repr((sorted(files.items()), parent)).encode()
                           ).hexdigest()
        self.commits[sha] = dict(files)
        return sha

    def __call__(self, *args, body=None, timeout=120.0):
        path = args[0].split(f'repos/{rvstore.REPO}', 1)[-1]
        method = args[args.index('--method') + 1] if '--method' in args \
            else 'GET'
        if path == '':
            return 0, dict(permissions=dict(pull=True, push=True))
        if path == f'/git/ref/heads/{rvstore.BRANCH}':
            return 0, dict(object=dict(sha=self.head))
        if path.startswith('/git/trees/') and method == 'GET':
            files = self.commits[path.split('/')[3].split('?')[0]]
            return 0, dict(tree=[dict(path=key, sha=val, type='blob')
                                 for key, val in files.items()])
        if path == '/git/blobs' and method == 'POST':
            data = base64.b64decode(body['content'])
            sha = hashlib.sha1(b'blob %d\0' % len(data) + data).hexdigest()
            self.blobs[sha] = data
            return 0, dict(sha=sha)
        if path.startswith('/git/blobs/'):
            return 0, dict(content=base64.b64encode(
                self.blobs[path.split('/')[-1]]).decode())
        if path.startswith('/git/commits/') and method == 'GET':
            return 0, dict(tree=dict(sha='tree:' + path.split('/')[-1]))
        if path == '/git/trees' and method == 'POST':
            base = self.commits[body['base_tree'][5:]]
            files = dict(base, **{one['path']: one['sha']
                                  for one in body['tree']})
            sha = 'tree' + hashlib.sha1(repr(sorted(files.items())).encode()
                                        ).hexdigest()
            self.commits[sha] = files
            return 0, dict(sha=sha)
        if path == '/git/commits' and method == 'POST':
            sha = self._commit(self.commits[body['tree']],
                               body['parents'][0])
            self.parent = body['parents'][0]
            return 0, dict(sha=sha)
        if path == f'/git/refs/heads/{rvstore.BRANCH}':
            if self.late:
                # another commit got there first
                self.late -= 1
                files = dict(self.commits[self.head], other='x')
                self.head = self._commit(files, self.head)
                return 1, None
            if self.parent != self.head:
                return 1, None
            self.head = body['sha']
            return 0, dict()
        raise AssertionError(path)


def _star(folder, name='GJ 1'):
    for rel, text in (('target.json', json.dumps(dict(main=name))),
                      ('manifest.json', json.dumps(dict(folder='/my/disk',
                                                        archives={}))),
                      ('rv/dace/raw_GJ1.csv', 'rjd,rv\n1,2\n'),
                      ('lbl/GJ_1_HARPSN_lbl.rdb', 'rjd\tvrad\n--\t--\n1\t2\n'),
                      ('phot/tess.csv', 'big')):
        path = os.path.join(folder, *rel.split('/'))
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as handle:
            handle.write(text)


def test_written_and_read(tmp_path, monkeypatch):
    repo = Repo()
    monkeypatch.setattr(rvstore, '_gh', repo)
    monkeypatch.setattr(rvstore, '_ACCESS', {})
    mine = tmp_path / 'mine' / 'GJ_1'
    _star(str(mine))
    ident = dict(main='GJ 1', aliases=['HIP 1'])
    # not kept there unless you choose to
    assert rvstore.push(ident, str(mine)) == 0
    rvstore.save_settings(True)
    assert rvstore.push(ident, str(mine)) == 4
    files = repo.commits[repo.head]
    assert sorted(files) == ['stars/GJ_1/lbl/GJ_1_HARPSN_lbl.rdb',
                             'stars/GJ_1/manifest.json',
                             'stars/GJ_1/rv/dace/raw_GJ1.csv',
                             'stars/GJ_1/target.json']
    # the local path out of the manifest
    assert '/my/disk' not in repo.blobs[files['stars/GJ_1/manifest.json']
                                        ].decode()
    # nothing new: nothing written
    assert rvstore.push(ident, str(mine)) == 0
    # read by another, under another name of the star: what the folder
    #   lacks, nothing replaced
    theirs = tmp_path / 'theirs' / 'HIP_1'
    (theirs / 'rv' / 'dace').mkdir(parents=True)
    (theirs / 'rv' / 'dace' / 'raw_GJ1.csv').write_text('mine')
    assert rvstore.pull(ident, str(theirs)) == 3
    assert (theirs / 'rv' / 'dace' / 'raw_GJ1.csv').read_text() == 'mine'
    assert (theirs / 'lbl' / 'GJ_1_HARPSN_lbl.rdb').read_text().endswith(
        '1\t2\n')
    # a commit another one got ahead of: asked again from it
    (mine / 'rv' / 'dace' / 'raw_GJ1.csv').write_text('rjd,rv\n1,3\n')
    repo.late = 1
    assert rvstore.push(ident, str(mine)) == 1
    assert 'other' in repo.commits[repo.head]
    # a star never to share: never written
    shared.save_settings(False, ['HIP 1'])
    (mine / 'rv' / 'dace' / 'raw_GJ1.csv').write_text('rjd,rv\n1,4\n')
    assert rvstore.push(ident, str(mine)) == 0
    # every star of a folder
    shared.save_settings(False, [])
    _star(str(tmp_path / 'mine' / 'GJ_2'), 'GJ 2')
    assert dict(rvstore.push_all(str(tmp_path / 'mine'))) == dict(
        GJ_1=1, GJ_2=4)


def test_no_access(tmp_path, monkeypatch):
    """without access (no gh, not a collaborator): as if there were none"""
    monkeypatch.setattr(rvstore, '_gh', lambda *a, **k: (-1, None))
    monkeypatch.setattr(rvstore, '_ACCESS', {})
    rvstore.save_settings(True)
    _star(str(tmp_path / 'GJ_1'))
    assert rvstore.access() == dict(read=False, write=False)
    assert rvstore.pull(dict(main='GJ 1'), str(tmp_path / 'x')) == 0
    assert rvstore.push(dict(main='GJ 1'), str(tmp_path / 'GJ_1')) == 0
