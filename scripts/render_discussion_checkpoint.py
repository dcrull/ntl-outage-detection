#!/usr/bin/env python3
"""Execute the checkpoint notebook and export a standalone HTML discussion copy."""
import json
from pathlib import Path
import sys

import nbformat
from nbclient import NotebookClient
from nbconvert import HTMLExporter

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts.discussion_checkpoint import file_hash


def main():
    path = ROOT / 'notebooks/06_discussion_checkpoint.ipynb'
    notebook = nbformat.read(path, as_version=4)
    starting_hash = file_hash(path)

    def progress(cell, cell_index, **kwargs):
        if cell.cell_type == 'code':
            print(f'Running cell {cell_index}: {cell.source.splitlines()[0]}', flush=True)

    client = NotebookClient(notebook, timeout=1800, resources={'metadata': {'path': str(ROOT)}},
                            on_cell_start=progress)
    try:
        client.execute()
    finally:
        # An interactive edit during rendering belongs to the user. Keep it and
        # write the executed version alongside it instead of overwriting it.
        if file_hash(path) != starting_hash:
            path = path.with_name(path.stem + '.rendered.ipynb')
            print('Working notebook changed during execution; saving rendered copy:', path, flush=True)
        nbformat.write(notebook, path)
    nbformat.validate(notebook)
    out = ROOT / 'data/published/discussion_checkpoint'
    manifest_path = out / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    if not all(row['converged'] for row in manifest['cluster_summaries']):
        manifest['status'] = 'needs_cluster_review'
        manifest_path.write_text(json.dumps(manifest, indent=2))
        raise RuntimeError('At least one clustering run reached its iteration limit; inspect before sharing')
    html, _ = HTMLExporter(exclude_input=True).from_notebook_node(notebook)
    html_path = out / 'discussion_checkpoint.html'
    html_path.write_text(html)
    manifest['output_sha256'][html_path.name] = file_hash(html_path)
    manifest['notebook_sha256'] = file_hash(path)
    manifest['renderer_sha256'] = file_hash(Path(__file__))
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print('Completed notebook and standalone HTML:', html_path, flush=True)


if __name__ == '__main__':
    main()
