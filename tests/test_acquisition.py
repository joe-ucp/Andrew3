"""A reacquisition must never silently replace a frozen source commitment."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from src import fetch_primary

class AcquisitionTests(unittest.TestCase):
    def test_retained_manifest_controls_source_and_date(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'evidence').mkdir();(root/'data').mkdir()
            source=root/'data/test.bin';source.write_bytes(b'original')
            old=dict(filename='data/test.bin',url='https://example.org/frozen',
                     sha256=hashlib.sha256(b'original').hexdigest(),retrieved_utc='frozen-acquisition')
            (root/'evidence/primary_sources.json').write_text(json.dumps([old]))
            with patch.object(fetch_primary,'ROOT',root):
                self.assertEqual(fetch_primary.fetch(old['url'],old['filename'],'test'),old)
                source.write_bytes(b'changed')
                with self.assertRaisesRegex(ValueError,'Source hash changed'):
                    fetch_primary.fetch(old['url'],old['filename'],'test')
                self.assertEqual(json.loads((root/'evidence/primary_sources.json').read_text()),[old])

if __name__=='__main__':unittest.main()
