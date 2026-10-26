"""Provider downloads must satisfy the existing commitments before installation."""
import hashlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile

from src.fetch_public_inputs import acquire_one,checked_write
from src import fetch_raw_counts
from src.source_inventory import source_scope,inventory_path,validate_source_locator,MANIFESTS,FIELDS

class PublicAcquisitionTests(unittest.TestCase):
    def test_changed_provider_bytes_are_not_installed(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            record=dict(filename='data/raw/test.bin',url='https://example.org/input',sha256=hashlib.sha256(b'original').hexdigest())
            with patch('src.fetch_public_inputs.read_provider',return_value=b'changed'):
                with self.assertRaisesRegex(ValueError,'Provider/source drift'):
                    acquire_one(root,[record])
            self.assertFalse((root/record['filename']).exists())

    def test_manual_import_checks_bytes_and_shared_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)/'checkout';root.mkdir()
            manual=Path(directory)/'downloads';manual.mkdir()
            (manual/'same.rpt').write_bytes(b'frozen')
            common=dict(url='https://example.org/same.rpt',sha256=hashlib.sha256(b'frozen').hexdigest())
            records=[dict(filename='data/raw/a/same.rpt',**common),dict(filename='data/raw/b/same.rpt',**common)]
            with patch('src.fetch_public_inputs.read_provider',side_effect=AssertionError('Unexpected network')):
                result=acquire_one(root,records,manual)
            self.assertTrue(all(r['status']=='imported_manual_verified' for r in result))
            self.assertTrue(all((root/r['filename']).read_bytes()==b'frozen' for r in records))
            self.assertEqual((manual/'same.rpt').read_bytes(),b'frozen')

    def test_archive_wrapper_does_not_replace_workbook_commitment(self):
        with tempfile.TemporaryDirectory() as directory:
            buf=io.BytesIO()
            with zipfile.ZipFile(buf,'w') as archive:
                archive.writestr('mmc2.xlsx',b'exact workbook bytes')
                archive.writestr('provider-notice.txt',b'Container may be regenerated')
            record=dict(filename='data/raw/scores/mmc2.xlsx',url='https://example.org/archive',
                archive_member='mmc2.xlsx',sha256=hashlib.sha256(b'exact workbook bytes').hexdigest())
            with patch('src.fetch_public_inputs.read_provider',return_value=buf.getvalue()):
                result=acquire_one(Path(directory),[record])
            self.assertEqual(result[0]['download_container_sha256'],hashlib.sha256(buf.getvalue()).hexdigest())
            self.assertEqual((Path(directory)/record['filename']).read_bytes(),b'exact workbook bytes')

    def test_retained_drift_and_unsafe_paths_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'input.bin').write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError,'Existing local input'):
                checked_write(root,'input.bin',b'original',hashlib.sha256(b'original').hexdigest())
            with self.assertRaises(ValueError):
                checked_write(root,'../escape.bin',b'original',hashlib.sha256(b'original').hexdigest())
            self.assertEqual((root/'input.bin').read_bytes(),b'changed')

    def test_legacy_raw_helper_preserves_manifest_and_rejects_remote_drift(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'evidence').mkdir()
            old=dict(filename='data/raw/counts/input.bin',url='https://example.org/input',
                sha256=hashlib.sha256(b'original').hexdigest(),retrieved_utc='frozen-date')
            manifest=root/'evidence/raw_count_sources.json';manifest.write_text(json.dumps([old]),encoding='utf-8')
            before=manifest.read_bytes()
            with patch.object(fetch_raw_counts,'ROOT',root),patch('urllib.request.urlopen',return_value=io.BytesIO(b'changed')):
                with self.assertRaisesRegex(ValueError,'Remote source drift'):
                    fetch_raw_counts.fetch(('input.bin',old['url']))
            self.assertFalse((root/old['filename']).exists())
            self.assertEqual(manifest.read_bytes(),before)

    def test_scope_catalogs_are_explicit_and_separate(self):
        root=Path('checkout')
        self.assertEqual(inventory_path(root,'full'),root/'sources.tsv')
        self.assertEqual(inventory_path(root,'computational'),root/'evidence/computational_sources.tsv')
        with self.assertRaises(ValueError):source_scope('skip-checks')

    def test_historical_locator_requires_a_real_commitment_and_checks_present_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);filename='data/raw/literature/article.xml'
            expected=hashlib.sha256(b'original').hexdigest()
            for manifest in MANIFESTS:
                path=root/manifest;path.parent.mkdir(parents=True,exist_ok=True)
                records=[dict(filename=filename,sha256=expected,bytes=8)] if manifest==MANIFESTS[2] else []
                path.write_text(json.dumps(records),encoding='utf-8')
            row=dict(filename=filename,sha256=expected,bytes=8,source_manifest=MANIFESTS[2],
                url_or_accession='https://example.org/article',retrieval_utc='frozen-date',
                retrieval_basis='recorded',access_note='public',role='primary article')
            with (root/'sources.tsv').open('w',encoding='utf-8',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=FIELDS,delimiter='\t');writer.writeheader();writer.writerow(row)
            with self.assertRaisesRegex(ValueError,'Unknown historical source locator'):
                validate_source_locator(root,'data/raw/literature/invented.xml','computational')
            self.assertEqual(validate_source_locator(root,filename,'computational')['status'],
                'historical_commitment_only_payload_not_fetched')
            with self.assertRaises(FileNotFoundError):validate_source_locator(root,filename,'full')
            (root/filename).write_bytes(b'changed!')
            with self.assertRaisesRegex(ValueError,'payload differs'):
                validate_source_locator(root,filename,'computational')
            (root/filename).write_bytes(b'original')
            self.assertEqual(validate_source_locator(root,filename,'full')['status'],'payload_hash_and_size_verified')

if __name__=='__main__':unittest.main()
