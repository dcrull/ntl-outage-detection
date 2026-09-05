import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from scripts.ingest_ting import fetch_batch, ingest, iso, profile


class IngestionTests(unittest.TestCase):
    def test_date_conversion(self):
        self.assertEqual(iso(0), '1970-01-01T00:00:00+00:00')
        self.assertIsNone(iso(None))

    def test_transfer_limit_splits_batch(self):
        with patch('scripts.ingest_ting.request', side_effect=[
            {'exceededTransferLimit': True},
            {'features': [{'attributes': {'OBJECTID': 1}}]},
            {'features': [{'attributes': {'OBJECTID': 2}}]},
        ]):
            self.assertEqual(len(fetch_batch('url', [1, 2], 'OBJECTID')), 2)

    def test_missing_and_duplicate_ids_fail(self):
        for values in ([1], [1, 1]):
            with self.subTest(values=values), patch('scripts.ingest_ting.request', return_value={
                'features': [{'attributes': {'OBJECTID': i}} for i in values]
            }):
                with self.assertRaises(RuntimeError):
                    fetch_batch('url', [1, 2], 'OBJECTID')

    def test_event_parts_are_not_deduplicated(self):
        features = [{'attributes': {'eventId': 'a'}}] * 2
        result = profile(features, [])
        self.assertEqual(result['records'], 2)
        self.assertEqual(result['distinct_nonnull_event_ids'], 1)

    def test_snapshot_outputs_and_refuses_overwrite(self):
        metadata = {'objectIdField': 'OBJECTID', 'geometryType': 'esriGeometryPolygon',
                    'fields': [{'name': 'OBJECTID', 'type': 'esriFieldTypeOID'}]}
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / 'snapshot'
            with patch('scripts.ingest_ting.request', side_effect=[
                {}, metadata, {'objectIds': [1]}, {'count': 1},
                {'features': [{'attributes': {'OBJECTID': 1}, 'geometry': {'rings': []}}]},
                {'objectIds': [1]},
            ]):
                ingest(output)
            manifest = json.loads((output / 'manifest.json').read_text())
            self.assertEqual(manifest['status'], 'complete')
            self.assertIn('attributes.csv', manifest['sha256'])
            with self.assertRaises(FileExistsError):
                ingest(output)


if __name__ == '__main__':
    unittest.main()
