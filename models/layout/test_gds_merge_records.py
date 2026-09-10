"""GDS记录级合并解析的边界测试，不依赖用户源版图。"""
import struct
import unittest
from build_15mm_ysj_merge import gds_records, structures


def record(kind,payload=b'',dtype=0):
    if len(payload)%2:payload+=b'\0'
    return struct.pack('>HBB',len(payload)+4,kind,dtype)+payload


class RecordTests(unittest.TestCase):
    def test_unsigned_record_length(self):
        raw=record(16,b'\0'*40000,3)
        self.assertEqual(list(gds_records(raw))[0][2],raw)

    def test_structure_bytes_preserved(self):
        cell=record(5,b'\0'*24,2)+record(6,b'ABC',6)+record(7)
        library=record(3,b'\0'*16,5)+cell+record(4)
        cells,units,end=structures(library)
        self.assertEqual(cells,{'ABC':cell})
        self.assertEqual(units,b'\0'*16)
        self.assertEqual(end,record(4))

    def test_reject_truncated_and_duplicate(self):
        with self.assertRaises(ValueError):list(gds_records(b'\x00\x08\x00\x00'))
        cell=record(5,b'\0'*24,2)+record(6,b'ABC',6)+record(7)
        with self.assertRaises(ValueError):structures(record(3,b'\0'*16,5)+cell+cell+record(4))


if __name__=='__main__':unittest.main()
