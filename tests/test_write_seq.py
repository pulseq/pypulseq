"""Tests for the write_seq module"""

import hashlib

import pypulseq as pp
import pytest


@pytest.mark.parametrize('v141_compat', [False, True])
def test_write_seq_line_endings_and_signature(tmp_path, v141_compat):
    seq = pp.Sequence()
    seq.add_block(pp.make_trapezoid('x', area=1000))
    seq.add_block(pp.make_delay(1e-3))
    seq.add_block(pp.make_trapezoid('y', area=-500), pp.make_adc(num_samples=50, duration=1e-3))

    file_name = tmp_path / 'test.seq'
    md5 = seq.write(file_name, v141_compat=v141_compat)

    content = file_name.read_bytes()
    # .seq files must be byte-identical across platforms: LF only, never CRLF
    assert b'\r' not in content

    # The signature hash must match the md5 of the file content up to (and excluding) the
    # newline preceding [SIGNATURE], as documented in the file itself.
    body, signature = content.split(b'\n[SIGNATURE]\n')
    assert hashlib.md5(body).hexdigest() == md5
    assert f'Hash {md5}'.encode() in signature
