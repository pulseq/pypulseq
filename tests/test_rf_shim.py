import numpy as np
import pypulseq as pp
import pytest


# Basic sequence with rf shim
@pytest.fixture
def create_test_data():
    seq = pp.Sequence()
    rf = pp.make_block_pulse(flip_angle=0.5 * np.pi, duration=1e-3)

    # Generate rf shim pattern for an 8ch transmit coil
    n_tx_channels = 8
    phases = np.arange(n_tx_channels, dtype=float) / n_tx_channels * 2 * np.pi
    amplitudes = np.ones(n_tx_channels, dtype=float)

    # Shims
    shim_vec = amplitudes * np.exp(1j * phases)

    seq.add_block(rf, pp.make_rf_shim(shim_vec))

    return seq, shim_vec


def test_roundtrip_shim(create_test_data, tmp_path):
    seq, expected = create_test_data

    # Create a temporary file path using tmp_path
    filename = tmp_path / 'test_sequence.seq'

    # Write the sequence to file
    seq.write(str(filename))

    # Read the sequence back
    seq_loaded = pp.Sequence(system=seq.system)
    seq_loaded.read(filename)

    # Extract shim vectors from the loaded sequence
    loaded_shims = []
    for block_counter in seq_loaded.block_events:
        block = seq_loaded.get_block(block_counter)
        if hasattr(block, 'rf_shim'):
            loaded_shims.append(block.rf_shim.shim_vector)

    # There should be one shim in this test
    assert len(loaded_shims) == 1

    # Compare loaded shim with original shim (magnitudes and phases are written with 6 significant digits)
    np.testing.assert_allclose(loaded_shims[0], expected, atol=1e-5)
