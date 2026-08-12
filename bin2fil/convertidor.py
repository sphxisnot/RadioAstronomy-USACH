"""
Reemplazo de sum_cpowers()/file_runthrough() (main_fixed.py) para la
Tarea 6 del board (Optimizacion del conversor). Misma logica de
conversion I/Q -> poder de canal que main_fixed.py (reconstruccion I/Q,
resta de DC offset, fftshift para orden ascendente de frecuencia), pero
vectorizada y procesada por batches en vez de:

    - abrir el archivo con np.fromfile(path, ...) una vez POR CADA chunk
      de `channels` muestras (bin2cpow original), y
    - hacer un solo FFT de `channels` puntos a la vez dentro de un loop
      de Python.

Esto ya deja el conversor rapido en Python puro; la Tarea 7
(reescritura en C) seguiria siendo relevante solo si esto no alcanza
para observaciones muy largas o para uso en tiempo real, pero para
post-procesamiento offline (que es el caso de uso actual, PRESTO
corre despues sobre el .fil ya escrito) esto deberia ser suficiente.
"""

import numpy as np


def convert_iq_to_filterbank_power(
    input_path,
    sample_rate,
    channels=32,
    group_chunks=20,
    dtype_bits=8,
    batch_groups=5000,
):
    """
    Convierte un archivo I/Q crudo completo en un array 2D de poder por
    canal (tiempo x canal), vectorizado y por batches.
    """
    np_dtype = np.uint8 if dtype_bits == 8 else np.int16
    dc_offset = 127.5 if dtype_bits == 8 else 0.0

    raw = np.memmap(input_path, dtype=np_dtype, mode="r")
    n_complex_total = raw.shape[0] // 2

    samples_per_group = channels * group_chunks
    n_groups_total = n_complex_total // samples_per_group
    tsamp_integrated = samples_per_group / sample_rate

    power = np.empty((n_groups_total, channels), dtype=np.float32)

    idx = 0
    complex_offset = 0
    while idx < n_groups_total:
        g = min(batch_groups, n_groups_total - idx)
        n_complex_batch = g * samples_per_group

        raw_batch = np.asarray(
            raw[complex_offset * 2: (complex_offset + n_complex_batch) * 2],
            dtype=np.float64,
        )
        raw_batch -= dc_offset
        iq = raw_batch[0::2] + 1j * raw_batch[1::2]

        blocks = iq.reshape(g, group_chunks, channels)
        fft_blocks = np.fft.fftshift(np.fft.fft(blocks, axis=2), axes=2)
        pow_blocks = np.abs(fft_blocks) ** 2          # (g, group_chunks, channels)
        integrated = pow_blocks.sum(axis=1)            # (g, channels)

        power[idx:idx + g] = integrated.astype(np.float32)

        idx += g
        complex_offset += n_complex_batch

    return power, tsamp_integrated


def write_filterbank_body(power, outfile_path, mode="ab"):

    with open(outfile_path, mode) as fil:
        power.astype(np.float32).tofile(fil)


def convert_and_write(
    input_path,
    outfile_path,
    sample_rate,
    channels=32,
    group_chunks=20,
    dtype_bits=8,
    batch_groups=5000,
):
    """
    Atajo: hace convert_iq_to_filterbank_power + write_filterbank_body en
    un solo llamado. Pensado como reemplazo directo de file_runthrough()
    """
    power, tsamp_integrated = convert_iq_to_filterbank_power(
        input_path, sample_rate, channels=channels, group_chunks=group_chunks,
        dtype_bits=dtype_bits, batch_groups=batch_groups,
    )
    write_filterbank_body(power, outfile_path)
    return power.shape[0], tsamp_integrated