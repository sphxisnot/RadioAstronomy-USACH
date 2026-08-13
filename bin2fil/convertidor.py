"""
conversion I/Q -> poder de canal que main_fixed.py (reconstruccion I/Q,
resta de DC offset, fftshift para orden ascendente de frecuencia), pero
vectorizada y procesada por batches.
La ganancia viene de dos cosas:
    1. Leer el archivo UNA vez (np.memmap) en vez de abrirlo/buscarlo
       ~1.9 millones de veces.
    2. Hacer la FFT de un batch de miles de chunks de una sola llamada
       a np.fft.fft(..., axis=-1) en vez de un chunk a la vez -- numpy
       vectoriza esto internamente (BLAS/FFTW por debajo), mucho mas
       rapido que el mismo trabajo repetido en un loop de Python.

Esto ya deja el conversor rapido en Python puro; la Tarea 7
(reescritura en C) seguiria siendo relevante solo si esto no alcanza
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

        # FIX confirmado con readfile de PRESTO: la convención sigproc es
        # banda DESCENDENTE (fch1 = canal más alto, foff negativo). Tras el
        # fftshift los canales quedan en orden ascendente de frecuencia, así
        # que hay que invertirlos para calzar con el header que escribe
        # write_header() en main_fixed.py. Sin este flip, cada canal queda
        # asociado a la frecuencia equivocada y la dedispersión corrige mal.
        integrated = integrated[:, ::-1]

        power[idx:idx + g] = integrated.astype(np.float32)

        idx += g
        complex_offset += n_complex_batch

    return power, tsamp_integrated


def write_filterbank_body(power, outfile_path, mode="ab"):
    """
    Escribe el cuerpo del filterbank (los datos, sin header) a partir del
    array (n_groups, channels) devuelto por convert_iq_to_filterbank_power.
    Usar DESPUES de escribir el header con write_header() de main_fixed.py
    sobre el mismo archivo, en modo append -- mismo flujo que
    file_runthrough() original.
    """
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
    en main_fixed.py
    """
    power, tsamp_integrated = convert_iq_to_filterbank_power(
        input_path, sample_rate, channels=channels, group_chunks=group_chunks,
        dtype_bits=dtype_bits, batch_groups=batch_groups,
    )
    write_filterbank_body(power, outfile_path)
    return power.shape[0], tsamp_integrated