import os
import struct
import numpy as np
from sampling_qol import ObsParameter
# ============================================================
# FIX (ver diagnóstico): bin2cpow tenía dos bugs que corrompían
# el filterbank de salida:
#1. Solo leía `channels` valores del archivo en vez de
#  `channels*2`, y los pasaba a la FFT como si fueran una
#   señal REAL.
# 2. La salida de np.fft.fft() viene en orden "DC, +1, +2, ...,
#    Nyquist-1, -Nyquist, ..., -1", no en orden ascendente de
#    frecuencia.

def bin2cpow(data, off=0, d_type=1, channels=32, dc_offset=None):
    """
    Args:
        data: ruta al archivo I/Q, O BIEN un np.memmap/ndarray ya abierto
            sobre el (ver nota de rendimiento abajo).
    """
    data_type = {1: np.uint8, 2: np.int16}  # tipo de dato dependiendo de la sdr
    np_dtype = data_type[d_type]

    # Cada muestra compleja (I, Q) ocupa 2 valores reales en el archivo, por eso hay
    # que leer channels*2 valores para formar `channels` muestras complejas
    n_iq_values = channels * 2

    if isinstance(data, (str, os.PathLike)):
        raw = np.fromfile(
            data, dtype=np_dtype, sep="", count=n_iq_values, offset=off
        )
    else:
        start = off // np.dtype(np_dtype).itemsize
        raw = np.asarray(data[start:start + n_iq_values])

    if raw.size < n_iq_values:
        raise EOFError(
            f"No quedan suficientes muestras en {data!r} para completar un chunk "
            f"({raw.size} leidas, se esperaban {n_iq_values})."
        )

    raw = raw.astype(np.float64)
    if dc_offset is None:
        dc_offset = 127.5 if d_type == 1 else 0.0
    raw -= dc_offset

    # Reconstruccion I/Q: el archivo viene intercalado I,Q,I,Q,...
    samples = raw[0::2] + 1j * raw[1::2]

    freq_data = np.fft.fft(samples)
    freq_data = np.fft.fftshift(freq_data)
    channel_pow = np.abs(freq_data) ** 2
    return channel_pow


def write_header(obsparams: ObsParameter, is_presto=False, group_chunks=4):
    file = (
        str(obsparams.file).removesuffix(".iq").removesuffix(".bin")
    )
    outfile = file + ".fil"
    with open(outfile, "wb") as fil:
        fil.write(struct.pack("<I", 12))
        fil.write(bytearray("HEADER_START", "ascii"))

        fil.write(struct.pack("<I", 9))
        fil.write(bytearray("data_type", "ascii"))
        fil.write(struct.pack("<I", 1))

        fil.write(struct.pack("<I", 4))
        fil.write(bytearray("nifs", "ascii"))
        fil.write(struct.pack("<I", 1))

        fil.write(struct.pack("<I", 12))
        fil.write(bytearray("telescope_id", "ascii"))
        fil.write(struct.pack("<I", 0))


        fil.write(struct.pack("<I", 5))
        fil.write(bytearray("nbits", "ascii"))
        fil.write(struct.pack("<I", 32))


        bandwidth_MHz = obsparams.sample_rate / 1e6
        channel_width_MHz = bandwidth_MHz / obsparams.channels
        # centro del canal más ALTO de la banda:
        fch1_MHz = obsparams.cfreq + bandwidth_MHz / 2 - channel_width_MHz / 2

        fil.write(struct.pack("<I", 4))
        fil.write(bytearray("foff", "ascii"))
        fil.write(struct.pack("<d", -channel_width_MHz))

        fil.write(struct.pack("<I", 4))
        fil.write(bytearray("fch1", "ascii"))
        fil.write(struct.pack("<d", fch1_MHz + int(is_presto) * 992e3))

        fil.write(struct.pack("<I", 6))
        fil.write(bytearray("nchans", "ascii"))
        fil.write(struct.pack("<I", obsparams.channels))

        tsamp_row = obsparams.channels * group_chunks / obsparams.sample_rate

        fil.write(struct.pack("<I", 5))
        fil.write(bytearray("tsamp", "ascii"))
        fil.write(struct.pack("<d", tsamp_row))

        fil.write(struct.pack("<I", 6))
        fil.write(bytearray("tstart", "ascii"))
        fil.write(
            struct.pack("<d", obsparams.obstime)
        )

        fil.write(struct.pack("<I", 11))
        fil.write(bytearray("source_name", "ascii"))
        fil.write(struct.pack("<I", len(obsparams.source)))
        fil.write(bytearray(obsparams.source, "ascii"))

        fil.write(struct.pack("<I", 7))
        fil.write(bytearray("src_raj", "ascii"))
        fil.write(struct.pack("<d", obsparams.ra))

        fil.write(struct.pack("<I", 7))
        fil.write(bytearray("src_dej", "ascii"))
        fil.write(struct.pack("<d", obsparams.dec))

        fil.write(struct.pack("<I", 10))
        fil.write(bytearray("HEADER_END", "ascii"))


def sum_cpowers(
    data, off, d_type=1, channels=32, group_chunks=4
):
    bytes_per_cycle = (
        channels * 2 * d_type
    )
    new_off = off

    sum_powers = np.zeros(
        channels, dtype=np.float64
    )

    for _ in range(group_chunks):
        power = bin2cpow(data=data, off=new_off, d_type=d_type, channels=channels)
        sum_powers += power
        new_off += bytes_per_cycle
    return sum_powers, new_off


def file_runthrough(data, outfile, d_type=1, channels=32, group_chunks=4):
    #Recorre el archivo binario completo, chunk por chunk, escribiendo el filterbank.

    bytes_per_cycle = channels * 2 * d_type
    bytes_per_piece = (
        bytes_per_cycle * group_chunks
    )

    file_length = os.path.getsize(data)
    total_pieces = (
        file_length // bytes_per_piece
    )  # numero de trozos, para el ultimo ciclo que recorre todo el archivo

    # Se abre el archivo UNA sola vez y se le pasa el memmap a bin2cpow, en vez
    # de dejar que cada chunk lo reabra por su cuenta
    raw_map = np.memmap(data, dtype={1: np.uint8, 2: np.int16}[d_type], mode="r")

    offset = 0

    with open(outfile, "ab") as fil:
        for piece in range(total_pieces):
            sum_powers, new_off = sum_cpowers(
                data=raw_map, off=offset, d_type=d_type, channels=channels,
                group_chunks=group_chunks,
            )
            offset = new_off
            sum_powers[::-1].astype(np.float32).tofile(fil)