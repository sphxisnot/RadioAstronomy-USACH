import os
import struct
import numpy as np
from sampling_qol import ObsParameter

# ============================================================
# FIX (ver diagnóstico): bin2cpow tenía dos bugs que corrompían
# el filterbank de salida:
#
#   1. Solo leía `channels` valores del archivo en vez de
#      `channels*2`, y los pasaba a la FFT como si fueran una
#      señal REAL. El offset entre llamadas SÍ avanzaba en
#      `channels*2*d_type` (ver bytes_per_cycle en sum_cpowers),
#      así que en la práctica se descartaba la mitad de cada
#      chunk sin nunca leerlo, y lo que sí se leía eran I y Q
#      intercalados tratados como una serie real -> espectro
#      corrupto (tono/leakage espurio, sobre todo si el SDR es
#      uint8 con DC offset ~127.5, ver punto 2).
#
#   2. La salida de np.fft.fft() viene en orden "DC, +1, +2, ...,
#      Nyquist-1, -Nyquist, ..., -1", no en orden ascendente de
#      frecuencia. El header (fch1/foff/nchans) escrito por
#      write_header asume implícitamente orden monótono de canal
#      -> frecuencia. Sin fftshift, PRESTO/sigpyproc aplican el
#      delay de dedispersión al canal equivocado en cada
#      frecuencia, lo que destruye selectivamente la señal al
#      plegar (no baja el SNR parejo, lo desalinea).


def bin2cpow(data, off=0, d_type=1, channels=32, dc_offset=None):
    """
    Returns:
        channel_pow (array-like): Array conteniendo el poder por canal de un solo sample,
            en orden ASCENDENTE de frecuencia (canal 0 = frecuencia más baja del bloque),
            consistente con la convención fch1 + i*foff usada en write_header.
    """
    data_type = {1: np.uint8, 2: np.int16}  # tipo de dato dependiendo de la sdr

    # Cada muestra compleja (I, Q) ocupa 2 valores reales en el archivo, por eso hay
    # que leer channels*2 valores para formar `channels` muestras complejas -- esto
    # tiene que calzar con bytes_per_cycle en sum_cpowers (channels*2*d_type).
    n_iq_values = channels * 2
    raw = np.fromfile(
        data, dtype=data_type[d_type], sep="", count=n_iq_values, offset=off
    )
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
    samples = raw[0::2] + 1j * raw[1::2]  # ahora sí son `channels` muestras complejas

    freq_data = np.fft.fft(samples)  # transformada de Fourier
    freq_data = np.fft.fftshift(freq_data)  # reordena a frecuencia ascendente (-Nyq..+Nyq)
    channel_pow = np.abs(freq_data) ** 2  # conversión a poder de canal
    return channel_pow


def write_header(obsparams: ObsParameter, is_presto=False, group_chunks=20):
    file = (
        str(obsparams.file).removesuffix(".iq").removesuffix(".bin")
    )  # me di cuenta que tengo que chequear si es bin o iq
    outfile = file + ".fil"
    with open(outfile, "wb") as fil:
        fil.write(struct.pack("<I", 12))
        fil.write(bytearray("HEADER_START", "ascii"))

        fil.write(struct.pack("<I", 9))
        fil.write(bytearray("data_type", "ascii"))
        fil.write(struct.pack("<I", 1))  # tipo de dato, 1 es filterbank

        fil.write(struct.pack("<I", 4))
        fil.write(bytearray("nifs", "ascii"))
        fil.write(struct.pack("<I", 1))

        fil.write(struct.pack("<I", 12))
        fil.write(bytearray("telescope_id", "ascii"))
        fil.write(struct.pack("<I", 0))  # 0 es el valor para datos artificiales

        # ============================================================
        # FIX CONFIRMADO CON PRESTO (readfile sobre tstFold04.fil):
        # nbits se tomaba del SDR crudo (8 para RTL-SDR, 16 para Airspy),
        # pero eso describe el formato de ENTRADA, no el de salida. El
        # conversor escribe poder de canal como float32 (4 bytes), asi que
        # nbits DEBE ser 32.
        #
        # Con nbits=8, PRESTO leia cada byte individual como una muestra:
        # reportaba 384000 espectros en vez de 96000 (4x de mas), duracion
        # de 120 s en vez de 30 s, y los VALORES leidos eran basura
        # (fragmentos de float32 interpretados como enteros de 8 bits).
        fil.write(struct.pack("<I", 5))
        fil.write(bytearray("nbits", "ascii"))
        fil.write(struct.pack("<I", 32))  # float32: el conversor escribe .astype(np.float32)

        # ============================================================
        # FIX 1 (confirmado con tstFold04.params): header_data() en
        # sampling_qol.py calculaba channel_width en Hz
        # (sample_rate/channels) y lo restaba directamente de cfreq, que
        # viene en MHz del .params -> fch1 salía -31595.904 (sin sentido
        # físico). Aquí se recalcula todo consistentemente en MHz.
        #
        # FIX 2 (confirmado con readfile de PRESTO): la convención sigproc
        # es que fch1 es el canal MÁS ALTO y foff es NEGATIVO (la banda
        # desciende). Al escribir foff positivo esperando orden ascendente,
        # PRESTO igual interpretó fch1 como el canal más alto y leyó la
        # banda como 397.024-399.008 MHz (centrada en 398.016) en vez de
        # 398.976-401.024 (centrada en 400) -> desplazada 1.984 MHz.
        #
        # Por eso ahora se escribe en convención sigproc estándar, y el
        # conversor debe entregar los canales en orden DESCENDENTE de
        # frecuencia (ver flip_to_descending en converter_optimized.py).
        bandwidth_MHz = obsparams.sample_rate / 1e6
        channel_width_MHz = bandwidth_MHz / obsparams.channels
        # centro del canal más ALTO de la banda:
        fch1_MHz = obsparams.cfreq + bandwidth_MHz / 2 - channel_width_MHz / 2

        fil.write(struct.pack("<I", 4))
        fil.write(bytearray("foff", "ascii"))
        fil.write(struct.pack("<d", -channel_width_MHz))  # NEGATIVO: banda descendente

        fil.write(struct.pack("<I", 4))
        fil.write(bytearray("fch1", "ascii"))  # canal más alto, en MHz
        fil.write(struct.pack("<d", fch1_MHz + int(is_presto) * 992e3))

        fil.write(struct.pack("<I", 6))
        fil.write(bytearray("nchans", "ascii"))
        fil.write(struct.pack("<I", obsparams.channels))  # cantidad de canales de freq

        # ============================================================
        # FIX CONFIRMADO CON DATOS REALES: obsparams.tsample = 1/sample_rate
        # es el tiempo entre muestras I/Q CRUDAS, pero cada fila que
        # realmente se escribe en el .fil integra `group_chunks` espectros
        # de `channels` muestras cada uno. Con el tsamp crudo, un programa
        # que lea este .fil (PRESTO, sigpyproc) cree que la observación
        # dura ~640x menos de lo real (confirmado: 0.047s en vez de 30s
        # reales en tstFold04.fil). tsamp real de cada fila:
        tsamp_row = obsparams.channels * group_chunks / obsparams.sample_rate

        fil.write(struct.pack("<I", 5))
        fil.write(bytearray("tsamp", "ascii"))
        fil.write(struct.pack("<d", tsamp_row))  # tiempo real entre filas escritas

        fil.write(struct.pack("<I", 6))
        fil.write(bytearray("tstart", "ascii"))
        fil.write(
            struct.pack("<d", obsparams.obstime)
        )  # tiempo de inicio de la medición

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
    data, off, d_type=1, channels=32
):  # hay que cambiar el nombre de las funciones por algo creativo
    """
    Suma 20 espectros consecutivos (ver nota HawkRAO abajo) para un mismo chunk.
    """

    bytes_per_cycle = (
        channels * 2 * d_type
    )  # channels muestras complejas = channels*2 valores reales en el archivo
    new_off = off

    sum_powers = np.zeros(
        channels, dtype=np.float64
    )  # mal nombre, variable para acumular las potencias en el ciclo

    for _ in range(
        20
    ):  # numero arbitrario que aparece en las notas tecnicas de HawkRAO
        power = bin2cpow(data=data, off=new_off, d_type=d_type, channels=channels)
        sum_powers += power
        new_off += bytes_per_cycle
    return sum_powers, new_off


def file_runthrough(data, outfile, d_type=1, channels=32):
    """
    Recorre el archivo binario completo, chunk por chunk, escribiendo el filterbank.
    """
    bytes_per_cycle = channels * 2 * d_type
    bytes_per_piece = (
        bytes_per_cycle * 20
    )  # numero de bytes por cada vez que se usa sum_cpowers()

    file_length = os.path.getsize(data)
    total_pieces = (
        file_length // bytes_per_piece
    )  # numero de trozos, para el ultimo ciclo que recorre todo el archivo

    offset = 0

    with open(outfile, "ab") as fil:
        for piece in range(total_pieces):
            sum_powers, new_off = sum_cpowers(
                data=data, off=offset, d_type=d_type, channels=channels
            )
            offset = new_off
            sum_powers.astype(np.float32).tofile(fil)