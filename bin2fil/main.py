import os
import struct

import numpy as np
from sampling_qol import ObsParameter, Source

# sanitizar input si es que se conserva esto (improbable si es que pasamos a C)
# añadir inferencia de tipo más adelante?
print("Elija una SDR:\n[1] RTL (8-bit unsigned int)\n[2] Airspy (16-bit signed int)")
sdr_type = int(input("Ponga 1 o 2:"))
# filesize = os.path.getsize(bin_file)
# obs_len_secs = filesize / (2 * channels * sample_rate)
# time_mjd = (os.path.getmtime(bin_file) / 864000) + 40587 - obs_len_secs / 864000


# definición de funciones, hecho por modularidad y porque nos hará la vida un poco más facil si es que pasamos a C
# en particular, definí esta porque hay que hacer el proceso muchas veces.
def bin2cpow(data, off=0, d_type=1, channels=32):
    """
    Función para convertir el archivo de formato binario a poder de canal. Usando la Transformada de Fourier
    Rápida (FFT).

    Args:
        data (str o pathObject): Nombre o path del archivo binario conteniendo la señal.
        off (int): Desfase o posición en el archivo, necesario para los bucles.
        d_type (int): 1 o 2. Tipo de dato contenido en el archivo binario, uint8 por defecto.
        channels (int): Cantidad de canales, idealmente una potencia de 2.
    Returns:
        channel_pow (array-like): Array conteniendo el poder por canal de un solo sample.
    """
    data_type = {1: np.uint8, 2: np.int16}  # tipo de dato dependiendo de la sdr
    samples = np.fromfile(
        data, dtype=data_type[d_type], sep="", count=channels * 2, offset=off
    )  # lectura de los datos, el separador es "empty" solo por precaución para que asuma binario, no sé si es
    # realmente necesario, quizá hacer offset un multiplicador a channel en vez de un número independiente
    freq_data = np.fft.fft(samples.astype(np.float64))  # transformada de Fourier
    channel_pow = np.abs(freq_data) ** 2  # conversión a poder de canal
    return channel_pow


def freq_transformation(file, chunksize):
    samples = filesize / (
        2 * channels
    )  # la cantidad de muestras contenidas en el archivo
    chunks = np.floor(
        samples / chunksize
    )  # la cantidad de trozos en los que dividiremos el proceso
    overflow_chunk_size = samples - chunks * chunksize  # el trozo sobrante
    for chunk in range(
        chunks + 1
    ):  # iteramos en todos los trozos (el +1 es para incluir el sobrante)
        freq_data = np.empty([20, channels])
        offset = 1  # tengo que cambiar esto
        for power in range(20):
            freq_data[power, :] = bin2cpow(file, off=offset * power, channels=channels)
        freq = freq_data.sum(axis=0)  # suma los 20 power samples
        # falta continuar el bucle
        # en este momento no se me occure bien sobre qué tengo que iterar y cuantas veces, pero no es complicado
    return freq


def write_header(obsparams: ObsParameter):
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

        fil.write(struct.pack("<I", 5))
        fil.write(bytearray("nbits", "ascii"))
        fil.write(
            struct.pack("<I", int(8 * sdr_type))
        )  # número de bits por muestra 8 para la rtl y 16 para la airspy

        fil.write(struct.pack("<I", 4))
        fil.write(bytearray("foff", "ascii"))
        fil.write(struct.pack("<d", obsparams.channel_width))  # ancho de cada canal

        fil.write(struct.pack("<I", 4))
        fil.write(bytearray("fch1", "ascii"))  # frecuencia central del primer canal
        fil.write(struct.pack("<d", obsparams.fch1))

        fil.write(struct.pack("<I", 6))
        fil.write(bytearray("nchans", "ascii"))
        fil.write(struct.pack("<I", obsparams.channels))  # cantidad de canales de freq

        fil.write(struct.pack("<I", 5))
        fil.write(bytearray("tsamp", "ascii"))
        fil.write(struct.pack("<d", obsparams.tsample))  # tiempo entre muestras

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
    .

    Args:
        data (str o pathObject): Nombre o path del archivo binario conteniendo la señal.
        off (int): Desfase o posición en el archivo, necesario para los bucles.
        d_type (int): 1 o 2. Tipo de dato contenido en el archivo binario, uint8 por defecto.
        channels (int): Cantidad de canales, idealmente una potencia de 2.
    Returns:
        sum_powers (array-like): Array con la suma de 20 espectros.
        new_off (int): Nuevo desfase o posición en el archivo, necesario para los bucles.
    """

    bytes_per_cycle = (
        channels * 2 * d_type
    )  ###hay que tener ojo aqui por si en algun futuro se cambia d_type por otra cosa que no sea 1, 2
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
    .

    Args:
        data (str o pathObject): Nombre o path del archivo binario conteniendo la señal.
        outfile (str o pathObject): Nombre o path del archivo filterbank de salida.
        d_type (int): 1 o 2. Tipo de dato contenido en el archivo binario, uint8 por defecto.
        channels (int): Cantidad de canales, idealmente una potencia de 2.
    """
    bytes_per_cycle = channels * 2 * d_type
    bytes_per_piece = (
        bytes_per_cycle * 20
    )  # numero de bytes por cada vez que se usa paso3()

    file_length = os.path.getsize(data)
    total_pieces = (
        file_length // bytes_per_piece
    )  # numero de trozos(?), para el ultimo ciclo que recorre todo el archivo

    offset = 0

    with open(outfile, "ab") as fil:
        for piece in range(total_pieces):
            sum_powers, new_off = sum_cpowers(
                data=data, off=offset, d_type=d_type, channels=channels
            )  # solo se usa la función paso3()
            offset = new_off
            sum_powers.astype(np.float32).tofile(fil)


# datos para prueba
vela_pulsar = Source(
    "J0835–4510", 083520.6, -451034.8
)  # creamos el objeto porque solo andamos mirando vela

# pequeña sección de tiempo porque estoy 60% seguro que la fecha de modificación de
# los archivos de drive será cuando los descargué, así que lo estoy haciendo a mano con el historial de drive
# y un poco de inferencia
from astropy.time import Time

mjd1 = Time("2024-05-04 12:30").mjd[0]
tst_fold04 = ObsParameter(2.048e6, mjd1, 400, vela_pulsar, "tstFold04.bin", 1)
tst_fold04.set_channels(32)
tst_fold04.header_data()
write_header(tst_fold04)
file_runthrough(tst_fold04.file, "tst_fold04.fil")
