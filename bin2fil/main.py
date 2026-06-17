import os
import struct

import numpy as np

# por ahora haré que el filename sea escrito por el usuario, considerar implementar un walk quizá
bin_file = input("Escribir nombre del archivo .bin a convertir:\n")
# sanitizar input si es que se conserva esto (improbable si es que pasamos a C)
# añadir inferencia de tipo más adelante?
print("Elija una SDR:\n[1] RTL (8-bit unsigned int)\n[2] Airspy (16-bit signed int)")
sdr_type = int(input("Ponga 1 o 2:"))
sample_rate = 2.048e6
filesize = os.path.getsize(bin_file)
obs_len_secs = filesize / (2 * channels * sample_rate)
time_mjd = (os.path.getmtime(bin_file) / 864000) + 40587 - obs_len_secs / 864000


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
    freq_data = np.fft.fft(samples)  # transformada de Fourier
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


def write_header(file):
    outfile = str(file).rstrip(".bin") + ".fil"
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
        fil.write(struct.pack("<d", channel_w))  # ancho de cada canal

        fil.write(struct.pack("<I", 4))
        fil.write(bytearray("fch1", "ascii"))  # frecuencia central del primer canal
        fil.write(struct.pack("<d", freq_ch1))

        fil.write(struct.pack("<I", 6))
        fil.write(bytearray("nchans", "ascii"))
        fil.write(struct.pack("<I", channels))  # cantidad de canales de freq

        fil.write(struct.pack("<I", 5))
        fil.write(bytearray("tsamp", "ascii"))
        fil.write(struct.pack("<d", tsample))  # tiempo entre muestras

        fil.write(struct.pack("<I", 6))
        fil.write(bytearray("tstart", "ascii"))
        fil.write(struct.pack("<d", time_mjd))  # tiempo de inicio de la medición

        fil.write(struct.pack("<I", 11))
        fil.write(bytearray("source_name", "ascii"))
        fil.write(struct.pack("<I", len(source_name)))
        fil.write(bytearray(source_name, "ascii"))

        fil.write(struct.pack("<I", 7))
        fil.write(bytearray("src_raj", "ascii"))
        fil.write(struct.pack("<d", source_ra))

        fil.write(struct.pack("<I", 7))
        fil.write(bytearray("src_dej", "ascii"))
        fil.write(struct.pack("<d", source_dec))

        fil.write(struct.pack("<I", 10))
        fil.write(bytearray("HEADER_END", "ascii"))
