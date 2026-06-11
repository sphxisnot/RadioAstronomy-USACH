import os
import struct

import numpy as np

# por ahora haré que el filename sea escrito por el usuario, considerar implementar un walk quizá
bin_file = input("Escribir nombre del archivo .bin a convertir:\n")
# sanitizar input si es que se conserva esto (improbable si es que pasamos a C)
# añadir inferencia de tipo más adelante?
print("Elija una SDR:\n[1] RTL (8-bit unsigned int)\n[2] Airspy (16-bit signed int)")
sdr_type = int(input("Ponga 1 o 2:"))


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
    channel_pow = freq_data.real**2 + freq_data.imag**2  # conversión a poder de canal
    return channel_pow


# hay que hacer un bucle que itere sobre los datos, no sé si un while con un contador o un for sobre el tamaño
# del archivo. De todas maneras sería usando información del archivo usando os, no leyendo el elemento entero
# (por el bien de la memoria y mi consciencia)
def write_header(file):
    outfile = str(file).rstrip(".bin") + ".fil"
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
        fil.write(struct.pack("<I", 0))  # 0 es el valor para datos artificiales

        fil.write(struct.pack("<I", 5))
        fil.write(bytearray("nbits", "ascii"))
        fil.write(
            struct.pack("<I", nbits)
        )  # número de bits por muestra 8 para la rtl y 16 para la airspy

        fil.write(struct.pack("<I", 4))
        fil.write(bytearray("foff", "ascii"))
        fil.write(struct.pack("<d", channel_w))  # ancho de cada canal

        fil.write(struct.pack("<I", 4))
        fil.write(bytearray("fch1", "ascii"))  # frecuencia central del primer canal
        fil.write(struct.pack("<d", freq_ch1))

        # fil.write()
