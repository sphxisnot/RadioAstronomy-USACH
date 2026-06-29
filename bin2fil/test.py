from main import file_runthrough, write_header
from sampling_qol import ObsParameter, Source
from sigpyproc.readers import FilReader

# datos para prueba
vela_pulsar = Source(
    "J08354510", 083520.6, -451034.8
)  # creamos el objeto porque solo andamos mirando vela

# pequeña sección de tiempo porque estoy 60% seguro que la fecha de modificación de
# los archivos de drive será cuando los descargué, así que lo estoy haciendo a mano con el historial de drive
# y un poco de inferencia
mjd1 = 60434
tst_fold04 = ObsParameter(2.048e6, mjd1, 400, vela_pulsar, "tstFold04.bin", 1)
tst_fold04.set_channels(32)
tst_fold04.header_data()
write_header(tst_fold04)
file_runthrough(tst_fold04.file, "tstFold04.fil")

outfilterbank = FilReader("tstFold04.fil")  # lee el filterbank
print(
    outfilterbank.header.fch1
)  # esto solo printea el fch1 (la frecuencia central del primer canal)
