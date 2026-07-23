from main import file_runthrough, write_header
from sampling_qol import load_data
from sigpyproc.readers import FilReader

# datos para prueba


# pequeña sección de tiempo porque estoy 60% seguro que la fecha de modificación de
# los archivos de drive será cuando los descargué, así que lo estoy haciendo a mano con el historial de drive
# y un poco de inferencia
tst_fold04 = load_data("tstFold04.params")
tst_fold04.set_channels(32)
tst_fold04.header_data()
write_header(tst_fold04)
file_runthrough(tst_fold04.file, "tstFold04.fil")

outfilterbank = FilReader("tstFold04.fil")  # lee el filterbank
print(
    outfilterbank.header
)  # esto solo printea el fch1 (la frecuencia central del primer canal)
