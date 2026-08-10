class Source(object):
    def __init__(self, name: str, right_ascension: float, declination: float):
        self.name = name
        self.ra = right_ascension
        self.dec = declination


VELA_PULSAR = Source(
    "J08354510", 083520.6, -451034.8
)  # creamos el objeto porque solo andamos mirando vela


class ObsParameter(object):
    def __init__( 
        self,
        sample_rate: float,
        obstime: float,
        center_frequency: float,
        rawfile: str,
        sdr: int,
        source: Source = VELA_PULSAR,
    ):
        self.sample_rate = sample_rate
        self.obstime = obstime
        self.source = source.name
        self.ra = source.ra
        self.dec = source.dec
        self.cfreq = center_frequency
        self.tsample = (1 / self.sample_rate) * 32 * 20 # 32 canales y 20 muestras calculadas en el paso 3 de la canalización de Hawkrao
        self.file = rawfile
        self.sdr = sdr

    def set_channels(self, channels):
        self.channels = channels

    def header_data(self):
        try:
            self.channel_width = -(self.sample_rate / self.channels) * 1e-6
            self.fch1 = (
                self.cfreq + (self.sample_rate / 2.0) * 1e-6 + (self.channel_width * 0.5)
            )

        except NameError:
            print(
                "No se ha especificado la cantidad de canales, defínala e intente de nuevo."
            )


def load_data(file):
    with open(file) as data:
        data.readline()
        lines = []
        for i in range(3):
            lines.append(float(data.readline().rstrip("\n").split(",")[-1]))
        lines.append(data.readline().rstrip("\n").split(",")[-1])
        lines.append(int(data.readline().rstrip("\n").split(",")[-1]) // 8)
    return ObsParameter(*lines)


def take_samples():
    pass
