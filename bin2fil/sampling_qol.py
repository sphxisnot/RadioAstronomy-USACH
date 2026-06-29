class Source(object):
    def __init__(self, name: str, right_ascension: float, declination: float):
        self.name = name
        self.ra = right_ascension
        self.dec = declination


class ObsParameter(object):
    def __init__(
        self,
        sample_rate: float,
        obstime: float,
        center_frequency: float,
        source: Source,
        rawfile: str,
        sdr: int,
    ):
        self.sample_rate = sample_rate
        self.obstime = obstime
        self.source = source.name
        self.ra = source.ra
        self.dec = source.dec
        self.cfreq = center_frequency
        self.tsample = 1 / self.sample_rate
        self.file = rawfile
        self.sdr = sdr

    def set_channels(self, channels):
        self.channels = channels

    def header_data(self):
        try:
            self.channel_width = self.sample_rate / self.channels
            self.fch1 = (
                self.cfreq + (self.sample_rate * 2e-6) - (self.channel_width * 0.5)
            )

        except NameError:
            print(
                "No se ha especificado la cantidad de canales, defínala e intente de nuevo."
            )


def take_samples():
    pass
