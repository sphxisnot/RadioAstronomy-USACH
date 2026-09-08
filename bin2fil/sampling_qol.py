import time
import subprocess
import numpy as np
import glob

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
        dm: float,
        pw50: float,
        sdr: int,
        rawfile: str,
        source: Source = VELA_PULSAR
    ):
        self.sample_rate = sample_rate
        self.obstime = obstime
        self.source = source.name
        self.ra = source.ra
        self.dec = source.dec
        self.cfreq = center_frequency
        self.file = glob.glob("**/"+rawfile,recursive=True)[0]
        self.sdr = sdr
        self.dm = dm
        self.pw50 = pw50
        self.delay = self.dm*8.3e3*sample_rate*1e-3/(self.cfreq**3)
        self.channels = int(np.max([np.pow(2,np.ceil(np.log2((4*self.delay)/self.pw50))),32]))
        self.channel_width = -(self.sample_rate / self.channels) * 1e-6
        self.fch1 = (
            self.cfreq + (self.sample_rate *0.5) * 1e-6 + (self.channel_width * 0.5)
        )
        self.tsample = (1 / self.sample_rate) * self.channels * 20
        self.outfile = str(rawfile).removesuffix(".iq").removesuffix(".bin")+".fil"


def load_data(file):
    with open(file,"r") as data:
        data.readline()
        lines = []
        for i in range(5):
            lines.append(float(data.readline().rstrip("\n").split(",")[-1]))
        lines.append(int(data.readline().rstrip("\n").split(",")[-1]))
        lines.append(data.readline().rstrip("\n").split(",")[-1])
        source=data.readline().rstrip("\n").split(",")[1:4]
        pulsar=Source(source[0],float(source[1]),float(source[2]))
        lines.append(pulsar)
    return ObsParameter(*lines)


def write_obs_data(outfile, sample_rate, time, tuned_freq, datafile, data_type,source):
    with open(outfile, "+w") as params:
        params.write("Name [Unit], Value")
        params.write(f"Sample Rate [Samples/Second], {sample_rate}")
        params.write(f"Obs Time [MJD], {time}")
        params.write(f"Tuned Frequency [Mhz], {tuned_freq * 1e-6}")
        params.write(f"Raw Data File, {datafile}")
        params.write(f"Data Type [bits], {data_type}")
        params.write(f"Source, {source}")

def live_sampling(sample_rate,tuned_frequency,integration_time,outfile):
    num_samples= integration_time*sample_rate
    cli_command = f"airspy_rx -r {outfile}.bin -a {sample_rate} -f {tuned_frequency} -n {num_samples}"
    subprocess.run(cli_command)

def listen_from_file(obsfile):
    with open(obsfile,"r+") as file:
        file.readline()
        obsparams=[]
    live_sampling(*obsparams)
