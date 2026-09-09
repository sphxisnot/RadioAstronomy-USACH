"""
Interfaz grafica (Tkinter, libreria estandar de Python, no requiere instalar
nada extra) para ingresar los parametros de la observacion a mano y correr
el pipeline de conversion I/Q -> filterbank sin tener que editar archivos
.params ni el codigo cada vez.

Uso:
    python pulsar_gui.py

Requiere que main.py, pulsar_main.py, convertidor.py, sampling_qol.py y
optimizador_channels.py esten en la misma carpeta (o en el PYTHONPATH),
tal como ya lo requiere test.py.
"""
import os
import queue
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext

from sampling_qol import ObsParameter, Source
from pulsar_main import write_header
from convertidor import convert_and_write
from optimizador_channels import optimal_channels, nearest_power_of_2, evaluate


# ------------------------------------------------------------------
# Lectura / escritura de .params compatible con sampling_qol.load_data
# ------------------------------------------------------------------
def load_params_into(fields, path):
    """Llena el diccionario `fields` (StringVars) leyendo un .params
    con el mismo formato que espera sampling_qol.load_data."""
    with open(path, "r") as f:
        f.readline()  # linea de encabezado, se descarta
        sample_rate = f.readline().rstrip("\n").split(",")[-1]
        obstime = f.readline().rstrip("\n").split(",")[-1]
        cfreq = f.readline().rstrip("\n").split(",")[-1]
        rawfile = f.readline().rstrip("\n").split(",")[-1]
        data_type_bits = f.readline().rstrip("\n").split(",")[-1]
        source_parts = f.readline().rstrip("\n").split(",")[1:4]

    fields["sample_rate"].set(sample_rate.strip())
    fields["obstime"].set(obstime.strip())
    fields["cfreq"].set(cfreq.strip())
    fields["rawfile"].set(rawfile.strip())
    fields["data_type_bits"].set(data_type_bits.strip())
    fields["source_name"].set(source_parts[0].strip())
    fields["ra"].set(source_parts[1].strip())
    fields["dec"].set(source_parts[2].strip())


def save_params_from(fields, path):
    """Guarda un .params con el mismo formato que lee sampling_qol.load_data,
    a partir de los valores actuales de los campos de la GUI."""
    with open(path, "w") as f:
        f.write("Name [Unit] , Value\n")
        f.write(f"Sample Rate [Samples/Second],{fields['sample_rate'].get()}\n")
        f.write(f"Obs Time [MJD],{fields['obstime'].get()}\n")
        f.write(f"Tuned Frequency [MHz],{fields['cfreq'].get()}\n")
        f.write(f"Raw Data File,{fields['rawfile'].get()}\n")
        f.write(f"Data Type,{fields['data_type_bits'].get()}\n")
        f.write(
            f"Source, {fields['source_name'].get()}, "
            f"{fields['ra'].get()}, {fields['dec'].get()}\n"
        )


# ------------------------------------------------------------------
# GUI
# ------------------------------------------------------------------
class PulsarGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Pulsar - Conversion I/Q a Filterbank")
        self.geometry("640x760")
        self.resizable(False, False)

        self.log_queue = queue.Queue()
        self.worker = None

        self.fields = {
            "sample_rate": tk.StringVar(value="10.0e6"),
            "obstime": tk.StringVar(value="60892.91"),
            "cfreq": tk.StringVar(value="400"),
            "rawfile": tk.StringVar(value="400psr_01.iq"),
            "data_type_bits": tk.StringVar(value="16"),
            "source_name": tk.StringVar(value="J08354510"),
            "ra": tk.StringVar(value="083520.6"),
            "dec": tk.StringVar(value="-451034.8"),
            "dm": tk.StringVar(value="67.77"),
            "pulse_width_ms": tk.StringVar(value="1.0"),
            "group_chunks": tk.StringVar(value="4"),
            "channels": tk.StringVar(value=""),  # vacio = automatico
        }

        self._build_layout()
        self.after(100, self._poll_log_queue)

    # ---------------- layout ----------------
    def _build_layout(self):
        pad = {"padx": 8, "pady": 4}

        top = ttk.Frame(self)
        top.pack(fill="x", **pad)
        ttk.Button(top, text="Cargar .params...", command=self._on_load_params).pack(
            side="left", padx=4
        )
        ttk.Button(top, text="Guardar como .params...", command=self._on_save_params).pack(
            side="left", padx=4
        )

        obs_frame = ttk.LabelFrame(self, text="Parametros de la observacion")
        obs_frame.pack(fill="x", **pad)

        self._row(obs_frame, 0, "Archivo de datos (.iq / .bin)", "rawfile", browse=True)
        self._row(obs_frame, 1, "Sample Rate [muestras/s]", "sample_rate")
        self._row(obs_frame, 2, "Obs Time [MJD]", "obstime")
        self._row(obs_frame, 3, "Frecuencia central [MHz]", "cfreq")

        ttk.Label(obs_frame, text="Tipo de dato").grid(row=4, column=0, sticky="w", **pad)
        dt_frame = ttk.Frame(obs_frame)
        dt_frame.grid(row=4, column=1, sticky="w")
        ttk.Radiobutton(
            dt_frame, text="8 bits (RTL-SDR)", value="8",
            variable=self.fields["data_type_bits"],
        ).pack(side="left")
        ttk.Radiobutton(
            dt_frame, text="16 bits (AirSpy)", value="16",
            variable=self.fields["data_type_bits"],
        ).pack(side="left")

        self._row(obs_frame, 5, "Nombre fuente", "source_name")
        self._row(obs_frame, 6, "RA (src_raj)", "ra")
        self._row(obs_frame, 7, "Dec (src_dej)", "dec")

        proc_frame = ttk.LabelFrame(self, text="Parametros de procesamiento")
        proc_frame.pack(fill="x", **pad)

        self._row(proc_frame, 0, "DM de la fuente [pc/cm^3]", "dm")
        self._row(proc_frame, 1, "Ancho del pulso esperado [ms]", "pulse_width_ms")
        self._row(proc_frame, 2, "Group chunks (integracion)", "group_chunks")
        self._row(
            proc_frame, 3,
            "Canales (vacio = automatico segun DM)", "channels",
        )

        ttk.Button(
            self, text="Evaluar parametros", command=self._on_evaluate
        ).pack(pady=4)
        ttk.Button(
            self, text="Procesar (generar .fil)", command=self._on_process
        ).pack(pady=4)

        log_frame = ttk.LabelFrame(self, text="Registro")
        log_frame.pack(fill="both", expand=True, **pad)
        self.log_widget = scrolledtext.ScrolledText(log_frame, height=18, state="disabled")
        self.log_widget.pack(fill="both", expand=True)

    def _row(self, parent, r, label, key, browse=False):
        ttk.Label(parent, text=label).grid(row=r, column=0, sticky="w", padx=8, pady=4)
        entry = ttk.Entry(parent, textvariable=self.fields[key], width=40)
        entry.grid(row=r, column=1, sticky="w", padx=4)
        if browse:
            ttk.Button(
                parent, text="Buscar...", command=lambda: self._browse_rawfile(key)
            ).grid(row=r, column=2, padx=4)

    def _browse_rawfile(self, key):
        path = filedialog.askopenfilename(
            title="Seleccionar archivo de datos crudos",
            filetypes=[("Datos IQ/BIN", "*.iq *.bin"), ("Todos los archivos", "*.*")],
        )
        if path:
            self.fields[key].set(path)

    # ---------------- log helpers ----------------
    def _log(self, msg):
        self.log_queue.put(msg)

    def _poll_log_queue(self):
        try:
            while True:
                msg = self.log_queue.get_nowait()
                self.log_widget.configure(state="normal")
                self.log_widget.insert("end", msg + "\n")
                self.log_widget.see("end")
                self.log_widget.configure(state="disabled")
        except queue.Empty:
            pass
        self.after(100, self._poll_log_queue)

    # ---------------- acciones ----------------
    def _on_load_params(self):
        path = filedialog.askopenfilename(
            title="Cargar .params", filetypes=[("Params", "*.params"), ("Todos", "*.*")]
        )
        if not path:
            return
        try:
            load_params_into(self.fields, path)
            self._log(f"Parametros cargados desde: {path}")
        except Exception as e:
            messagebox.showerror("Error al cargar", str(e))

    def _on_save_params(self):
        path = filedialog.asksaveasfilename(
            title="Guardar .params", defaultextension=".params",
            filetypes=[("Params", "*.params")],
        )
        if not path:
            return
        try:
            save_params_from(self.fields, path)
            self._log(f"Parametros guardados en: {path}")
        except Exception as e:
            messagebox.showerror("Error al guardar", str(e))

    def _read_form(self):
        """Valida y castea los campos del formulario. Lanza ValueError con
        un mensaje entendible si algo no es valido."""
        try:
            sample_rate = float(self.fields["sample_rate"].get())
            obstime = float(self.fields["obstime"].get())
            cfreq = float(self.fields["cfreq"].get())
            rawfile = self.fields["rawfile"].get().strip()
            data_type_bits = int(self.fields["data_type_bits"].get())
            source_name = self.fields["source_name"].get().strip()
            ra = float(self.fields["ra"].get())
            dec = float(self.fields["dec"].get())
            dm = float(self.fields["dm"].get())
            pulse_width_ms = float(self.fields["pulse_width_ms"].get())
            group_chunks = int(self.fields["group_chunks"].get())
        except ValueError as e:
            raise ValueError(f"Revisa los campos numericos: {e}")

        if not rawfile:
            raise ValueError("Debes indicar el archivo de datos crudos.")
        if not os.path.isfile(rawfile):
            raise ValueError(f"No se encuentra el archivo de datos: {rawfile}")
        if data_type_bits not in (8, 16):
            raise ValueError("El tipo de dato debe ser 8 o 16 bits.")

        channels_str = self.fields["channels"].get().strip()
        forced_channels = int(channels_str) if channels_str else None

        sdr = data_type_bits // 8  # 1 -> 8 bits, 2 -> 16 bits (igual que load_data)
        source = Source(source_name, ra, dec)
        obs = ObsParameter(
            sample_rate=sample_rate,
            obstime=obstime,
            center_frequency=cfreq,
            rawfile=rawfile,
            sdr=sdr,
            source=source,
        )
        return obs, dm, pulse_width_ms, group_chunks, forced_channels, data_type_bits

    def _decide_channels(self, obs, dm, group_chunks, forced_channels):
        if forced_channels is not None:
            channels = forced_channels
            self._log(f"Canales fijados a mano: {channels}")
        else:
            n_opt = optimal_channels(dm, obs.sample_rate, obs.cfreq, group_chunks)
            channels = nearest_power_of_2(n_opt)
            self._log(f"Canales calculados automaticamente para DM={dm}:")
            self._log(f"  optimo teorico: {n_opt:.1f}  ->  usando {channels}")
        return channels

    def _on_evaluate(self):
        try:
            obs, dm, pulse_width_ms, group_chunks, forced_channels, _ = self._read_form()
        except ValueError as e:
            messagebox.showerror("Parametros invalidos", str(e))
            return

        channels = self._decide_channels(obs, dm, group_chunks, forced_channels)
        m = evaluate(channels, dm, obs.sample_rate, obs.cfreq, group_chunks)
        self._log(f"  Ancho de canal    : {m['chan_width_MHz']:.4f} MHz")
        self._log(f"  Smearing intra-ch : {m['t_dm_ms']:.3f} ms")
        self._log(f"  Resolucion (tsamp): {m['tsamp_ms']:.3f} ms")
        efectivo = f"  Ancho efectivo    : {m['effective_ms']:.3f} ms"
        if m["effective_ms"] < pulse_width_ms:
            self._log(efectivo + f"  (< {pulse_width_ms} ms del pulso: OK)")
        else:
            self._log(efectivo + f"  (!) mayor que el pulso de {pulse_width_ms} ms")
            self._log("      El pulso saldra ensanchado. Considerar bajar group chunks.")

    def _on_process(self):
        if self.worker and self.worker.is_alive():
            messagebox.showwarning("En proceso", "Ya hay una conversion en curso.")
            return
        try:
            form = self._read_form()
        except ValueError as e:
            messagebox.showerror("Parametros invalidos", str(e))
            return

        self.worker = threading.Thread(target=self._process_worker, args=form, daemon=True)
        self.worker.start()

    def _process_worker(self, obs, dm, pulse_width_ms, group_chunks, forced_channels, data_type_bits):
        try:
            channels = self._decide_channels(obs, dm, group_chunks, forced_channels)

            m = evaluate(channels, dm, obs.sample_rate, obs.cfreq, group_chunks)
            self._log(f"  Ancho de canal    : {m['chan_width_MHz']:.4f} MHz")
            self._log(f"  Smearing intra-ch : {m['t_dm_ms']:.3f} ms")
            self._log(f"  Resolucion (tsamp): {m['tsamp_ms']:.3f} ms")
            if m["effective_ms"] >= pulse_width_ms:
                self._log(
                    f"  (!) Ancho efectivo {m['effective_ms']:.3f} ms es mayor "
                    f"que el pulso de {pulse_width_ms} ms. El pulso saldra ensanchado."
                )

            obs.set_channels(channels)
            obs.header_data()

            outfil = str(obs.file).removesuffix(".iq").removesuffix(".bin") + ".fil"
            if os.path.exists(outfil):
                self._log(f"Aviso: {outfil} ya existia, se agregaran datos al final (modo 'ab').")

            self._log(f"Escribiendo header en: {outfil}")
            write_header(obs, group_chunks=group_chunks)

            self._log("Procesando archivo (esto puede tardar segun el tamano)...")
            n_rows, tsamp = convert_and_write(
                obs.file, outfil,
                sample_rate=obs.sample_rate,
                channels=channels,
                group_chunks=group_chunks,
                dtype_bits=data_type_bits,
            )

            self._log(f"\nGenerado: {outfil}")
            self._log(f"  {n_rows:,} filas x {channels} canales")
            self._log(f"  tsamp    = {tsamp * 1e6:.1f} us")
            self._log(f"  Duracion = {n_rows * tsamp:.1f} s")
            self._log(f"\nVerificar con: readfile {outfil}")
        except Exception as e:
            self._log(f"ERROR: {e}")
            messagebox.showerror("Error durante el procesamiento", str(e))


if __name__ == "__main__":
    app = PulsarGUI()
    app.mainloop()
