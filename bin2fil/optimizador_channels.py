"""
optimizador_channels.py

Calcula el numero optimo de canales de FFT para una observacion de pulsar.

  - MAS canales -> menos smearing por dispersion dentro de cada canal
        t_DM = 8.3e3 * DM * (B/nchan) / nu^3     [s]   (baja con nchan)

  - MAS canales -> peor resolucion temporal
        tsamp = nchan * group_chunks / sample_rate [s] (sube con nchan)

Como uno baja y el otro sube, el ancho efectivo del pulso
(suma en cuadratura) tiene un minimo. Ocurre cuando t_DM = tsamp:

        nchan_opt = sqrt( 8.3e3 * DM * B * sample_rate / (nu^3 * N) )

donde B y nu van en MHz, sample_rate en muestras/s, N = group_chunks.

Uso:
    python optimizador_channels.py --dm 67.77 --sample-rate 10e6 --cfreq 400
    python optimizador_channels.py --dm 67.77 --sample-rate 10e6 --cfreq 400 \\
        --group-chunks 4 --pulse-width 1.0
"""

import argparse
import math


def dm_smearing_s(dm, bandwidth_MHz, nu_MHz):
    """Retraso por dispersion a traves de un ancho de banda dado, en segundos."""
    return 8.3e3 * dm * bandwidth_MHz / nu_MHz**3


def optimal_channels(dm, sample_rate, cfreq_MHz, group_chunks=20):
    """
    Numero de canales que minimiza el ancho efectivo del pulso.

    Devuelve el valor exacto (float). En la practica conviene redondear a
    la potencia de 2 mas cercana, porque la FFT es mas eficiente asi.
    """
    B = sample_rate / 1e6  # ancho de banda en MHz = sample rate en MHz
    n2 = 8.3e3 * dm * B * sample_rate / (cfreq_MHz**3 * group_chunks)
    return math.sqrt(n2)


def nearest_power_of_2(n, round_up=False):

    if n <= 0:
        raise ValueError(
            f"El numero de canales debe ser positivo, se recibio {n!r}. "
            "Un n_opt <= 0 suele significar dm=0 o parametros mal cargados."
        )

    lo = 2**math.floor(math.log2(n))
    up = 2**math.ceil(math.log2(n))

    chosen = up if round_up else (lo if n < math.sqrt(lo * up) else up)

    return max(1, int(chosen))


def evaluate(nchan, dm, sample_rate, cfreq_MHz, group_chunks, total_samples=None):
    """Evalua una configuracion concreta y devuelve sus metricas."""
    B = sample_rate/1e6
    chan_width = B/nchan
    t_dm = dm_smearing_s(dm, chan_width, cfreq_MHz)
    tsamp = nchan*group_chunks/sample_rate
    effective = math.sqrt(t_dm**2 + tsamp**2)

    out = {
        "nchan": nchan,
        "chan_width_MHz": chan_width,
        "t_dm_ms": t_dm * 1e3,
        "tsamp_ms": tsamp * 1e3,
        "effective_ms": effective * 1e3,
    }
    if total_samples:
        rows = total_samples / (nchan * group_chunks)
        out["rows"] = rows
        out["fil_size_GB"] = rows*nchan*4/1e9
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dm", type=float, required=True,
                    help="Dispersion Measure en pc/cm3 (Vela: 67.77)")
    ap.add_argument("--sample-rate", type=float, required=True,
                    help="Muestras complejas por segundo (ej. 10e6)")
    ap.add_argument("--cfreq", type=float, required=True,
                    help="Frecuencia central en MHz (ej. 400)")

    ap.add_argument("--group-chunks", type=int, default=4,
                    help="Espectros integrados por fila de salida "
                         "(default 4, igual que GROUP_CHUNKS en test.py)")
    ap.add_argument("--pulse-width", type=float, default=1.0,
                    help="Ancho del pulso en ms (Vela: ~1.0)")
    ap.add_argument("--file-size-gb", type=float, default=None,
                    help="Tamano del .iq en GB, para estimar el .fil resultante")
    ap.add_argument("--bits", type=int, choices=[8, 16], default=16,
                    help="Bits por componente I/Q del archivo de entrada")
    args = ap.parse_args()

    B = args.sample_rate/1e6
    total_samples = None
    if args.file_size_gb:
        bytes_per_complex = 2*(args.bits // 8)
        total_samples = args.file_size_gb*1e9/bytes_per_complex

    print("="*66)
    print("CALCULO DE CANALES OPTIMOS")
    print("="*66)
    print(f"  DM                : {args.dm} pc/cm3")
    print(f"  Sample rate       : {args.sample_rate:.3e} muestras/s")
    print(f"  Ancho de banda    : {B} MHz")
    print(f"  Frecuencia central: {args.cfreq} MHz")
    print(f"  group_chunks (N)  : {args.group_chunks}")
    print(f"  Ancho del pulso   : {args.pulse_width} ms")

    total_smear = dm_smearing_s(args.dm, B, args.cfreq)*1e3
    print(f"\n  Dispersion a traves de la banda completa: {total_smear:.1f} ms")
    if total_smear > 50:
        print("  -> Muy grande: la dedispersion es imprescindible, no opcional.")

    n_opt = optimal_channels(args.dm, args.sample_rate, args.cfreq, args.group_chunks)
    n_pow2 = nearest_power_of_2(n_opt)
    print(f"\n  Optimo teorico    : {n_opt:.1f} canales")
    print(f"  Potencia de 2      : {n_pow2} canales  <-- RECOMENDADO")

    print("\n" + "="*66)
    print("COMPARACION DE CONFIGURACIONES")
    print("="*66)
    header = f"{'Canales':>8} {'Ancho ch':>10} {'Smear DM':>10} {'tsamp':>9} {'Efectivo':>10}"
    if total_samples:
        header += f" {'.fil':>9}"
    print(header)
    print(f"{'':>8} {'(MHz)':>10} {'(ms)':>10} {'(ms)':>9} {'(ms)':>10}" +
          (f" {'(GB)':>9}" if total_samples else ""))
    print("-"*len(header))

    candidates = sorted({32, 64, 128, 256, 512, 1024, n_pow2})
    for n in candidates:
        m = evaluate(n, args.dm, args.sample_rate, args.cfreq,
                     args.group_chunks, total_samples)
        mark = " <--" if n == n_pow2 else ""
        ok = "" if m["effective_ms"] < args.pulse_width else "  (!)"
        line = (f"{n:>8} {m['chan_width_MHz']:>10.4f} {m['t_dm_ms']:>10.3f} "
                f"{m['tsamp_ms']:>9.3f} {m['effective_ms']:>10.3f}")
        if total_samples:
            line += f" {m['fil_size_GB']:>9.2f}"
        print(line + mark + ok)

    print("\n  (!) = ancho efectivo mayor que el pulso: se pierde resolucion")

    # Efecto de group_chunks
    print("\n" + "="*66)
    print("EFECTO DE group_chunks (N)")
    print("="*66)
    print("Bajar N mejora la resolucion temporal y permite mas canales,")
    print("a costa de un .fil mas grande (menos integracion).\n")
    print(f"{'N':>4} {'n_opt':>8} {'pot.2':>7} {'Efectivo (ms)':>15}" +
          (f" {'.fil (GB)':>11}" if total_samples else ""))
    print("-"*50)
    for N in [1, 2, 4, 10, 20, 40]:
        no = optimal_channels(args.dm, args.sample_rate, args.cfreq, N)
        np2 = nearest_power_of_2(no)
        m = evaluate(np2, args.dm, args.sample_rate, args.cfreq, N, total_samples)
        line = f"{N:>4} {no:>8.1f} {np2:>7} {m['effective_ms']:>15.3f}"
        if total_samples:
            line += f" {m['fil_size_GB']:>11.2f}"
        print(line)

    print("\n" + "="*66)
    print("NOTA")
    print("="*66)
    print("Este calculo optimiza la RESOLUCION del perfil, no la sensibilidad.")
    print("Si el sistema no tiene suficiente ganancia o tiempo de integracion,")
    print("mas canales no van a producir una deteccion. El tiempo de")
    print("observacion importa: la SNR crece como sqrt(tiempo).")


if __name__ == "__main__":
    main()