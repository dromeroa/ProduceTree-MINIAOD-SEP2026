import os
import glob
import subprocess
import time
import psutil
from multiprocessing import Pool

# --- CONFIGURACION ---
FOLDER_PATH     = '/code/Signal-1000'
OUTPUT_DIR      = "output_signal-1000_chunks"
CORES           = 14
MAX_RAM_PERCENT = 70.0
FILES_PER_JOB   = 4

# Heredar el entorno CMSSW completo (cmsRun necesita CMSSW_BASE, ROOT, etc.)
CMSSW_ENV = os.environ.copy()


def run_process(args):
    job_id, files_chunk = args
    if not files_chunk:
        return "Tarea {} ignorada: chunk vacio.".format(job_id)

    # --- CONTROL DE RAM ---
    wait_attempts = 0
    while psutil.virtual_memory().percent > MAX_RAM_PERCENT:
        wait_attempts += 1
        if wait_attempts > 20:
            return "Tarea {} abortada: RAM sostenidamente alta ({}%).".format(
                job_id, round(psutil.virtual_memory().percent, 1))
        print("[Tarea {}] RAM al {}%. Esperando 15s...".format(
            job_id, round(psutil.virtual_memory().percent, 1)))
        time.sleep(15)

    # --- RUTAS DE SALIDA ---
    # sonic_cfg.py usa TFileService con parametro outputCSV (no outputFile/ROOT)
    files_arg      = ",".join(files_chunk)
    output_csv     = os.path.join(OUTPUT_DIR, "output_{}.csv".format(job_id))
    log_filename   = os.path.join(OUTPUT_DIR, "log_{}.txt".format(job_id))

    # Asegurar que OUTPUT_DIR existe (Python 2.7 compatible)
    if not os.path.exists(OUTPUT_DIR):
        try:
            os.makedirs(OUTPUT_DIR)
        except OSError:
            pass  # race condition: otro proceso ya lo creo

    cmd = [
        "cmsRun", "sonic_cfg.py",
        "inputFiles={}".format(files_arg),
        "outputCSV={}".format(output_csv),   # <-- parametro correcto
    ]

    print("[Tarea {}] Iniciando con {} archivos. RAM: {}%".format(
        job_id, len(files_chunk), round(psutil.virtual_memory().percent, 1)))

    try:
        with open(log_filename, "w") as logfile:
            ret = subprocess.call(
                cmd,
                stdout=logfile,
                stderr=subprocess.STDOUT,
                env=CMSSW_ENV
            )

        # --- Verificar exit code ---
        if ret != 0:
            msg = "[Tarea {}] FALLO cmsRun (exit code {}). Ver: {}".format(
                job_id, ret, log_filename)
            print(msg)
            return msg

        # --- Verificar que el CSV existe y tiene contenido ---
        if not os.path.exists(output_csv):
            msg = "[Tarea {}] FALLO: cmsRun OK pero no creo {}".format(
                job_id, output_csv)
            print(msg)
            return msg

        size_kb = os.path.getsize(output_csv) / 1024.0
        if size_kb < 0.1:
            msg = "[Tarea {}] ADVERTENCIA: CSV vacio o casi vacio ({:.2f} KB): {}".format(
                job_id, size_kb, output_csv)
            print(msg)
            return msg

        msg = "[Tarea {}] OK -> {} ({:.1f} KB). RAM: {}%".format(
            job_id, output_csv, size_kb, round(psutil.virtual_memory().percent, 1))
        print(msg)
        return msg

    except Exception as e:
        msg = "[Tarea {}] EXCEPCION: {}".format(job_id, str(e))
        print(msg)
        return msg


if __name__ == '__main__':

    # --- Verificacion temprana del entorno ---
    try:
        proc = subprocess.Popen(
            ["which", "cmsRun"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout, _ = proc.communicate()
        cmsrun_path = stdout.decode("utf-8").strip()
    except Exception:
        cmsrun_path = ""

    if not cmsrun_path:
        print("FATAL: 'cmsRun' no esta en PATH. Ejecuta primero: cmsenv")
        exit(1)

    if not os.path.exists("sonic_cfg.py"):
        print("FATAL: 'sonic_cfg.py' no encontrado en: {}".format(os.getcwd()))
        exit(1)

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # --- Buscar archivos ROOT de entrada ---
    raw_files = sorted(glob.glob(os.path.join(FOLDER_PATH, "*.root")))
    all_files = ["file:{}".format(f) for f in raw_files]

    if not all_files:
        print("FATAL: No se encontraron archivos .root en {}".format(FOLDER_PATH))
        exit(1)

    # --- Crear chunks ---
    pool_args = []
    for i in range(0, len(all_files), FILES_PER_JOB):
        chunk  = all_files[i : i + FILES_PER_JOB]
        job_id = i // FILES_PER_JOB
        pool_args.append((job_id, chunk))

    print("\n--- INFO ---")
    print("Archivos .root encontrados   : {}".format(len(all_files)))
    print("Tareas totales               : {}".format(len(pool_args)))
    print("Procesos simultaneos (Cores) : {}".format(CORES))
    print("Limite de RAM                : {}%".format(MAX_RAM_PERCENT))
    print("Directorio de salida         : {}".format(os.path.abspath(OUTPUT_DIR)))
    print("Salida por tarea             : {}/output_N.csv".format(OUTPUT_DIR))
    print("cmsRun                       : {}".format(cmsrun_path))
    print("------------\n")

    # --- Ejecucion con Pool ---
    pool = Pool(processes=CORES)
    try:
        results = pool.map(run_process, pool_args)
    except KeyboardInterrupt:
        print("\nInterrumpido. Cerrando pool...")
        pool.terminate()
        pool.join()
        exit(1)
    pool.close()
    pool.join()

    # --- Resumen ---
    print("\n========== RESUMEN ==========")
    ok     = [r for r in results if r and "OK"        in r]
    failed = [r for r in results if r and ("FALLO"    in r
                                        or "ERROR"    in r
                                        or "EXCEPCION" in r
                                        or "ADVERTENCIA" in r)]
    print("Tareas exitosas  : {}".format(len(ok)))
    print("Tareas con error : {}".format(len(failed)))
    if failed:
        print("\nDetalle de errores:")
        for f in failed:
            print("  {}".format(f))

    # --- Fusion de CSVs ---
    print("\nVerificando archivos para fusion...")
    csv_files = sorted(glob.glob("{}/output_*.csv".format(OUTPUT_DIR)))
    print("CSVs generados: {}".format(len(csv_files)))

    if csv_files:
        output_merged = "data_completa.csv"
        print("Fusionando {} CSVs -> {}".format(len(csv_files), output_merged))
        header_written = False
        rows_total = 0
        try:
            with open(output_merged, "w") as fout:
                for csv_path in csv_files:
                    with open(csv_path, "r") as fin:
                        lines = fin.readlines()
                    if not lines:
                        print("  ADVERTENCIA: {} esta vacio, saltando.".format(csv_path))
                        continue
                    if not header_written:
                        # Escribir cabecera solo del primer archivo
                        fout.writelines(lines)
                        header_written = True
                        rows_total += max(0, len(lines) - 1)
                    else:
                        # Saltar la primera linea (cabecera) de los demas
                        fout.writelines(lines[1:])
                        rows_total += max(0, len(lines) - 1)

            size_mb = os.path.getsize(output_merged) / (1024.0 * 1024.0)
            print("Fusion completada: {} ({:.1f} MB, ~{} filas de datos)".format(
                output_merged, size_mb, rows_total))
        except Exception as e:
            print("ERROR durante la fusion: {}".format(e))
    else:
        print("Error: No se generaron CSVs de salida.")
        logs = sorted(glob.glob("{}/log_*.txt".format(OUTPUT_DIR)))
        if logs:
            print("\nPrimeras 40 lineas de {}:".format(logs[0]))
            print("-" * 50)
            with open(logs[0]) as f:
                for i, line in enumerate(f):
                    if i >= 40: break
                    print(line.rstrip())
