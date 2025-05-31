import os
import shutil
import subprocess

SOURCE_DIR = os.path.join(os.path.dirname(__file__), 'SOURCE')
INPUT_DIR = os.path.join(SOURCE_DIR, '..', 'INPUT')
OUTPUT_DIR = os.path.join(SOURCE_DIR, '..', 'OUTPUT')
RESULTS_DIR = os.path.join(SOURCE_DIR, '..', 'RESULTS')

INPUT_FILE = os.path.join(INPUT_DIR, 'input.dat')
CONFIG_SPIN_SRC = os.path.join(OUTPUT_DIR, 'CONFIG', 'config.spin')
CONFIG_SPIN_DEST = os.path.join(INPUT_DIR, 'CONFIG', 'config.spin')
SIMULATOR_EXE = os.path.join(SOURCE_DIR, 'simulator.exe')

# Temperature list to iterate
TEMPERATURES = [2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 0.8, 0.7, 0.6, 0.5, 0.3]


def modify_input_file(temp, restart):
    with open(INPUT_FILE, 'r') as f:
        lines = f.readlines()

    with open(INPUT_FILE, 'w') as f:
        for line in lines:
            if line.startswith('TEMP'):
                f.write(f'TEMP                   {temp}\n')
            elif line.startswith('RESTART'):
                f.write(f'RESTART                {restart}\n')
            else:
                f.write(line)


def copy_config_spin():
    if os.path.exists(CONFIG_SPIN_SRC):
        shutil.copy(CONFIG_SPIN_SRC, CONFIG_SPIN_DEST)


def backup_output(temp):
    temp_str = f'{temp:.1f}'.replace('.', 'p')
    dir_name = f'out {temp:.1f} h=0 MRT'
    dest = os.path.join(RESULTS_DIR, temp_str)
    if os.path.exists(dest):
        shutil.rmtree(dest)
    shutil.copytree(OUTPUT_DIR, dest)


def build_simulator():
    print("Building simulator...")
    result = subprocess.run(['make'], cwd=SOURCE_DIR)
    if result.returncode != 0:
        raise RuntimeError("Build failed.")


def run_simulator():
    print("Running simulator...")
    result = subprocess.run([SIMULATOR_EXE], cwd=SOURCE_DIR)
    if result.returncode != 0:
        raise RuntimeError("Simulator run failed.")


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    build_simulator()

    for i, temp in enumerate(TEMPERATURES):
        restart_flag = 0 if i == 0 else 1
        print(f"\n=== Running simulation at T = {temp} (RESTART={restart_flag}) ===")

        # Set input file parameters
        modify_input_file(temp, restart_flag)

        # Copy previous config.spin if restarting
        if restart_flag == 1:
            copy_config_spin()

        # Run the simulation
        run_simulator()

        # Backup the OUTPUT directory
        backup_output(temp)

    print("\nAll simulations completed.")


if __name__ == '__main__':
    main()
