import sys
import subprocess
from pathlib import Path

ROOT_DIR = Path()
OUTPUT_FILE_DIR = None


def clear_input_buffer():
    if sys.platform == 'win32':
        import msvcrt
        while msvcrt.kbhit():
            msvcrt.getch()
    else:
        import termios
        termios.tcflush(sys.stdin, termios.TCIFLUSH)


def list_benchmark_tests(root_dir: Path, output_file_dir: Path) -> None:
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "--collect-only", "-m", "benchmark", "-q", "--disable-warnings"],
        cwd=root_dir,
        capture_output=True,
        text=True
    )

    # Combine stdout and stderr, then filter for test lines
    output = result.stdout + result.stderr
    lines = output.strip().split('\n')
    test_functions = []

    for line in lines:
        if '::' in line and line.strip().startswith('tests'):
            test_functions.append(line.strip())

    if test_functions:
        print("Available benchmarks:")
        print("=" * 70)
        print(" 0 - Terminate task")
        for i, test_func in enumerate(test_functions):
            print(f" {i + 1} - {test_func}")
        print("=" * 70)

        clear_input_buffer()

        try:
            selection = int(input("Enter the index of the benchmark to run: "))
            if selection == 0:
                print("\nTask terminated.")
                sys.exit(1)
            elif 0 < selection <= len(test_functions):
                selected_test = test_functions[selection - 1]

                output_file_dir.mkdir(parents=True, exist_ok=True)
                output_file = output_file_dir / "test_to_profile.txt"
                with open(output_file, 'w') as f:
                    f.write(selected_test)

                print(f"\nSelected: {selected_test}")
                print(f"Written to: {output_file}")
                sys.exit(0)
            else:
                print(f"\nInvalid index. Please enter a number between 0 and {len(test_functions)}")
                sys.exit(1)
        except ValueError:
            print("\nInvalid input. Please enter a valid integer.")
            sys.exit(1)

    else:
        sys.exit(1)


if __name__ == "__main__":

    if len(sys.argv) > 2:
        ROOT_DIR = Path(sys.argv[1]).resolve()
        OUTPUT_FILE_DIR = Path(sys.argv[2]).resolve()
    else:
        raise ValueError("Missing arguments. Usage: list_benchmarks.py <project_root> <output_file>")

    list_benchmark_tests(ROOT_DIR, OUTPUT_FILE_DIR)
