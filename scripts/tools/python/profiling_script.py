import sys
import cProfile
import pstats
import subprocess
import importlib
from pathlib import Path

ROOT_DIR = Path()
OUTPUT_PATH = Path()


def find_and_import_test(test_name: str, root_dir: Path):
    if "::" not in test_name:
        raise ValueError(f"Test name must be in format 'path/to/file.py::test_function', got: {test_name}")

    file_path_str, func_name = test_name.split("::", 1)
    file_path = root_dir / file_path_str

    if not file_path.exists():
        raise ValueError(f"Test file not found: {file_path}")

    relative_path = file_path.relative_to(root_dir)
    module_name = str(relative_path.with_suffix("")).replace("/", ".").replace("\\", ".")

    if str(root_dir) not in sys.path:
        sys.path.insert(0, str(root_dir))

    try:
        module = importlib.import_module(module_name)
        if hasattr(module, func_name):
            return getattr(module, func_name), module_name
        else:
            raise ValueError(f"Function '{func_name}' not found in module '{module_name}'")
    except ImportError as e:
        raise ValueError(f"Could not import module '{module_name}': {e}")


def create_mock_benchmark():
    class MockBenchmark:
        def __call__(self, func):
            # Run the function just once for profiling
            func()

    return MockBenchmark()


def profile_test(test_name: str, root_dir: Path, output_path: Path) -> None:
    output_path.mkdir(parents=True, exist_ok=True)
    sanitized_name = test_name.replace("::", "_").replace("/", "_").replace("\\", "_")
    profile_file = output_path / f"{sanitized_name}.prof"
    print(f"Profiling test: {test_name}")
    print(f"Output: {profile_file}")
    print()

    try:
        test_func, module_name = find_and_import_test(test_name, root_dir)
        print(f"Found test in: {module_name}")
    except ValueError as e:
        print(f"✗ Error: {e}")
        sys.exit(1)

    profiler = cProfile.Profile()
    mock_benchmark = create_mock_benchmark()

    print("\nRunning profiling...")

    try:
        profiler.enable()
        test_func(mock_benchmark)
        profiler.disable()

        profiler.dump_stats(str(profile_file))
        print(f"\n✓ Profile saved to: {profile_file}")

        stats = pstats.Stats(profiler)
        stats.strip_dirs()
        stats.sort_stats('cumulative')
        print("\nTop 10 functions by cumulative time:")
        print("=" * 70)
        stats.print_stats(10)

        print("\nLaunching snakeviz...")
        subprocess.run([sys.executable, "-m", "snakeviz", str(profile_file)])

    except Exception as e:
        profiler.disable()
        print(f"\n✗ Profiling failed: {e}")
        sys.exit(1)


if __name__ == "__main__":

    if len(sys.argv) < 4:
        print("Usage: python profiling_script.py <root_dir> <output_path> <test_function_name>")
        sys.exit(1)

    ROOT_DIR = Path(sys.argv[1]).resolve()
    OUTPUT_PATH = Path(sys.argv[2]).resolve()
    TEST_NAME = sys.argv[3]

    profile_test(TEST_NAME, ROOT_DIR, OUTPUT_PATH)
