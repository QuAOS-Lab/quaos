from types import SimpleNamespace
from sympleq.applications.randomized_benchmarking.experiments import cost_aware_surface_method as casm


def run():
    meta = {
        "settings": {
            "q_values": [27, 28],
            "backend_model": "sympleq",
            "rng_seed": 2025,
            "max_qubit_window": 5,
            "min_distinct_q_coverage": 30,
            "acquisition_q_resolution": 30,
            "acquisition_ratio_points": 15,
            "grid_resolution": [11, 11, 7, 7, 1, 5, 1],
            "boundary_fit_resolution": [17, 17, 9, 9, 1, 7, 1],
            "max_cost_per_run": 35.0,
        }
    }

    settings = SimpleNamespace(
        q_values=(27, 28),
        backend_model="sympleq",
        rng_seed=2025,
        max_qubit_window=5,
        min_distinct_q_coverage=30,
        acquisition_q_resolution=30,
        acquisition_ratio_points=15,
        grid_resolution=(11, 11, 7, 7, 1, 5, 1),
        boundary_fit_resolution=(17, 17, 9, 9, 1, 7, 1),
        max_cost_per_run=35.0,
    )

    print("Expect True:", casm._checkpoint_q_values_match(meta, settings))

    # Change max_cost_per_run should fail
    settings.max_cost_per_run = 40.0
    print("Expect False (max_cost_per_run changed):", casm._checkpoint_q_values_match(meta, settings))

    settings.max_cost_per_run = 35.0
    # Change min_distinct_q_coverage should fail
    settings.min_distinct_q_coverage = 31
    print("Expect False (min_distinct_q_coverage changed):", casm._checkpoint_q_values_match(meta, settings))


if __name__ == "__main__":
    run()
