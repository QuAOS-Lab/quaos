# H-Wrapper FLE Variant

This folder is the independent 2D FLE variant with a protected native-H wrapper.

Circuit structure:

```text
S V S on every qubit
single random 1q scrambler layer
random Clifford/native body
inverse body
inverse scrambler
S V S on every qubit
measure
```

The H wrapper is protected from deterministic identity elimination. The nominal
config gate counts include the wrapper gates so the Quantinuum backend's gate
count checks remain valid.

Run:

```powershell
C:\tmp\simpleq312\Scripts\python.exe src\sympleq\applications\randomized_benchmarking\experiments\GP_Levelset_estimation\H_wrapper\run_FLE.py
```

Plot one checkpoint:

```powershell
C:\tmp\simpleq312\Scripts\python.exe src\sympleq\applications\randomized_benchmarking\experiments\GP_Levelset_estimation\H_wrapper\plot_score_fle.py PATH\TO\measurement_XXX_globalsur_YYYY.json
```

Output folders go under:

```text
Personal\FLE\H_wrapper\...
```
