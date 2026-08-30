# RF models

This folder contains the source definition and post-processing for the custom
GSG travelling-wave electrodes.

- `hfss_cpw_config.json`: stack placeholders, parameter grid and pass/fail
  targets. Replace every unconfirmed foundry value before accepting a solve.
- `build_hfss_cpw.py`: PyAEDT builder for 0.5-mm, 1.0-mm and 10-mm HFSS Driven
  Terminal CPW models.
- `postprocess_tline.py`: two-length ABCD extraction of `Z0`, propagation loss
  and microwave effective index from unrenormalized Touchstone files.

Check the environment without launching AEDT:

```powershell
.\.venv\Scripts\python.exe models\rf\build_hfss_cpw.py check
.\.venv\Scripts\python.exe models\rf\postprocess_tline.py --self-test
```

Generated AEDT projects, meshes, Touchstone data and field exports belong under
`results/hfss/` and remain ignored by Git.
