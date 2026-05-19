# Place your `fieldMinMax.dat` files here

The four expected filenames are:

- `coarse_fieldMinMax.dat`
- `medium_fieldMinMax.dat`
- `fine_fieldMinMax.dat`
- `extra-fine_fieldMinMax.dat`

Then run:

```bash
python ../verify_abstract.py \
    --coarse     coarse_fieldMinMax.dat \
    --medium     medium_fieldMinMax.dat \
    --fine       fine_fieldMinMax.dat \
    --extra-fine extra-fine_fieldMinMax.dat
```

These files are not version-controlled by default (see `.gitignore`).
If you want to commit them — for archival reproducibility — remove
the `*.dat` entry from `.gitignore`.
