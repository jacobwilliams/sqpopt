# SQPOPT options dialog

A Qt dialog (via [qtpy](https://github.com/spyder-ide/qtpy), so it works with
PySide6 or PyQt) for editing all of SQPOPT's settings, for inclusion in other
programs. The settings are the solver options, the Hessian, QP solver, line
search, merit function, filter, funnel, and trust-region settings. They don't
include the problem definition, the callbacks, or the bounds. The dialog groups
them by topic, validates them with the same rules as the Fortran code, and returns a
Python dict.

```python
from sqpopt_options import edit_options

settings = edit_options(parent=main_window)     # None if cancelled
settings = edit_options(values=settings)        # edit previous settings
```

or, to keep the dialog:

```python
from sqpopt_options import SqpoptOptionsDialog

dlg = SqpoptOptionsDialog(parent, values=settings)
if dlg.exec():
    settings = dlg.values()                  # all options
    changes  = dlg.values(changed_only=True) # only those changed from the defaults
```

## The settings dict

The dict mirrors the Fortran objects the options are set on, with the Fortran
names:

```python
settings['options']['max_iter']                  # options%max_iter
settings['hessian']['damping']                   # hessian%damping
settings['qp_solver']['dense_qp']['opt_tol']     # qp_solver%dense_qp%opt_tol
settings['linesearch']['filter']['gamma_theta']  # linesearch%filter%gamma_theta
settings['trust_region']['enabled']              # trust_region%enabled
```

Selector options (`hessian_mode`, `qp_solver_mode`, `linesearch_mode`,
`merit_mode`, `penalty_update`, `restoration_mode`, `null_space`) hold the
integer value of the Fortran constant (e.g. `sqpopt_qp_auto` = `0`). With
`values(enum_names=True)`, they hold the constant's name instead. Initial
values may use either form, and may be nested or flat with Fortran references
as keys (`{'qp_solver%max_step': 5.0}`). Options not given start at their
defaults.

The dialog can also save and load the settings as JSON, and copy them as a
Python dict, as JSON, or as Fortran assignment statements
(`qp_solver%max_step = 5.0_wp`). The dict and JSON give selector options as
their integer values (the menus show the names); loading accepts either
form.

## Without Qt

`sqpopt_options.schema` has the option definitions (Fortran names, types,
defaults, limits, and descriptions, grouped by topic) and the dict helpers,
and doesn't need Qt:

| function | |
|---|---|
| `default_values()` | all options at their defaults |
| `merge_values(values)` | a complete dict from partial, flat, or named values |
| `validate(values)` | `{fortran reference: problem}` (empty if valid) |
| `changed_values(values)` | only the options that differ from the defaults |
| `flatten(values)`, `unflatten(flat)` | to and from `{'qp_solver%max_step': 2.0, ...}` |
| `to_enum_names(values)`, `to_fortran_values(values)` | selector options as names or integers |
| `fortran_assignments(values)` | Fortran assignment statements |

## Running it

From this folder:

```sh
pixi run python -m sqpopt_options                        # show the dialog, print the result
pixi run python -m sqpopt_options --changed-only --format fortran
pixi run python -m unittest discover -s tests            # the tests
```

`tests/test_schema.py` checks the option definitions against the Fortran
sources in `../src`: every option must be a field of its Fortran type, with
the same default, and every field of those types must be an option or
listed as internal state. When an option is added to the Fortran code, add it
to `schema.py` too; the test fails until you do.
