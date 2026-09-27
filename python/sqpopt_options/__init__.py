"""An options dialog for SQPOPT, for inclusion in other programs.

Usage (Qt, via qtpy)::

    from sqpopt_options import edit_options

    settings = edit_options(parent=main_window)        # None if cancelled
    settings = edit_options(values=settings)           # edit previous settings

or, to keep the dialog around::

    from sqpopt_options import SqpoptOptionsDialog

    dlg = SqpoptOptionsDialog(parent, values=settings)
    if dlg.exec():
        settings = dlg.values()

The result is a dict mirroring the Fortran objects the options are set on,
with the Fortran names (see `sqpopt_options.schema`)::

    settings['options']['max_iter']                 # options%max_iter
    settings['qp_solver']['dense_qp']['opt_tol']    # qp_solver%dense_qp%opt_tol
    settings['linesearch']['filter']['gamma_theta'] # linesearch%filter%gamma_theta

The schema and the dict helpers (`default_values`, `merge_values`,
`validate`, `flatten`, ...) don't need Qt, and are importable on their own
from `sqpopt_options.schema`.
"""

from .schema import (TOPICS, OPTIONS, SQPOPT_INFINITY, default_values, merge_values, changed_values, validate,
                     flatten, unflatten, to_enum_names, to_fortran_values, python_literal, fortran_assignments)

__all__ = ['SqpoptOptionsDialog', 'edit_options', 'TOPICS', 'OPTIONS', 'SQPOPT_INFINITY', 'default_values',
           'merge_values', 'changed_values', 'validate', 'flatten', 'unflatten', 'to_enum_names',
           'to_fortran_values', 'python_literal', 'fortran_assignments']


def __getattr__(name):
    # the dialog is imported on first use, so that the schema can be used without Qt
    if name in ('SqpoptOptionsDialog', 'edit_options'):
        from . import dialog
        return getattr(dialog, name)
    raise AttributeError(name)
