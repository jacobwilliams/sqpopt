"""Tests of sqpopt_options.schema (no Qt needed).

The schema is checked against the Fortran sources: every option must be a
field of the Fortran type it is set on, with the same default, and every
field of those types must be either an option or listed here as internal
state (so that an option added to the Fortran code is not forgotten here).

Run from the repository root:  pixi run python -m unittest discover -s python/tests
"""

import pathlib
import re
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from sqpopt_options import schema  # noqa: E402

SRC = pathlib.Path(__file__).resolve().parents[2] / 'src'

#: the Fortran type each component path is, and its fields that are not user options
TYPES = {
    ('options',): ('sqpopt_options_module.f90', 'sqpopt_options_type', set()),
    ('hessian',): ('sqpopt_hessian_module.f90', 'sqpopt_hessian_type',
                   {'n', 'max_history', 'n_history', 'use_sr1', 'first', 'gamma', 'gamma0', 'exact', 'shift', 'shift_dominant',
                    'mid_valid', 'mid_ok', 'diag_valid'}),
    ('qp_solver',): ('sqpopt_qp_solver_module.f90', 'sqpopt_qp_solver_type',
                     {'mode', 'step_scale', 'capped', 'n_short', 'n_elastic', 'n_iter', 'negative_curvature', 'n_working',
                      'n_slacks', 'time', 'out_of_memory', 'direct', 'direct_used', 'direct_outcome', 'direct_changes',
                      'n_solves', 'n_direct', 'unconstrained_used', 'n_unconstrained'}),
    ('qp_solver', 'dense_qp'): ('sqpopt_qp_dense_module.f90', 'sqpopt_dense_qp_type',
                                {'n_iter', 'negative_curvature', 'force_weight', 'n_working', 'n_slacks'}),
    ('qp_solver', 'sparse_qp'): ('sqpopt_qp_reduced_hessian_module.f90', 'sqpopt_reduced_hessian_qp_type',
                                 {'n_iter', 'negative_curvature', 'force_weight', 'n_working', 'n_slacks'}),
    ('linesearch',): ('sqpopt_linesearch_module.f90', 'sqpopt_linesearch_type',
                      {'mode', 'nm_count', 'watchdog_ready', 'watchdog_relaxed_remaining',
                       'watchdog_cooldown_remaining', 'watchdog_w_opt', 'used_soc', 'used_nonmonotone',
                       'used_relaxed'}),
    ('linesearch', 'merit'): ('sqpopt_merit_module.f90', 'sqpopt_merit_type',
                              {'mode', 'penalty_update', 'joint_active', 'penalty_floor'}),
    ('linesearch', 'filter'): ('sqpopt_filter_module.f90', 'sqpopt_filter_type',
                               {'ready', 'theta_max', 'theta_min'}),
    ('linesearch', 'funnel'): ('sqpopt_funnel_module.f90', 'sqpopt_funnel_type', {'ready', 'width'}),
    ('trust_region',): ('sqpopt_trust_region_module.f90', 'sqpopt_trust_region_type', {'ready', 'radius', 'used_soc'}),
}

FIELD = re.compile(r'^\s*(integer|real\(wp\)|logical)\s*::\s*(\w+)\s*=\s*([^!]+?)\s*(!.*)?$', re.IGNORECASE)


#: the Fortran type of each option kind
FORTRAN_TYPE = {'int': 'integer', 'choice': 'integer', 'float': 'real(wp)', 'bool': 'logical'}


def fortran_fields(filename: str, type_name: str) -> dict[str, tuple[str, str]]:
    """the scalar fields with a default value of a Fortran derived type: {name: (type, default as written)}"""
    text = (SRC / filename).read_text()
    m = re.search(rf'type\s*,\s*public\s*::\s*{type_name}\b(.*?)^\s*(contains|end type)', text,
                  re.IGNORECASE | re.DOTALL | re.MULTILINE)
    assert m, f'{type_name} not found in {filename}'
    fields = {}
    for line in m.group(1).splitlines():
        f = FIELD.match(line)
        if f:
            fields[f.group(2).lower()] = (f.group(1).lower(), f.group(3).strip())
    return fields


def fortran_default(text: str):
    """the Python value of a Fortran default as written in the source"""
    t = text.lower()
    if t == '.true.':
        return True
    if t == '.false.':
        return False
    if t == 'output_unit':
        return 6
    if t == '-sqpopt_infinity':
        return -schema.SQPOPT_INFINITY
    for choices in schema.ALL_CHOICES:
        for c in choices:
            if c.name == t:
                return c.value
    t = t.replace('_wp', '')
    return float(t) if any(ch in t for ch in '.e') else int(t)


class TestAgainstFortran(unittest.TestCase):

    def test_options_match_fortran(self):
        for prefix, (filename, type_name, internal) in TYPES.items():
            fields = fortran_fields(filename, type_name)
            options = {o.name: o for o in schema.OPTIONS.values() if o.path[:-1] == prefix}
            for name, o in options.items():
                with self.subTest(option=o.fortran):
                    self.assertIn(name, fields, f'{o.fortran} is not a field of {type_name}')
                    ftype, default = fields[name]
                    self.assertEqual(fortran_default(default), o.default, f'{o.fortran}: default')
                    # (the kind must match the Fortran type: the generated setter of the Python bindings
                    # converts by kind, and e.g. `100 == 100.0` would hide a mismatch in the default test)
                    self.assertEqual(ftype, FORTRAN_TYPE[o.kind], f'{o.fortran}: kind {o.kind!r} for a {ftype}')
                    self.assertIsInstance(o.default, {'int': int, 'float': float, 'bool': bool, 'choice': int}[o.kind])
            for name in fields:
                with self.subTest(field=f'{type_name}%{name}'):
                    self.assertTrue(name in options or name in internal,
                                    f'{type_name}%{name} is neither an option nor listed as internal state')

    def test_every_component_path_is_known(self):
        for o in schema.OPTIONS.values():
            self.assertIn(o.path[:-1], TYPES, o.fortran)


class TestValues(unittest.TestCase):

    def test_defaults_are_valid(self):
        self.assertEqual(schema.validate(schema.default_values()), {})

    def test_each_option_once(self):
        refs = [o.fortran for t in schema.TOPICS for o in t.options()]
        self.assertEqual(len(refs), len(set(refs)))

    def test_nested_layout(self):
        v = schema.default_values()
        self.assertEqual(v['qp_solver']['dense_qp']['max_iter'], 100)
        self.assertEqual(v['linesearch']['filter']['gamma_theta'], 1e-5)
        self.assertEqual(v['options']['qp_solver_mode'], 0)
        self.assertIs(v['trust_region']['enabled'], False)

    def test_flatten_roundtrip(self):
        v = schema.default_values()
        self.assertEqual(schema.unflatten(schema.flatten(v)), v)

    def test_merge(self):
        v = schema.merge_values({'options': {'max_iter': 500, 'linesearch_mode': 'sqpopt_linesearch_funnel'}})
        self.assertEqual(v['options']['max_iter'], 500)
        self.assertEqual(v['options']['linesearch_mode'], 5)
        self.assertEqual(v['options']['ktol'], 1e-6)
        v = schema.merge_values({'qp_solver%max_step': 5})
        self.assertEqual(v['qp_solver']['max_step'], 5.0)
        self.assertIsInstance(v['qp_solver']['max_step'], float)
        with self.assertRaises(KeyError):
            schema.merge_values({'options': {'no_such_option': 1}})

    def test_changed_and_enum_names(self):
        v = schema.merge_values({'options': {'qp_solver_mode': 3}, 'linesearch': {'funnel': {'beta': 0.5}}})
        self.assertEqual(schema.changed_values(v),
                         {'options': {'qp_solver_mode': 3}, 'linesearch': {'funnel': {'beta': 0.5}}})
        named = schema.to_enum_names(v)
        self.assertEqual(named['options']['qp_solver_mode'], 'sqpopt_qp_reduced_hessian')
        self.assertEqual(named['options']['print_level'], 0)  # (no Fortran constant)
        self.assertEqual(schema.to_fortran_values(named), v)

    def test_validation(self):
        v = schema.merge_values({'options': {'ktol': 0.0}})
        self.assertIn('options%ktol', schema.validate(v))
        v = schema.merge_values({'qp_solver': {'dense_qp': {'elastic_weight_max': 1.0}}})
        self.assertIn('qp_solver%dense_qp%elastic_weight_max', schema.validate(v))
        v = schema.merge_values({'trust_region': {'enabled': True, 'eta1': 0.9, 'eta2': 0.5}})
        self.assertIn('trust_region%eta2', schema.validate(v))
        v = schema.merge_values({'trust_region': {'eta1': 0.9, 'eta2': 0.5}})   # (not checked unless enabled)
        self.assertEqual(schema.validate(v), {})
        v = schema.default_values()
        v['options']['max_iter'] = 1.5
        self.assertIn('options%max_iter', schema.validate(v))

    def test_fortran_assignments(self):
        v = schema.merge_values({'options': {'linesearch_mode': 1, 'scaling': False, 'ktol': 1e-9},
                                 'qp_solver': {'max_step': 5.0}})
        lines = schema.fortran_assignments(v).splitlines()
        self.assertIn('options%linesearch_mode = sqpopt_linesearch_armijo', lines)
        self.assertIn('options%scaling = .false.', lines)
        self.assertIn('options%ktol = 1e-9_wp', lines)
        self.assertIn('qp_solver%max_step = 5.0_wp', lines)

    def test_python_literal(self):
        v = schema.default_values()
        self.assertEqual(eval(schema.python_literal(v)), v)


if __name__ == '__main__':
    unittest.main()
