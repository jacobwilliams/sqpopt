"""Tests of the SQPOPT options dialog, run offscreen (no window is shown).

Run from the repository root:  pixi run python -m unittest discover -s python/tests
"""

import os
import pathlib
import sys
import unittest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from qtpy import QtWidgets  # noqa: E402

from sqpopt_options import schema, SqpoptOptionsDialog  # noqa: E402

app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


class TestDialog(unittest.TestCase):

    def setUp(self):
        self.dlg = SqpoptOptionsDialog()

    def tearDown(self):
        self.dlg.deleteLater()

    def editor(self, ref):
        return self.dlg._editors[ref]

    def test_defaults(self):
        self.assertEqual(self.dlg.values(), schema.default_values())
        self.assertEqual(self.dlg.values(changed_only=True), {})
        self.assertTrue(self.dlg._ok.isEnabled())
        self.assertEqual(set(self.dlg._editors), set(schema.OPTIONS))

    def test_edit_each_kind(self):
        self.editor('options%max_iter').widget.setValue(500)
        self.editor('options%ktol').widget.setText('1e-9')
        self.editor('options%scaling').widget.setChecked(False)
        combo = self.editor('options%linesearch_mode').widget
        combo.setCurrentIndex(combo.findData(1))
        self.assertEqual(self.dlg.values(changed_only=True),
                         {'options': {'max_iter': 500, 'ktol': 1e-9, 'scaling': False, 'linesearch_mode': 1}})
        self.assertEqual(self.dlg.values(enum_names=True)['options']['linesearch_mode'], 'sqpopt_linesearch_armijo')
        self.assertTrue(self.editor('options%max_iter').reset.isEnabled())
        self.editor('options%max_iter').reset.click()
        self.assertEqual(self.dlg.values()['options']['max_iter'], 100)

    def test_initial_values(self):
        dlg = SqpoptOptionsDialog(values={'qp_solver': {'sparse_qp': {'null_space': 'sqpopt_null_space_lsqr'}},
                                          'linesearch%funnel%beta': 0.5})
        self.assertEqual(dlg.values(changed_only=True),
                         {'qp_solver': {'sparse_qp': {'null_space': 2}}, 'linesearch': {'funnel': {'beta': 0.5}}})
        dlg.restore_defaults()
        self.assertEqual(dlg.values(changed_only=True), {})

    def test_invalid_disables_ok(self):
        e = self.editor('linesearch%sigma')
        e.widget.setText('1.5')
        self.assertIn('linesearch%sigma', self.dlg.validate())
        self.assertFalse(self.dlg._ok.isEnabled())
        self.assertTrue(e.widget.styleSheet())
        e.widget.setText('abc')
        self.assertEqual(self.dlg.validate()['linesearch%sigma'], 'not a number')
        e.widget.setText('0.2')
        self.assertTrue(self.dlg._ok.isEnabled())
        # a cross-option rule:
        self.editor('qp_solver%dense_qp%elastic_weight_max').widget.setText('1')
        self.assertIn('qp_solver%dense_qp%elastic_weight_max', self.dlg.validate())
        self.assertFalse(self.dlg._ok.isEnabled())

    def test_search_filter(self):
        self.dlg.search.setText('funnel%beta')
        visible = [ref for ref, e in self.dlg._editors.items() if not e.widget.isHidden()]
        self.assertEqual(visible, ['linesearch%funnel%beta'])
        visible_topics = [self.dlg.topic_list.item(i).text() for i in range(self.dlg.topic_list.count())
                          if not self.dlg.topic_list.item(i).isHidden()]
        self.assertEqual(visible_topics, ['Funnel'])
        self.dlg.search.setText('')
        self.dlg.changed_only.setChecked(True)
        self.assertTrue(all(e.widget.isHidden() for e in self.dlg._editors.values()))
        self.editor('options%max_iter').widget.setValue(7)
        self.assertFalse(self.editor('options%max_iter').widget.isHidden())

    def test_relevance_notices(self):
        filt = next(b for p in self.dlg._pages for b in p.boxes if b.section.title == 'Filter')
        merit = next(b for p in self.dlg._pages for b in p.boxes if b.section.title == 'Penalty parameter')
        self.assertTrue(filt.notice.isHidden())          # (filter is the default)
        self.assertFalse(merit.notice.isHidden())
        combo = self.editor('options%linesearch_mode').widget
        combo.setCurrentIndex(combo.findData(1))         # armijo
        self.assertFalse(filt.notice.isHidden())
        self.assertTrue(merit.notice.isHidden())

    def test_save_uses_integer_values(self):
        import json
        import tempfile
        from unittest import mock
        combo = self.editor('options%linesearch_mode').widget
        combo.setCurrentIndex(combo.findData(5))    # funnel
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / 'opts.json'
            with mock.patch.object(QtWidgets.QFileDialog, 'getSaveFileName', return_value=(str(path), '')):
                self.dlg._save()
            saved = json.loads(path.read_text())
            self.assertEqual(saved['options']['linesearch_mode'], 5)
            self.assertEqual(saved['qp_solver']['sparse_qp']['null_space'], 1)
            # (and loading it back, or a file with constant names, gives the same settings)
            dlg = SqpoptOptionsDialog(values=saved)
            self.assertEqual(dlg.values(), self.dlg.values())

    def test_copy_formats(self):
        self.editor('qp_solver%max_step').widget.setText('5')
        self.dlg._copy(schema.fortran_assignments(self.dlg.values()))
        self.assertEqual(QtWidgets.QApplication.clipboard().text(), 'qp_solver%max_step = 5.0_wp')


if __name__ == '__main__':
    unittest.main()
