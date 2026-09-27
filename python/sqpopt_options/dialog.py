"""The SQPOPT options dialog (Qt, via qtpy)."""

from __future__ import annotations

import json
from typing import Any

from qtpy import QtCore, QtGui, QtWidgets

from . import schema
from .schema import Option, format_number

__all__ = ['SqpoptOptionsDialog', 'edit_options']

_ERROR_STYLE = 'border: 1px solid #d0342c; border-radius: 3px;'
_EDITOR_WIDTH = 160  #: width of the number editors (a selector's combo box may be wider)


def _range_hint(o: Option) -> str:
    """e.g. '(0, 1)', '> 0', '>= 1', '0 = automatic'"""
    parts = []
    if o.minimum is not None and o.maximum is not None:
        lo = '(' if o.min_open else '['
        hi = ')' if o.max_open else ']'
        parts.append(f'{lo}{format_number(o.minimum)}, {format_number(o.maximum)}{hi}')
    elif o.minimum is not None:
        parts.append(f"{'>' if o.min_open else '>='} {format_number(o.minimum)}")
    elif o.maximum is not None:
        parts.append(f"{'<' if o.max_open else '<='} {format_number(o.maximum)}")
    for v, meaning in o.special.items():
        parts.append(f'{format_number(v)} = {meaning}')
    return ';  '.join(parts)


def _describe_value(o: Option, v: Any) -> str:
    if o.kind == 'bool':
        return '.true.' if v else '.false.'
    if o.kind == 'choice':
        c = next((c for c in o.choices if c.value == v), None)
        return (c.name or f'{c.value} ({c.label})') if c else str(v)
    if o.kind == 'float' and v == -schema.SQPOPT_INFINITY:
        return '-sqpopt_infinity (-1e20)'
    return format_number(v)


class _Editor(QtCore.QObject):
    """the widgets of one option's row: label, editor, hint, and reset button"""

    changed = QtCore.Signal()
    focused = QtCore.Signal(object)

    def __init__(self, o: Option, label: str, parent: QtWidgets.QWidget):
        super().__init__(parent)
        self.option = o
        self.label = QtWidgets.QLabel(label, parent)
        mono = QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.SystemFont.FixedFont)
        self.label.setFont(mono)
        self.label.setToolTip(o.fortran)
        self.hint = QtWidgets.QLabel(_range_hint(o), parent)
        self.hint.setEnabled(False)   # (drawn in the disabled text color: a secondary note)
        self.reset = QtWidgets.QToolButton(parent)
        self.reset.setText('↺')
        self.reset.setAutoRaise(True)
        self.reset.setToolTip(f'Reset to the default ({_describe_value(o, o.default)})')
        self.reset.clicked.connect(lambda: self.set_value(o.default))
        self._valid = True

        if o.kind == 'bool':
            w = QtWidgets.QCheckBox(parent)
            w.toggled.connect(self._on_change)
        elif o.kind == 'choice':
            w = QtWidgets.QComboBox(parent)
            for c in o.choices:
                w.addItem(c.label + (f'   ({c.name})' if c.name else f'   ({c.value})'), c.value)
            # (a long entry is elided in the box, but shown in full in the popup:)
            w.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
            w.setMinimumContentsLength(34)
            w.view().setMinimumWidth(w.view().sizeHintForColumn(0) + 24)
            w.currentIndexChanged.connect(self._on_change)
        elif o.kind == 'int':
            w = QtWidgets.QSpinBox(parent)
            lo = int(o.minimum) if o.minimum is not None else (0 if 0 in o.special else -2**31)
            w.setRange(lo, 2**31 - 1)
            if lo in o.special:
                w.setSpecialValueText(f'{o.special[lo]} ({lo})')
            w.setFixedWidth(_EDITOR_WIDTH)
            w.valueChanged.connect(self._on_change)
        else:
            w = QtWidgets.QLineEdit(parent)
            w.setFixedWidth(_EDITOR_WIDTH)
            w.textChanged.connect(self._on_change)
        w.setToolTip(o.doc)
        w.installEventFilter(self)
        self.label.installEventFilter(self)
        self.widget = w
        self.set_value(o.default)

    # ---- value

    def value(self) -> Any:
        """the current value (a float option with unparseable text returns the text)"""
        w, o = self.widget, self.option
        if o.kind == 'bool':
            return w.isChecked()
        if o.kind == 'choice':
            return w.currentData()
        if o.kind == 'int':
            return w.value()
        text = w.text().strip()
        try:
            return float(text)
        except ValueError:
            return text

    def set_value(self, v: Any) -> None:
        w, o = self.widget, self.option
        if o.kind == 'bool':
            w.setChecked(bool(v))
        elif o.kind == 'choice':
            i = w.findData(v)
            if i >= 0:
                w.setCurrentIndex(i)
        elif o.kind == 'int':
            w.setValue(int(v))
        else:
            w.setText(format_number(float(v)))
        self._update_state()

    def is_changed(self) -> bool:
        return self.value() != self.option.default

    def error(self) -> str | None:
        v = self.value()
        if isinstance(v, str):
            return 'not a number'
        return self.option.check(v)

    # ---- display

    def set_extra_error(self, msg: str | None) -> None:
        """show a cross-option problem (e.g. radius0 < radius_min) on this row"""
        self._show_error(msg or self.error())

    def _show_error(self, msg: str | None) -> None:
        self._valid = msg is None
        self.widget.setStyleSheet('' if self._valid else _ERROR_STYLE)
        tip = self.option.doc if self._valid else f'<b>{msg}</b><br>{self.option.doc}'
        self.widget.setToolTip(tip)

    def _update_state(self) -> None:
        changed = self.is_changed()
        f = self.label.font()
        f.setBold(changed)
        self.label.setFont(f)
        self.reset.setEnabled(changed)
        self._show_error(self.error())

    def _on_change(self, *_):
        self._update_state()
        self.changed.emit()

    def rows_widgets(self) -> tuple[QtWidgets.QWidget, ...]:
        if self.option.kind == 'choice':
            return (self.label, self.widget, self.reset)   # (its hint is not shown)
        return (self.label, self.widget, self.hint, self.reset)

    def set_visible(self, visible: bool) -> None:
        for w in self.rows_widgets():
            w.setVisible(visible)

    def matches(self, text: str) -> bool:
        t = text.lower()
        return not t or t in self.option.fortran.lower() or t in self.option.doc.lower()

    def eventFilter(self, obj, event):
        if event.type() in (QtCore.QEvent.Type.FocusIn, QtCore.QEvent.Type.Enter):
            self.focused.emit(self)
        return False


class _SectionBox(QtWidgets.QGroupBox):

    def __init__(self, section: schema.Section, label_width: int, parent: QtWidgets.QWidget):
        # the Fortran component the options are set on, in the title if they share one (other than
        # `options`); otherwise, the options not on the section's main component get qualified labels:
        prefixes = [o.path[:-1] for o in section.options]
        main = max(set(prefixes), key=prefixes.count)
        title = section.title
        if len(set(prefixes)) == 1 and main != ('options',):
            title += f"   ({'%'.join(main)}%…)"
        super().__init__(title, parent)
        self.section = section
        lay = QtWidgets.QGridLayout(self)
        lay.setColumnStretch(2, 1)
        lay.setColumnMinimumWidth(0, label_width)   # (the same in every section, so the editors line up)
        lay.setColumnMinimumWidth(1, _EDITOR_WIDTH)
        lay.setHorizontalSpacing(12)
        row = 0
        self.notice = QtWidgets.QLabel(self)
        self.notice.setWordWrap(True)
        self.notice.setStyleSheet('color: #b26a00; font-style: italic;')
        self.notice.hide()
        lay.addWidget(self.notice, row, 0, 1, 4)
        row += 1
        if section.note:
            note = QtWidgets.QLabel(section.note, self)
            note.setWordWrap(True)
            note.setEnabled(False)
            lay.addWidget(note, row, 0, 1, 4)
            row += 1
        self.editors: list[_Editor] = []
        for o in section.options:
            e = _Editor(o, o.name if o.path[:-1] == main else o.fortran, self)
            lay.addWidget(e.label, row, 0)
            if o.kind == 'choice':
                # (a selector has no hint: its combo box may use that column too)
                lay.addWidget(e.widget, row, 1, 1, 2, QtCore.Qt.AlignmentFlag.AlignLeft)
                e.hint.hide()
            else:
                lay.addWidget(e.widget, row, 1)
                lay.addWidget(e.hint, row, 2)
            lay.addWidget(e.reset, row, 3)
            self.editors.append(e)
            row += 1

    def update_relevance(self, values: dict) -> None:
        reason = self.section.relevance(values) if self.section.relevance else None
        self.notice.setText(f'These options are {reason}.' if reason else '')
        self.notice.setVisible(bool(reason))


class _TopicPage(QtWidgets.QScrollArea):

    def __init__(self, topic: schema.Topic, label_width: int, parent: QtWidgets.QWidget):
        super().__init__(parent)
        self.topic = topic
        self.setWidgetResizable(True)
        self.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        inner = QtWidgets.QWidget(self)
        lay = QtWidgets.QVBoxLayout(inner)
        title = QtWidgets.QLabel(topic.title, inner)
        f = title.font()
        f.setPointSizeF(f.pointSizeF()*1.3)
        f.setBold(True)
        title.setFont(f)
        lay.addWidget(title)
        summary = QtWidgets.QLabel(topic.summary, inner)
        summary.setWordWrap(True)
        summary.setEnabled(False)
        lay.addWidget(summary)
        self.boxes = [_SectionBox(s, label_width, inner) for s in topic.sections]
        for b in self.boxes:
            lay.addWidget(b)
        lay.addStretch(1)
        self.setWidget(inner)

    def editors(self) -> list[_Editor]:
        return [e for b in self.boxes for e in b.editors]


class SqpoptOptionsDialog(QtWidgets.QDialog):
    """a dialog for editing all the SQPOPT settings (not the problem definition, callbacks, or bounds).

    `values` (optional) are the initial settings: a nested dict as returned by `values()`, possibly with only
    some options, or a flat dict with Fortran references as keys (see `schema.merge_values`); the other
    options start at their defaults. `values()` returns the edited settings."""

    #: emitted whenever any option is changed
    valuesChanged = QtCore.Signal()

    def __init__(self, parent: QtWidgets.QWidget | None = None, values: dict | None = None,
                 title: str = 'SQPOPT Options'):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(980, 720)

        # ---- search bar
        self.search = QtWidgets.QLineEdit(self)
        self.search.setPlaceholderText('Search options (name or description)…')
        self.search.setClearButtonEnabled(True)
        self.changed_only = QtWidgets.QCheckBox('Changed only', self)
        top = QtWidgets.QHBoxLayout()
        top.addWidget(self.search, 1)
        top.addWidget(self.changed_only)

        # ---- topics and pages
        self.topic_list = QtWidgets.QListWidget(self)
        self.topic_list.setMaximumWidth(210)
        self.pages = QtWidgets.QStackedWidget(self)
        self._pages: list[_TopicPage] = []
        mono = QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.SystemFont.FixedFont)
        mono.setBold(True)
        label_width = max(QtGui.QFontMetrics(mono).horizontalAdvance(o.name) for o in schema.OPTIONS.values()) + 8
        # (a qualified label, e.g. hessian%damping, may be wider; its row then just widens the column)
        for t in schema.TOPICS:
            page = _TopicPage(t, label_width, self.pages)
            self.pages.addWidget(page)
            self._pages.append(page)
            item = QtWidgets.QListWidgetItem(t.title, self.topic_list)
            item.setToolTip(t.summary)
        self.topic_list.currentRowChanged.connect(self.pages.setCurrentIndex)
        self.topic_list.currentRowChanged.connect(lambda _: self._show_description(None))
        split = QtWidgets.QSplitter(self)
        split.addWidget(self.topic_list)
        split.addWidget(self.pages)
        split.setStretchFactor(1, 1)
        split.setChildrenCollapsible(False)

        self._editors: dict[str, _Editor] = {}
        for page in self._pages:
            for e in page.editors():
                self._editors[e.option.fortran] = e
                e.changed.connect(self._on_changed)
                e.focused.connect(self._show_description)

        # ---- description and status
        self.description = QtWidgets.QLabel(self)
        self.description.setWordWrap(True)
        self.description.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.description.setAlignment(QtCore.Qt.AlignmentFlag.AlignTop | QtCore.Qt.AlignmentFlag.AlignLeft)
        self.description.setMinimumHeight(84)
        self.description.setFrameShape(QtWidgets.QFrame.Shape.StyledPanel)
        self.description.setMargin(8)
        self.status = QtWidgets.QLabel(self)
        self.status.setStyleSheet('color: #d0342c;')
        self.status.setWordWrap(True)

        # ---- buttons
        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok | QtWidgets.QDialogButtonBox.StandardButton.Cancel |
            QtWidgets.QDialogButtonBox.StandardButton.RestoreDefaults, self)
        self._ok = buttons.button(QtWidgets.QDialogButtonBox.StandardButton.Ok)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        buttons.button(QtWidgets.QDialogButtonBox.StandardButton.RestoreDefaults).clicked.connect(
            self.restore_defaults)
        load = buttons.addButton('Load…', QtWidgets.QDialogButtonBox.ButtonRole.ActionRole)
        load.clicked.connect(self._load)
        save = buttons.addButton('Save…', QtWidgets.QDialogButtonBox.ButtonRole.ActionRole)
        save.clicked.connect(self._save)
        copy = buttons.addButton('Copy', QtWidgets.QDialogButtonBox.ButtonRole.ActionRole)
        menu = QtWidgets.QMenu(copy)
        menu.addAction('As a Python dict', lambda: self._copy(schema.python_literal(self.values())))
        menu.addAction('As a Python dict (changed options only)',
                       lambda: self._copy(schema.python_literal(self.values(changed_only=True))))
        menu.addAction('As Fortran assignments (changed options only)',
                       lambda: self._copy(schema.fortran_assignments(self.values())))
        menu.addAction('As JSON', lambda: self._copy(json.dumps(self.values(), indent=2)))
        copy.setMenu(menu)

        lay = QtWidgets.QVBoxLayout(self)
        lay.addLayout(top)
        lay.addWidget(split, 1)
        lay.addWidget(self.description)
        lay.addWidget(self.status)
        lay.addWidget(buttons)

        self.search.textChanged.connect(self._apply_filter)
        self.changed_only.toggled.connect(self._apply_filter)

        self.set_values(values)
        self.topic_list.setCurrentRow(0)
        self._show_description(None)

    # ---- the settings

    def values(self, changed_only: bool = False, enum_names: bool = False) -> dict:
        """the settings as a nested dict with the Fortran names (see `sqpopt_options.schema`).

        `changed_only`: only the options that differ from their defaults. `enum_names`: selector options as the
        names of their Fortran constants (e.g. ``'sqpopt_qp_auto'``) instead of their integer values. (Invalid
        entries are returned as typed; `validate()` lists them.)"""
        v = schema.unflatten({ref: e.value() for ref, e in self._editors.items()})
        if changed_only:
            v = schema.changed_values(v)
        if enum_names:
            v = schema.to_enum_names(v)
        return v

    def set_values(self, values: dict | None) -> None:
        """set the options present in `values` (nested or flat, see `schema.merge_values`); the others are
        reset to their defaults"""
        full = schema.merge_values(values)
        for ref, e in self._editors.items():
            e.blockSignals(True)
            e.set_value(schema.get_value(full, e.option.path))
            e.blockSignals(False)
        self._on_changed()

    def restore_defaults(self) -> None:
        self.set_values(None)

    def validate(self) -> dict[str, str]:
        """the invalid options, as ``{fortran reference: message}``"""
        errors = {ref: msg for ref, e in self._editors.items() if (msg := e.error())}
        if not errors:
            errors = schema.validate(self.values())
        return errors

    # ---- internals

    def _on_changed(self) -> None:
        errors = self.validate()
        for ref, e in self._editors.items():
            e.set_extra_error(errors.get(ref))
        values = self.values()
        for page in self._pages:
            for b in page.boxes:
                b.update_relevance(values)
        for i, page in enumerate(self._pages):
            item = self.topic_list.item(i)
            n_changed = sum(e.is_changed() for e in page.editors())
            n_errors = sum(e.option.fortran in errors for e in page.editors())
            item.setText(page.topic.title + (f'  ({n_changed})' if n_changed else ''))
            f = item.font()
            f.setBold(n_changed > 0)
            item.setFont(f)
            item.setForeground(QtGui.QBrush(QtGui.QColor('#d0342c')) if n_errors else QtGui.QBrush())
        if errors:
            ref, msg = next(iter(errors.items()))
            more = f'  (and {len(errors) - 1} more)' if len(errors) > 1 else ''
            self.status.setText(f'{ref} {msg}{more}')
        else:
            self.status.setText('')
        self._ok.setEnabled(not errors)
        if self.changed_only.isChecked():
            self._apply_filter()
        self.valuesChanged.emit()

    def _apply_filter(self) -> None:
        text = self.search.text().strip()
        only_changed = self.changed_only.isChecked()
        first_visible = None
        for i, page in enumerate(self._pages):
            any_page = False
            for b in page.boxes:
                any_box = False
                for e in b.editors:
                    visible = e.matches(text) and (e.is_changed() or not only_changed)
                    e.set_visible(visible)
                    any_box |= visible
                b.setVisible(any_box)
                any_page |= any_box
            self.topic_list.item(i).setHidden(not any_page)
            if any_page and first_visible is None:
                first_visible = i
        current = self.topic_list.currentRow()
        if first_visible is not None and (current < 0 or self.topic_list.item(current).isHidden()):
            self.topic_list.setCurrentRow(first_visible)

    def _show_description(self, editor: _Editor | None) -> None:
        if editor is None:
            self.description.setText('Hover over or select an option to see its description. '
                                     'Changed options are shown in bold; ↺ resets one to its default.')
            return
        o = editor.option
        kind = {'int': 'integer', 'float': 'real', 'bool': 'logical', 'choice': 'integer selector'}[o.kind]
        hint = _range_hint(o)
        self.description.setText(
            f'<b><code>{o.fortran}</code></b> &nbsp; <i>{kind}, default {_describe_value(o, o.default)}'
            + (f', {hint}' if hint else '') + f'</i><br>{o.doc}')

    def _copy(self, text: str) -> None:
        QtWidgets.QApplication.clipboard().setText(text)

    def _load(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Load SQPOPT options', '', 'JSON (*.json)')
        if not path:
            return
        try:
            with open(path) as f:
                self.set_values(json.load(f))
        except (OSError, ValueError, KeyError) as exc:
            QtWidgets.QMessageBox.warning(self, 'Load SQPOPT options', f'Could not load {path}:\n{exc}')

    def _save(self) -> None:
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save SQPOPT options', 'sqpopt_options.json',
                                                        'JSON (*.json)')
        if path:
            with open(path, 'w') as f:
                json.dump(self.values(), f, indent=2)


def edit_options(values: dict | None = None, parent: QtWidgets.QWidget | None = None,
                 changed_only: bool = False, enum_names: bool = False, **kwargs) -> dict | None:
    """show the dialog modally, starting from `values`, and return the edited settings (see
    `SqpoptOptionsDialog.values` for `changed_only` and `enum_names`), or ``None`` if it was cancelled.
    Creates the ``QApplication`` if there is none yet."""
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])  # noqa: F841 (kept alive)
    dlg = SqpoptOptionsDialog(parent, values=values, **kwargs)
    if dlg.exec():
        return dlg.values(changed_only=changed_only, enum_names=enum_names)
    return None
