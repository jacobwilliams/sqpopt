"""Run the SQPOPT options dialog on its own, and print the result.

    python -m sqpopt_options [--load FILE.json] [--changed-only] [--enum-names] [--format python|json|fortran]
"""

import argparse
import json
import sys

from . import schema


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog='python -m sqpopt_options', description=__doc__.splitlines()[0])
    p.add_argument('--load', metavar='FILE', help='initial settings (JSON, nested or flat)')
    p.add_argument('--changed-only', action='store_true', help='print only the options changed from the defaults')
    p.add_argument('--enum-names', action='store_true',
                   help='give selector options as their Fortran constant names instead of integers')
    p.add_argument('--format', choices=('python', 'json', 'fortran'), default='python',
                   help='output format (default: a Python dict)')
    args = p.parse_args(argv)

    values = None
    if args.load:
        with open(args.load) as f:
            values = json.load(f)

    from .dialog import edit_options
    result = edit_options(values, changed_only=args.changed_only and args.format != 'fortran',
                          enum_names=args.enum_names)
    if result is None:
        return 1
    if args.format == 'json':
        print(json.dumps(result, indent=2))
    elif args.format == 'fortran':
        print(schema.fortran_assignments(schema.merge_values(result), only_changed=args.changed_only))
    else:
        print(schema.python_literal(result))
    return 0


if __name__ == '__main__':
    sys.exit(main())
