#!/usr/bin/env python3
"""Exact reviewed public surface gate, with a useful diff for deliberate changes."""
import difflib
import sys
from pathlib import Path


def check(baseline, actual):
    expected = Path(baseline).read_text().splitlines()
    observed = Path(actual).read_text().splitlines()
    if not expected or expected != observed:
        print('\n'.join(difflib.unified_diff(expected, observed, fromfile=str(baseline), tofile=str(actual))))
        raise ValueError('Public API changed. Review the diff and explicitly update docs/public-api.txt.')


if __name__ == '__main__':
    try:
        if len(sys.argv) != 3:
            raise ValueError('Usage: check-api.py BASELINE ACTUAL')
        check(sys.argv[1], sys.argv[2])
    except (ValueError, OSError) as error:
        sys.exit(str(error))
