#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Interfluo
#
# FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
# for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

"""Shared boilerplate for the examples in this directory.

An installed FlatLand needs none of this — `import flatland` just works. These
examples run straight out of a source checkout, so they put the repository's
`python/` directory on sys.path first and give a clear message if the shared
library has not been built yet.
"""

import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def load():
    """Import flatland from this checkout, or explain what to do about it."""
    sys.path.insert(0, os.path.join(REPO, "python"))
    try:
        import flatland
    except ImportError:
        sys.exit("cannot find the flatland package in %s/python" % REPO)
    except Exception as exc:                       # the library itself is missing
        sys.exit("%s\n\nBuild it first:\n    cd %s && make lib" % (exc, REPO))
    return flatland


def mesh_path(*parts):
    """Path to one of the meshes committed under examples/."""
    return os.path.join(REPO, "examples", *parts)


def rule(title):
    print("\n" + title)
    print("-" * len(title))
