#!/usr/bin/env python
"""Standard messages that are used in notifications."""

# Messages to be passed
# A parameter value has changed. This may trigger recomputation of other
# values.
VALUE_CHANGED = 1 << 0
# Parameter gets fixed or varied. Message must be passed to the residual so it
# can recollect variables.
VARY_CHANGED = 1 << 1
