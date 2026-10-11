#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      Complex Modeling Initiative
#                   (c) 2016 Brookhaven Science Associates,
#                   Brookhaven National Laboratory.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################
"""Unit tests for the weakrefcallable module."""

import pickle

import pytest

from diffpy.srfit import FitContribution
from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.util.weakrefcallable import WeakBoundMethod, weak_ref

# These tests build their objects inside each test instead of in fixtures,
# because pytest keeps fixture values alive and the tests need the object
# behind a weak reference to be garbage collected.


def _fallback_example(wbm, *args, **kwargs):
    return (wbm, args, kwargs)


def _make_contribution():
    contribution = FitContribution("f")
    contribution.set_equation("7")
    return contribution


# weak_ref stores the fallback it is given, or None without one.
@pytest.mark.parametrize(
    "input_fallback, expected_fallback",
    [
        # C1: A fallback function is given.
        # Expected: The wrapper stores it.
        (_fallback_example, _fallback_example),
        # C2: No fallback is given.
        # Expected: The wrapper's fallback is None.
        (None, None),
    ],
)
def test_weak_ref_fallback(input_fallback, expected_fallback):
    contribution = _make_contribution()
    weak_method = weak_ref(contribution._flush, fallback=input_fallback)
    actual_fallback = weak_method.fallback
    assert actual_fallback is expected_fallback


def test_call_live_object():
    # C1: The wrapped method's object is alive and the wrapper is called.
    # Expected: The bound method runs; here it flushes the cached value.
    contribution = _make_contribution()
    weak_flush = weak_ref(contribution._eq._flush)
    actual_evaluated_value = contribution.evaluate()
    actual_cached_value_before_flush = contribution._eq._value
    expected_evaluated_value = 7
    assert actual_evaluated_value == expected_evaluated_value
    assert actual_cached_value_before_flush == expected_evaluated_value
    weak_flush(())
    actual_cached_value_after_flush = contribution._eq._value
    expected_cached_value_after_flush = None
    assert actual_cached_value_after_flush == expected_cached_value_after_flush


def test_call_returns_method_result():
    # C1: The wrapped method returns a value and its object is alive.
    # Expected: The wrapper returns that value.
    parameter = Parameter("x", value=3)
    weak_get_value = weak_ref(parameter.get_value)
    actual_value = weak_get_value()
    expected_value = 3
    assert actual_value == expected_value


def test_call_dead_object_with_fallback():
    # C1: The wrapped method's object is deleted, then the wrapper is called.
    # Expected: The fallback is called with the wrapper and the arguments.
    contribution = _make_contribution()
    weak_flush = weak_ref(contribution._eq._flush, fallback=_fallback_example)
    del contribution
    assert weak_flush._wref() is None
    actual_call = weak_flush("any", "argument", foo=37)
    expected_call = (weak_flush, ("any", "argument"), {"foo": 37})
    assert actual_call[0] is expected_call[0]
    assert actual_call[1:] == expected_call[1:]


def test_call_dead_object_without_fallback_bad():
    # C1: The wrapped method's object is deleted and there is no fallback.
    # Expected: Calling the wrapper raises ReferenceError.
    parameter = Parameter("x", value=3)
    weak_get_value = weak_ref(parameter.get_value)
    del parameter
    with pytest.raises(ReferenceError):
        weak_get_value()


def test_hash():
    # C1: A wrapper's object is deleted, and the wrapper is pickled twice.
    # Expected: The hash is unchanged by the deletion, and both pickled
    # copies hash equally.
    contribution = FitContribution("f1")
    weak_flush = weak_ref(contribution._flush)
    expected_hash = hash(weak_flush)
    del contribution
    actual_wrapped_object = weak_flush._wref()
    expected_wrapped_object = None
    assert actual_wrapped_object == expected_wrapped_object
    actual_hash = hash(weak_flush)
    assert actual_hash == expected_hash
    first_copy = pickle.loads(pickle.dumps(weak_flush))
    second_copy = pickle.loads(pickle.dumps(weak_flush))
    actual_copy_hash = hash(second_copy)
    expected_copy_hash = hash(first_copy)
    assert actual_copy_hash == expected_copy_hash


def test_eq():
    # Two wrappers are equal when they wrap the same method of the same
    # object. A pickled copy holds an empty reference, so it only equals the
    # original once the original's object has died.
    contribution = FitContribution("f1")
    first = weak_ref(contribution._flush)
    second = weak_ref(contribution._flush)
    # C1: Two wrappers of the same live method.
    # Expected: They are equal.
    actual_is_equal = first == second
    expected_is_equal = True
    assert actual_is_equal == expected_is_equal
    # C2: A pickled copy of a wrapper whose object is alive.
    # Expected: The copy has an empty reference and is not equal.
    copy = pickle.loads(pickle.dumps(first))
    actual_wrapped_object = copy._wref()
    expected_wrapped_object = None
    assert actual_wrapped_object == expected_wrapped_object
    actual_is_equal = first == copy
    expected_is_equal = False
    assert actual_is_equal == expected_is_equal
    # C3: The original's object is deleted.
    # Expected: The original now equals the pickled copy.
    del contribution
    actual_wrapped_object = first._wref()
    expected_wrapped_object = None
    assert actual_wrapped_object == expected_wrapped_object
    actual_is_equal = first == copy
    expected_is_equal = True
    assert actual_is_equal == expected_is_equal
    # C4: A pickled copy of the pickled copy.
    # Expected: It has an empty reference and equals both.
    copy_of_copy = pickle.loads(pickle.dumps(copy))
    actual_wrapped_object = copy_of_copy._wref()
    expected_wrapped_object = None
    assert actual_wrapped_object == expected_wrapped_object
    actual_equals_copy = copy_of_copy == copy
    actual_equals_first = copy_of_copy == first
    expected_is_equal = True
    assert actual_equals_copy == expected_is_equal
    assert actual_equals_first == expected_is_equal


def test_pickling():
    # C1: A set holding a wrapper, the wrapped object and the wrapper are
    # pickled together. Unpickling calls the wrapper's __hash__.
    # Expected: The unpickled set contains the unpickled wrapper, which
    # refers to the unpickled object.
    contribution = _make_contribution()
    weak_flush = weak_ref(contribution._eq._flush, fallback=_fallback_example)
    holder = {weak_flush}
    data = pickle.dumps([holder, contribution._eq, weak_flush])
    actual_holder, actual_equation, actual_wrapper = pickle.loads(data)
    assert actual_wrapper in actual_holder
    assert actual_equation is actual_wrapper._wref()


def test_observable_drops_dead_observer():
    # C1: An observer's object is deleted, then the observable notifies.
    # Expected: The dead observer is dropped from the observable.
    contribution = _make_contribution()
    x = contribution.new_parameter("x", 5)
    contribution.set_equation("3 * x")
    assert contribution.evaluate() == 15
    observer = next(iter(x._observers))
    assert isinstance(observer, WeakBoundMethod)
    # Changing x resets the contribution's cached value.
    x.set_value(x.value + 1)
    assert contribution._eq._value is None
    assert contribution.evaluate() == 18
    del contribution
    assert observer in x._observers
    x.set_value(x.value + 1)
    actual_observer_count = len(x._observers)
    expected_observer_count = 0
    assert actual_observer_count == expected_observer_count
